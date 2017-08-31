#!/usr/bin/env python


#Epignetic pipeline for Illumina 450k and 850k

#Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys
import csv

#Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

#MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.readset import *
from bfx.design import * 

from bfx import rmarkdown
from bfx import qc
from bfx import dmp
from bfx import norm
from bfx import enrichment
from bfx import pca
from bfx import model
from pipelines import common
import utils

log = logging.getLogger(__name__)

class Epigene(common.MUGQICPipeline):
	"""
	Epigene 450k Pipeline
	=====================

	This pipeline is a modified implementation of the open source R pacakge minifi (http://bioconductor.org/packages/release/bioc/html/minfi.html)
	for the Illumina 450k platform. The pipeline takes raw .idat files and performs a series of optional preprocessing steps on it. Several 
	normalization methods are also available to be used before the data undergoes analysis for differentially methylated positions/regions.
	The output contains a list CpG sites and their methylation levels, a list of differentially methylated positions and a list of 
	differentially methylated regions, based on the contrasts selected in the design file. This step can be optionally run using the open source R package limma (http://bioconductor.org/packages/release/bioc/html/limma.html).
	
	The pipeline can then take the list of cpgs either given by the user or produced at the end of differential analysis and run enrichment analysis using the open source package LOLA (https://bioconductor.org/packages/release/bioc/html/LOLA.html).

	The final steps(optional) of the pipeline involve applying the classification models listed in the models design file on each sample and 	building a model for the sample data provided.
	
	At the end of the analysis, the pipeline can produce an html report outlining all the steps performed, parameters used and includes links to each data file. 
	"""

	RG_DIR = "report/data/rgchannel_set"
	GR_DIR = "report/data/gr_set"
	QC_DIR = "report/data/quality_control"
	CPG_DIR = "report/data/cpg"
	REP_DIR = "report"

	def __init__(self):
		self.argparser.add_argument("-r", "--samples", help="sample file", type=file)
		self.argparser.add_argument("-d", "--design", help="design file", type=file)
		super(Epigene, self).__init__()


	@property
	def samples(self):
		if not hasattr(self, "_samples"):
			if self.args.samples:
				samples = []
				illumina450k_sample_file=self.args.samples.name
				log.info("Parsing Illumina 450k sample sheet file " + illumina450k_sample_file + "...")
        			sample_csv = csv.DictReader(open(illumina450k_sample_file, "rb"), delimiter = ',')
        			for line in sample_csv:
                			sample = Illumina450kSample(line["Sample_Name"],
                                                                        os.path.join(os.path.dirname(illumina450k_sample_file), line["Sentrix_ID"], line["Sentrix_ID"] + "_" + line["Sentrix_Position"]))
                			sample._batch = line.get(config.param('column_names','batch'), None)
                			sample._cell_type = line.get("Tissue", None)
                			sample._sex = line.get(config.param('column_names','sex'), None)
                			sample._status = line.get(config.param('column_names','group'), None)
                			samples.append(sample)
				self._samples = samples
			else:
				self.argparser.error("arguement -r/--samples is required!")

		return self._samples 	

	#parse the contrasts supplied with the design file
	@property
	def contrasts(self):
		if self.args.design:
			return parse_design_file(self.args.design.name, self.samples)
		
		else:
			return []
			

	@property
	def steps(self):
		"""
		These are steps in the pipeline and they are listed in the same order they're processed.
		Any future steps added must also be mentioned here."""
		return [
                	self.read_idat,
			self.quality_control,
			self.normalize,
			self.probe_filter,
			self.cell_counts,
			self.batch_effect_removal,
			self.differential_methylated_pos,
			self.differential_methylated_regions,
			self.principal_comp,
			self.position_enrichment_analysis,
			self.region_enrichment_analysis,
			self.apply_models,
			self.build_model
			]

	def read_idat(self):	
		"""
		Read raw .idat files together with the experiment sample sheet to generate an annotated 
		minfi RGChannelSet object
		"""
		basenames = []

		jobs = []
		
		for sample in self.samples:
			basenames.append(sample.basename)

		jobs.append(
			Job(
					[self.args.samples.name],
					[os.path.join(self.RG_DIR,"rgset.rds"),os.path.join(self.REP_DIR, "data","samples.txt")],
					[
						["read_idat","module_R"],
						["read_idat","module_mugqic_R_packages"]
					],
					command="""\
	mkdir -p {rgchannel_dir}
	mkdir -p {grset_dir}
	mkdir -p report/data
	R --no-save --no-restore <<-EOF
	suppressPackageStartupMessages(library(minfi))
	targets <- read.csv("{samplesheet}")
	targets['Basename'] <- c{basenames} 
	rgset <- read.metharray.exp(target = targets)
	saveRDS(rgset, file = "{output_file}")
	nrows <- sapply("{samplesheet}", function(f) nrow(read.csv(f)) )
	sum <- sum(nrows)
	f <- file("report/data/samples.txt", "w")
	cat(nrows, file=f, sep="\n")
	close(f)

	EOF""".format(
			rgchannel_dir=self.RG_DIR,
			grset_dir = self.GR_DIR,
			samplesheet=self.args.samples.name,
			basenames=tuple(basenames),
			output_file=os.path.join(self.RG_DIR,"rgset.rds")
			), 
			name="read_idat")
		)

		# generate contrast table to be displayed in the final html report
		fill_in_entry = '| {contrast_name} | [{number_controls}]({control_link}) | [{number_cases}]({case_link}) |'
		try:
    			os.stat("output")
		except:
    			os.mkdir("output")  
		for contrast in self.contrasts:
			controls = [sample.name for sample in contrast.controls]
			number_controls = len(controls)
			cases = [sample.name for sample in contrast.treatments]
			number_cases = len(cases)
			with open("output/"+contrast.name+"_controls.txt", 'w+') as controlfile:
				for control in controls:
					controlfile.write(control+"\n")
			with open("output/"+contrast.name+"_cases.txt",'w+') as casefile:
				for case in cases:
					casefile.write(case + "\n")
			
			entry = fill_in_entry.format(
					contrast_name = contrast.name,
					number_controls = number_controls,
					number_cases = number_cases,
					control_link = "data/"+contrast.name+"_controls.txt",
					case_link = "data/"+contrast.name+"_cases.txt"
					)
			with open("output/"+contrast.name+"contrast.txt", 'w+') as table:
				table.write(entry + "\n")

		## generate partial report - introduction

		report_file = os.path.join("report", "Epigene.introduction.md")
		jobs.append(
			Job(
				[self.args.samples.name, os.path.join(self.REP_DIR, "data","samples.txt")],
				[report_file],
				[
					['read_idat', 'module_R'],
					['read_idat', 'module_mugqic_R_packages'],
					['read_idat', 'module_pandoc']
				],
				command="""\
mkdir -p report/data && \\
cp *controls.txt *cases.txt report/data/ && \\
cp *contrast.txt report/ && \\
tail -n +2 {samplesheet} | cut -d , -f2 | sort | uniq -c > report/data/group.txt && \\
group_list_md=`cat report/data/group.txt` && \\
sample_md=`cat report/data/samples.txt` && \\
pandoc \\
--to markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable data_type="{data_type}" \\
--variable date="{date}" \\
--variable sample_number="$sample_md" \\
--variable group_list="$group_list_md" \\
--variable contrast_table="`cat report/*contrast.txt`" \\
  {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}
mkdir -p report/images/ && \\
cp "{input_file}" report/images 
""".format(
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				samplesheet=self.args.samples.name,
				data_type=config.param("report","data_type"),
				date=datetime.datetime.now().strftime("%Y-%m-%d"),
				input_file=config.param("report", "diagram"),
				report_file=report_file
				),
			report_files=[report_file],
			name="intro_report")
			)


		return jobs


	def quality_control(self):
		"""
		Checks quality of data before any differential analysis by generating
		different types of QC plots
		"""
		jobs = []

		if config.param("quality_control","quality_control_plots").lower() == "y":
			job = concat_jobs([
						Job(command="mkdir -p " + self.QC_DIR),
						qc.quality_control_plots(
						os.path.join(self.RG_DIR,"rgset.rds"),
						os.path.join(self.QC_DIR,"qc_plot.png"),
						self.QC_DIR,
						config.param('column_names','group'))])
			job.name="quality_control_plots"

			jobs.append(job)


		return jobs
	

	def normalize(self):

		"""
		Normalize the data from a selection of various normalization methods
		"""		

		jobs = []				

		if config.param("normalize","method").lower() == "funnorm":
			jobs.append(concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.func(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			], name="funnorm_normalization")
			)
		elif config.param("normalize","method").lower() == "illumina":
			jobs.append(concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.illumina(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			], name="illumina_normalization")
			)
		elif config.param("normalize","method").lower() == "swan":
			jobs.append(concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.swan(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			], name="swan_normalization")
			)
		elif config.param("normalize","method").lower() == "noob":
			jobs.append(concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.noob(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			], name="noob_normalization")
			)
		elif config.param("normalize","method").lower() == "raw":
			jobs.append(concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.raw(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			], name="raw_normalization")
			)
		elif config.param("normalize","method").lower() == "quantile":
			jobs.append(concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.quantile(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			], name="quantile_normalization")
			)

		## generating a density bean plot after normalization

		grset_file = os.path.join(self.GR_DIR,"grset_norm.rds")
		densityBeanPlot_norm_file = os.path.join(self.QC_DIR, "densityBeanPlot_norm.png")
		jobs.append(
			Job(
				[grset_file],
				[densityBeanPlot_norm_file],
				[
					['normalize', 'module_R'],
					['normalize', 'module_mugqic_R_packages']
				],
				command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
library(beanplot)
rgset <- readRDS("{grset_file}")
png(file = "{densityBeanPlot_norm_file}", units = "in", width = 8.5, height = 8.5, res=300)
names <- pData(rgset)\$Sample_Name
groups <- pData(rgset)\${group}
gset <- getBeta(rgset)
par(cex.axis=0.7)
densityBeanPlot(gset, sampName = names, sampGroups = groups, main = "Post-normalization density bean plot")
dev.off()
EOF""".format(
		grset_file=grset_file,
		densityBeanPlot_norm_file=densityBeanPlot_norm_file,
		group=config.param('column_names','group')
		),
	name="densityBeanPlot_norm")
		)


		## generating partial report for normalization step
		
		grset_file = os.path.join(self.GR_DIR, "grset_norm.rds")
		report_file = os.path.join(self.REP_DIR, "Epigene.norm.md")
		jobs.append(
			Job(
				[grset_file],
				[report_file],
				[
					['normalize', 'module_R'],
					['normalize', 'module_mugqic_R_packages'],
					['normalize', 'module_pandoc']
				],
				command="""\
mkdir -p report && \\
pandoc \\
--to markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable method="{method}" \\
  {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				grset_file=grset_file,
				method=config.param("normalize","method"),
				report_file=report_file
				),
			report_files=[report_file],
			name="normalize_report")
			)

		return jobs



	def probe_filter(self):
		"""
		This step identifies poor quality probes/CpG sites in the sample that should be excluded from 
		downstream analysis. Various metrics can be chosen - please see documentation for specifics.
		The filtered data is a minfi Genomic Ratio object
		"""		

		jobs = []	
		
		#Find probes that failed to meet the metrics chosen
		#Note: There should be a neater/more efficient way of doing this. Ideally, the jobs do not
		#	   need to be appended together but rather concatenated as a single job. However, this
		#	   does not seem to be working when R runs interatively. Appending is the temporary solution 
		if config.param("probe_filter","find_failures").lower() == "y":
			jobs.append(concat_jobs([
							Job(command="mkdir -p " + self.QC_DIR),
							Job(command="mkdir -p " + self.QC_DIR),
							qc.find_failures(
								os.path.join(self.RG_DIR,"rgset.rds"),
								os.path.join(self.QC_DIR, "failed_probes.csv"),
								config.param("find_failures","padjust_method"),
								config.param("find_failures","pvalue"),
								config.param("find_failures","percent_failed_threshold"))
							], name="finding_failures")
			)

			csv_file = os.path.join(self.QC_DIR,"failed_probes.csv")
			report_file = os.path.join("report", "Epigene.failed_probes.md")

			jobs.append(
							Job(
								[csv_file],
								[report_file],
								[
									['probe_filter', 'module_R'],
									['probe_filter', 'module_mugqic_R_packages'],
									['probe_filter', 'module_pandoc']
								],
								command="""\
grep --regexp="$" --count "{csv_file}" > report/data/failed_removed.txt && \\
mkdir -p report && \\
pandoc \\
--to markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable percent_failed_threshold="{percent_failed_threshold}" \\
--variable pvalue="{pvalue}" \\
--variable failed_removed=`cat report/data/failed_removed.txt` \\
  {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				csv_file=csv_file,
				percent_failed_threshold=config.param("find_failures", "percent_failed_threshold"),
				pvalue=config.param("find_failures", "pvalue"),
				report_file=report_file
				),
			report_files=[report_file],
			name="failed_probes_report")
			)


		if ((config.param("probe_filter","cross_reactive_probes").lower() == "y") and (config.param("report","data_type").lower() == "450k")):  
			jobs.append(concat_jobs([
							Job(command="mkdir -p " + self.QC_DIR),
							qc.cross_reactive_probes(
								config.param("cross_reactive_probes","crpfile450K"),
								os.path.join(self.QC_DIR, "cross_reactive_probes.csv"))
						], name="cross_reactive_probes")
			)

			csv_file = os.path.join(self.QC_DIR,"cross_reactive_probes.csv")
			report_file = os.path.join("report", "Epigene.cross_reactive_probes.md")

			jobs.append(
							Job(
								[csv_file],
								[report_file],
								[
									['probe_filter', 'module_R'],
									['probe_filter', 'module_mugqic_R_packages'],
									['probe_filter', 'module_pandoc']
								],
								command="""\
mkdir -p report && \\
pandoc \\
--to markdown \\
--template {report_template_dir}/{basename_report_file} \\
  {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				report_file=report_file
				),
			report_files=[report_file],
			name="cross_reactive_probes_report")
			)


		if ((config.param("probe_filter","cross_reactive_probes").lower() == "y") and (config.param("report","data_type").lower() == "850k")):  
			jobs.append(concat_jobs([
							Job(command="mkdir -p " + self.QC_DIR),
							Job(command="mkdir -p " + self.QC_DIR),
							qc.cross_reactive_probes(
								config.param("cross_reactive_probes","crpfile850K"),
								os.path.join(self.QC_DIR, "cross_reactive_probes.csv"))
						], name="cross_reactive_probes")
			)

		if config.param("probe_filter","find_snps_chen").lower() == "y":
			jobs.append(concat_jobs([
							Job(command="mkdir -p " + self.QC_DIR),
							qc.find_snps_chen(
								config.param("find_snps_chen","pCpGfile"),
								os.path.join(self.QC_DIR, "snp_probes_chen.csv"),
								config.param("find_snps_chen","minimal_distance_to_SNPs"),
								config.param("find_snps_chen","minimal_minor_allele_frequency"))
							], name="find_snps_chen")
			)

		if config.param("probe_filter","find_extremes").lower() == "y":
			jobs.append(concat_jobs([
						Job(command="mkdir -p " + self.QC_DIR),
						Job(command="mkdir -p " + self.QC_DIR),
						qc.find_extremes(
							os.path.join(self.RG_DIR,"rgset.rds"),
							os.path.join(self.QC_DIR,"failed_probes.csv"))
						], name="find_extremes")
			)
		if config.param("probe_filter","find_chromosomes").lower() == "y":
			jobs.append(concat_jobs([
						Job(command="mkdir -p " + self.QC_DIR),
						Job(command="mkdir -p " + self.QC_DIR),
						qc.find_chromosomes(
							os.path.join(self.RG_DIR,"rgset.rds"),
							os.path.join(self.QC_DIR,"chr_probes.csv"),
							config.param("find_chromosomes","chromosomes",type = "list"))
						], name="find_chromosomes")
			)
		if config.param("probe_filter","find_snps_minfi").lower() == "y":
			jobs.append(concat_jobs([
						Job(command="mkdir -p " + self.QC_DIR),
						Job(command="mkdir -p " + self.QC_DIR),
						qc.find_snps_minfi(
							os.path.join(self.RG_DIR,"rgset.rds"),
							os.path.join(self.QC_DIR,"snp_probes_minfi.csv"))
						], name="find_snps_minfi")
			)
			csv_file = os.path.join(self.QC_DIR,"snp_probes_minfi.csv")
			report_file = os.path.join("report", "Epigene.snp_minfi.md")

			jobs.append(
							Job(
								[csv_file],
								[report_file],
								[
									['probe_filter', 'module_R'],
									['probe_filter', 'module_mugqic_R_packages'],
									['probe_filter', 'module_pandoc']
								],
								command="""\
grep --regexp="$" --count "{csv_file}" > report/data/snp_removed.txt && \\
mkdir -p report && \\
pandoc \\
--to markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable snp_removed=`cat report/data/snp_removed.txt` \\
  {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				csv_file=csv_file,
				report_file=report_file
				),
			report_files=[report_file],
			name="snp_probes_minfi_report")
			)
		
		#Remove chosen probes from normalized data		
		jobs.append(concat_jobs([
			Job(command="mkdir -p " + self.QC_DIR),
			qc.merge_remove(
				os.path.join(self.GR_DIR,"grset_norm.rds"),
				os.path.join(self.GR_DIR,"grset_norm_qc.rds"),
				self.QC_DIR
			)
		], name="probe_filter"))


		## generating partial report for probe filtering step
		
		grset_file = os.path.join(self.GR_DIR, "grset_norm_qc.rds")
		report_file = os.path.join(self.REP_DIR, "Epigene.probe_filter.md")
		jobs.append(
			Job(
				[grset_file],
				[report_file],
				[
					['probe_filter', 'module_R'],
					['probe_filter', 'module_mugqic_R_packages'],
					['probe_filter', 'module_pandoc']
				],
				command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
grset <- readRDS("{grset_file}")
total <- nrow(grset) - 1
f <- file("report/data/probes.txt", "w")
cat(total, file=f, sep="\n")
close(f)
EOF
mkdir -p report && \\
pandoc \\
--to markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable total=`cat report/data/probes.txt` \\
  {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				grset_file=grset_file,
				report_file=report_file
				),
			report_files=[report_file],
			name="probe_filter_report")
			)

						
		return jobs


	def batch_effect_removal(self):
		"""
		Reduces batch effects between different batches of the sample
		"""
		jobs = []
		grset_file = os.path.join(self.GR_DIR,"grset_norm_qc.rds")
		output_file = os.path.join(self.GR_DIR,"matrix.rds")

		if config.param("batch_effect_removal","batch_effect_removal") == "y":
			batch_col = config.param("column_names", "batch")

			job = Job(
							[grset_file],
							[output_file],
							[
								["batch_effect_removal","module_R"],
								["batch_effect_removal","module_mugqic_R_packages"]
							],
							command="""\
	mkdir -p {gr_dir} && \\
	R --no-save --no-restore <<-EOF
	suppressPackageStartupMessages(library(minfi))
	library(sva)
	gset <- readRDS("{grset_file}")
	edata <- getBeta(gset)
	pheno <- pData(gset)
	batch <- pheno\${batch_col}
	modcombat <- model.matrix(~1, data=pheno)
	combat_edata <- ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plot=FALSE)
	saveRDS(combat_edata, "{output_file}")
	EOF
	mkdir -p report && \\
	pandoc \\
	--to markdown \\
	--template {report_template_dir}/{basename_report_file} \\
	--variable batch_column = {batch_col} \\
	  {report_template_dir}/{basename_report_file} \\
	{report_template_dir}/{basename_report_file} \\
	> {report_file}""".format(
					gr_dir=self.GR_DIR,
					grset_file=grset_file,
					batch_col = batch_col,
					output_file=output_file
					),
					name="batch_effect_removal")

			jobs.append(job)
		else:

                        jobs.append(Job([grset_file],
					[output_file],
					command = "cp " + grset_file + " " +  output_file,
					 name="no_batch_effect"))
				

		return jobs


			

	def cell_counts(self):
		"""
		Estimate confounding levels between phenotype and cell type composition. The cell type
		composition of each sample is returned. This step should only be included for blood samples 
		"""

		jobs = []

		for sample in self.samples:
			if not hasattr(sample, 'cell_type'):
				log.info("Warning: Sample " + sample.name + "should be a blood sample!")
		job = Job(
					[os.path.join(self.RG_DIR,"rgset.rds")],
					[os.path.join(self.QC_DIR,"cell_counts.csv")],
					[
						["cell_counts","module_R"],
						["cell_counts","module_mugqic_R_packages"]
					],
					command="""\
	mkdir -p {qc_dir}
	R --no-save --no-restore <<-EOF
	suppressPackageStartupMessages(library(minfi))
	library(FlowSorted.Blood.450k)
	rgset <- readRDS("{rgset_file}")
	rgset\${sex} <- as.character(rgset\${sex})
	cellCounts <- estimateCellCounts(rgset)
	rownames(cellCounts) <- pData(rgset)\$Sample_Name
	write.csv(as.data.frame(cellCounts), file="{output_file}")

	EOF""".format(
			qc_dir=self.QC_DIR,
			rgset_file=os.path.join(self.RG_DIR,"rgset.rds"),
			output_file=os.path.join(self.QC_DIR,"cell_counts.csv"),
			sex=config.param('column_names','sex')
			),
			name="cell_counts")

		jobs.append(job)


		## generating partial report for cell counts

		cell_counts_file = os.path.join(self.QC_DIR, "cell_counts.csv")
		report_file = os.path.join(self.REP_DIR, "Epigene.cell_counts.md")
		jobs.append(
			Job(
				[cell_counts_file],
				[report_file],
				[['cell_counts', 'module_pandoc']],
				command="""\
mkdir -p report && \\
cat {cell_counts_file} | sed 's/,/\t/g' | cut -f1-7 | sed 's/\"//g' > report/temp_cell_counts.tsv && \\
cell_counts_table_md=`head -7 report/temp_cell_counts.tsv | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print "Sample Name|CD8+ T-cells|CD4+ T-cells|Natural Killer cells|B-cells|Monocytes|Granulocytes|"; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, sprintf("%.4f", $2), sprintf("%.4f", $3), sprintf("%.4f", $4), sprintf("%.4f", $5), sprintf("%.4f", $6), sprintf("%.4f", $7)}}}}' ;` && \\
rm report/temp_cell_counts.tsv && \\ 
pandoc \\
--to markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable cell_counts_table="$cell_counts_table_md" \\
  {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				cell_counts_file=cell_counts_file,
				report_file=report_file
				),
			report_files=[report_file],
			name="cell_counts_report")
			)

		return jobs


	def differential_methylated_pos(self):
		"""
		This step finds CpG sites that are differentially methylated by sample group (i.e control vs case).
		There can be multiple control-case pairs(contrasts) and this information is provided in the design file.
		This step can be run either using the minfi package or the limma package or both to compare results.
	
		When run using the minfi package, this step is essentially a wraparound the minfi function dmpFinder and uses a F-test to generate p-values for each site. The output is a table of differentially methylated CpG sites with associated information
		such as chromsomal position, average methylation difference (from controls), p-value, and whether it is
		confounded by blood composition (for blood samples only).

		When run using the limma package, this step uses the limma functions lmFit and eBayes. In addition, there is an option to specify
		which other factors (like age, sex etc) to consider when running differential analysis.
		"""
		jobs = []
		grset_file = os.path.join(self.GR_DIR,"grset_norm_qc.rds")
		input_rds_file = os.path.join(self.GR_DIR,"matrix.rds")
		minfi_data_files = []
		limma_data_files = []
		fill_in_entry = '**Table for contrast {contrast_name} (partial table;** [download full table]({link}))'

			# find differentially methylated positiosn using minfi package
		if config.param("differential_methylated_pos","minfi").lower() == "y":
			for contrast in self.contrasts:
				dmps_file = os.path.join(self.CPG_DIR,contrast.name+ "_minfi_dmps.csv")
            # Abort analysis if not enough samples (Will cause dmpFinder to throw an error)
				if len(contrast.controls) == 0 or contrast.treatments == 0 or len(contrast.controls)+len(contrast.treatments) < 2:
					log.warn("Insufficient sample size to compare case and control. Skipping contrast: " + contrast.name)
					continue

				minfi_data_file = os.path.join(self.REP_DIR, contrast.name + "_minfi_data.txt")
				report_entry = fill_in_entry.format(
													contrast_name = contrast.name,
													link = "data/cpg/"+os.path.join(os.path.basename(dmps_file))
												)
				job = concat_jobs([
					Job(command="mkdir -p " + self.CPG_DIR),
					dmp.minfi_differential_methylated_pos(
						contrast,
						grset_file,
						input_rds_file,
						dmps_file,
						config.param("minfi_differential_methylated_pos","padjust_method"),
						config.param("minfi_differential_methylated_pos","pvalue"),
						config.param("minfi_differential_methylated_pos","delta_beta_threshold"))
					])	

				job.name="minfi_differential_methylated_pos."+contrast.name

				jobs.append(job)

				## generating partial report for differentially methylated positions
				minfi_data_files.append(minfi_data_file)
				job = 	Job(
						[dmps_file],
						[minfi_data_file],
						command="""\
		grep --regexp="$" --count "{dmps_file}" > report/data/{contrast_name}_minfi_total_dmp.txt && \\
		cat {dmps_file} | sed 's/,/\t/g' | cut -f1-8 | sed 's/\"//g' > report/{contrast_name}_temp_dmps.tsv && \\
		dmps_table_md=`head -7 report/{contrast_name}_temp_dmps.tsv | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print "CpG Site|Avg Case Beta|Avg ControlBeta|Avg Delta Beta|Chromosome|Gene|P-value|Adj.P-value|"; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, sprintf("%.4f", $2), sprintf("%.4f", $3), sprintf("%.4f", $4), $5, $6, sprintf("%.4f", $7),sprintf("%.4f", $8)}}}}' ;` && \\
		rm report/{contrast_name}_temp_dmps.tsv && \\
		echo "{report_entry} \n\n$dmps_table_md\n" >> {data_file} && \\
		""".format(
				report_entry = report_entry,	
				dmps_file = dmps_file,
				data_file = minfi_data_file,
				contrast_name = contrast.name
			),
					name = contrast.name+"_minfi_dmps_table"
					)

				jobs.append(job)
		
			
			report_file = os.path.join("report", "Epigene.minfi_dmps.md")	
			job=Job(
					minfi_data_files,
					[report_file],
					[['differential_methylated_pos','module_pandoc']],
					command="""\
pandoc \\
--to markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable padjust_method="{padjust_method}" \\
--variable pvalue="{pvalue}" \\
--variable delta_beta_threshold="{delta_beta_threshold}" \\
--variable data="`cat {rep_dir}/*minfi_data.txt`" \\
  {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
				rep_dir = self.REP_DIR,
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				padjust_method=config.param("minfi_differential_methylated_pos","padjust_method"),
				pvalue=config.param("minfi_differential_methylated_pos","pvalue"),
				delta_beta_threshold=config.param("minfi_differential_methylated_pos","delta_beta_threshold"),
				report_file=report_file
			),
			report_files=[report_file],
			name="minfi_differential_methylated_pos_report")
			jobs.append(job)
			

			# Find dmps using the limma package
		if config.param('differential_methylated_pos', 'limma').lower() == 'y':
			for contrast in self.contrasts:
				dmps_file = os.path.join(self.CPG_DIR, contrast.name + "_limma_dmps.csv")
            
				# Abort analysis if not enough samples (Will cause dmpFinder to throw an error)
				if len(contrast.controls) == 0 or contrast.treatments == 0 or len(contrast.controls)+len(contrast.treatments) < 2 :
					log.warn("Insufficient sample size to compare case and control. Skipping contrast: " + contrast.name)
					continue
				



				limma_data_file = os.path.join(self.REP_DIR, contrast.name + "_limma_data.txt")
				report_entry = fill_in_entry.format(
                                                    contrast_name = contrast.name,
                                                    link = "data/cpg/"+os.path.join(os.path.basename(dmps_file))
                                                )

				#Add additional parameters to the design matrix based on options set in configuration file.
				if config.param("limma_differential_methylated_pos", "additional_params") == 'y':
					params =  config.param("limma_differential_methylated_pos", "design_matrix_columns")
					design_matrix_columns = "+Age+Sex" if params == "Age,Sex" else "+"+params
				else:
					design_matrix_columns = ""	
	
				job = concat_jobs([
					Job(command="mkdir -p " + self.CPG_DIR),
					dmp.limma_differential_methylated_pos(
						contrast,
						grset_file,
						input_rds_file,
						dmps_file,
						config.param("limma_differential_methylated_pos", "padjust_method"),
						config.param("limma_differential_methylated_pos", "pvalue"),
						config.param("limma_differential_methylated_pos","delta_beta_threshold"),
						design_matrix_columns
					)])

				job.name="limma_differential_methylated_pos."+contrast.name

				jobs.append(job)

				#generate data tables
				limma_data_files.append(limma_data_file)
				jobs.append(
						Job(
							[dmps_file],
							[limma_data_file],
							command="""\
	grep --regexp="$" --count "{dmps_file}" > report/data/{contrast_name}_limma_total_dmp.txt && \\
        cat {dmps_file} | sed 's/,/\t/g' | cut -f1-10 | sed 's/\"//g' > report/{contrast_name}_temp_limma_dmps.tsv && \\
        dmps_table_md=`head -7 report/{contrast_name}_temp_limma_dmps.tsv | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print "CpG Site|Gene|Chromosome|logFC|AveExpr|avgDeltaBeta|t|P.Value|adj.P.Val|B|"; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1,$8,$9, sprintf("%.4f", $2), sprintf("%.4f", $3),sprintf("%.4f", $10), sprintf("%.4f", $4),sprintf("%.4f", $5),sprintf("%.4f", $6), sprintf("%.4f", $7)}}}}' ;` && \\
        rm report/{contrast_name}_temp_limma_dmps.tsv && \\
        echo "{report_entry} \n\n$dmps_table_md\n" >> {data_file} && \\
		""".format(
				dmps_file = dmps_file,
				contrast_name = contrast.name,
				report_entry = report_entry,
				data_file = limma_data_file	
					),
							name = contrast.name+"_limma_dmps_table"
							))
				 ## generating partial report for differentially methylated positions

				report_file = os.path.join("report", "Epigene.limma_dmps.md")
				jobs.append(
					Job(
						limma_data_files,
						[report_file],
						[['differential_methylated_pos', 'module_pandoc']],
						command="""\
		pandoc \\
		--to markdown \\
		--template {report_template_dir}/{basename_report_file} \\
		--variable limma_padjust="{padjust_method}" \\
		--variable limma_pvalue="{pvalue}" \\
		--variable limma_delta_beta="{delta_beta}" \\
		--variable limma_data="`cat report/*limma_data.txt`" \\
		  {report_template_dir}/{basename_report_file} \\
		{report_template_dir}/{basename_report_file} \\
		> {report_file}""".format(
						report_template_dir=self.report_template_dir,
						basename_report_file=os.path.basename(report_file),
						padjust_method=config.param("limma_differential_methylated_pos","padjust_method"),
						pvalue=config.param("limma_differential_methylated_pos","pvalue"),
						report_file=report_file,
						delta_beta=config.param("limma_differential_methylated_pos","delta_beta_threshold")
					),
					report_files=[report_file],
					name="limma_differential_methylated_pos_report")
					)

		return jobs

	
	def differential_methylated_regions(self):
		"""
		This step finds differentially methalated regions for each contrast using the bumphunting algorithm in minfi bumphunter.
		The algorithm looks for genomic regions that are differentially methylated between two conditions.
		"""

		jobs = []
		grset_file = os.path.join(self.GR_DIR,"grset_norm_qc.rds")
		data_files = []
		fill_in_entry = '**Table for contrast {contrast_name} (partial table;** [download full table]({link}))'
		for contrast in self.contrasts:
			dmrs_file = os.path.join(self.CPG_DIR, contrast.name+"_dmrs.csv")
			data_file = os.path.join(self.REP_DIR, contrast.name+"_dmrs.txt")
			report_entry = fill_in_entry.format(contrast_name=contrast.name, link= "data/cpg/"+os.path.join(os.path.basename(dmrs_file)))

			job = concat_jobs([
				Job(command="mkdir -p " + self.CPG_DIR),
				dmp.find_dmr(
					contrast,
					grset_file,
					dmrs_file,
					config.param("differential_methylated_regions","nullMethod"),
					config.param("differential_methylated_regions","permutations"),
					config.param("differential_methylated_regions","delta_beta_threshold"),
					config.param("differential_methylated_regions", "pvalue"),
					config.param("differential_methylated_regions", "LCutoff"),
					config.param("differential_methylated_regions", "cores"))
				])

			job.name="differential_methylated_regions."+contrast.name
			data_files.append(data_file)
			jobs.append(job)
			jobs.append(Job(
						[dmrs_file],
						[data_file],
						command="""\
		cat {dmrs_file} | sed 's/,/\t/g' | cut -f2-5,10,12,13,16 | sed 's/\"//g' > {contrast_name}_temp_dmrs.tsv && \\
    	dmrs_table_md=`head -7 {contrast_name}_temp_dmrs.tsv | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print "chr|start|end|gene|value|L|p.value|fwer|";print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}}else {{print $1, sprintf("%\\47d", $2), sprintf("%\\47d", $3),$8, sprintf("%.4f", $4),$5, sprintf("%.4f", $6), sprintf("%.4f", $7)}}}}' ;` && \\ 
    	rm {contrast_name}_temp_dmrs.tsv && \\
		echo "{report_entry} \n\n$dmrs_table_md\n" >> {data_file} && \\
		""".format(
					contrast_name = contrast.name,
					report_entry = report_entry,
					data_file = data_file,
					dmrs_file = dmrs_file,
					),
						name = contrast.name + "_dmrs_table"))


		## generating partial report for differentially methylated regions

		report_file = os.path.join("report", "Epigene.dmrs.md")
		jobs.append(
				Job(
					data_files,
					[report_file],
					[['differential_methylated_regions', 'module_pandoc']],
					command="""\
		mkdir -p report && \\
		pandoc \\
		--to markdown \\
		--template {report_template_dir}/{basename_report_file} \\
		--variable nullMethod="{nullMethod}" \\
		--variable permutations="{permutations}" \\
		--variable delta_beta_threshold="{delta_beta_threshold}" \\
		--variable data="`cat report/*dmrs.txt`" \\
		{report_template_dir}/{basename_report_file} \\
		{report_template_dir}/{basename_report_file} \\
		> {report_file}""".format(
					report_template_dir=self.report_template_dir,
					basename_report_file=os.path.basename(report_file),
					nullMethod=config.param("differential_methylated_regions","nullMethod"),
					permutations=config.param("differential_methylated_regions","permutations"),
					delta_beta_threshold=config.param("differential_methylated_regions","delta_beta_threshold"),
					report_file=report_file
				),
				report_files=[report_file],
				name="differential_methylated_regions_report")
				)

		return jobs

		
	def principal_comp(self):
		"""
		This step generates a principal component analysis plot of the first to PC in the data. Samples are
		differentiated by their group (i.e controls or cases)
		"""
		jobs=[]
		png_files = []
		input_rds_file = os.path.join(self.GR_DIR, "grset_norm_qc.rds")
		for contrast in self.contrasts:
			#generate plots for post normalization data
			input_file = os.path.join(self.GR_DIR, "grset_norm_qc.rds")
			output_file = os.path.join("report", contrast.name+"_norm_pcaplot")
			png_files.append(output_file + "-1.png")
			jobname = "pca_norm."+contrast.name
			jobs.append(
				pca.pca_plot(input_rds_file,input_file, output_file, contrast,jobname))
		
			#for minfi generated dmps	
			if config.param("differential_methylated_pos","minfi") == 'y':
				input_file = os.path.join(self.CPG_DIR, contrast.name+"_minfi_dmps.csv")
				output_file = os.path.join("report",contrast.name+"_minfi_pcaplot")
				png_files.append(output_file + "-1.png")
				jobname = "pca_minfi."+contrast.name
				jobs.append(
					pca.pca_plot(input_rds_file,input_file, output_file, contrast,jobname))
	
			#for limma generated dmps
			if config.param("differential_methylated_pos","limma") == 'y':
				input_file = os.path.join(self.CPG_DIR, contrast.name+"_limma_dmps.csv")
				output_file = os.path.join("report",contrast.name+"_limma_pcaplot")
				png_files.append(output_file + "-1.png")
				jobname = "pca_limma."+contrast.name
				jobs.append(
					pca.pca_plot(input_rds_file,input_file, output_file, contrast,jobname))

		
		pca_file = ("report/pcaplot.pdf")
		report_file = os.path.join("report", "Epigene.pcaplot.md")

		
		job = Job(
				png_files,
				[report_file],
				[['principal_comp', 'module_pandoc']],
				command="""\
mkdir -p report && \\
mkdir -p report/images && \\
cp report/*pcaplot-1.png report/images/ && \\
pandoc \\
--to markdown \\
--template {report_template_dir}/{basename_report_file} \\
  {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				pca_file=pca_file,
				report_file=report_file
			),
			report_files=[report_file],
			name="principal_comp_report")
		jobs.append(job)
			

		
		return jobs

	def position_enrichment_analysis(self):
		"""
		This step uses the R package LOLA to test for overlap of positions either identified in the differential_methylated_pos step or provided by the user
		against region sets selected previously. These region sets can be selected from the LOLAcore database, the default database used
		by this pipeline or any other bed files selected by the user. The information on how to develop your own region sets can be found
		here: http://databio.org/regiondb
		"""

		report_file = 'report/Epigene.enrichment_analysis.md'
		report_data = 'data/enrichment_analysis'
		grset_file = os.path.join(self.GR_DIR,"grset_norm_qc.rds")
		fill_in_entry= "| Positions ({source}) | {contrast_name} | [download link]({results_link}) |"
		jobs = []
		data_files = []
		
		#if the cpgs should be taken from a file provided by the user
		if config.param('enrichment_analysis','own_file') == 'y':
			analysis_file = os.path.join(self.REP_DIR, "data","enrichment_analysis","cpg_position_enrichment_analysis.csv")
			cpgs_file = config.param('enrichment_analysis','cpgs_file')
			grset_file = os.path.join(self.GR_DIR,"grset_norm.rds")
			jobs.append(
				enrichment.pos_enr(
                                    grset_file,
                                    cpgs_file,
                                    config.param('position_enrichment_analysis', 'cluster_cpu').split('=')[-1],
                                    config.param('position_enrichment_analysis', 'LOLA_dir'),
                                    tuple(config.param('position_enrichment_analysis', 'LOLA_bed_files', type='dirpathlist')),
                                    config.param('position_enrichment_analysis', 'assembly'),
                                    tuple(config.param('position_enrichment_analysis', 'collection', type='list')),
                                    analysis_file))

		else:

				for contrast in self.contrasts:
					if config.param("differential_methylated_pos", "minfi") == 'y':

							data_file = os.path.join(self.REP_DIR, contrast.name + "_minfi_pos_enr.txt")
							analysis_file = os.path.join(self.REP_DIR, "data","enrichment_analysis",contrast.name+"_minfi_position_enrichment_analysis.csv")

							report_entry = fill_in_entry.format(
							contrast_name = contrast.name,
							source = "minfi",
							results_link=os.path.join(report_data, os.path.basename(analysis_file))),

							dmps_file = os.path.join(self.CPG_DIR,contrast.name+"_minfi_dmps.csv")

							jobs.append(
									enrichment.pos_enr(
										grset_file,
										dmps_file,
										config.param('position_enrichment_analysis', 'cluster_cpu').split('=')[-1],
										config.param('position_enrichment_analysis', 'LOLA_dir'),
										tuple(config.param('position_enrichment_analysis', 'LOLA_bed_files', type='dirpathlist')),
										config.param('position_enrichment_analysis', 'assembly'),
										tuple(config.param('position_enrichment_analysis', 'collection', type='list')),
										analysis_file,
										contrast.name, "minfi"))
							#generate data tables
							data_files.append(data_file)
							jobs.append(
									Job(
										[analysis_file],
										[data_file],
										command="""\
						cat {analysis_file} | sed 's/,/\t/g' | cut -f2-5,12-16 | sed 's/\"//g' > {contrast_name}_minfi_temp_pos_enr.tsv && \\
						pos_enr_table_md=`head -7 {contrast_name}_limma_temp_pos_enr.tsv | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print "dbSet|collection|pValueLog|b|c|d|description|";print "-----|-----:|-----:|-----:|-----:|-----:|-----:"}}else{{print $1, sprintf($2), sprintf( $3), sprintf($5), sprintf( $6), sprintf($7), $8}}}}' ;` && \\ 
						echo "{report_entry}" >> {data_file} && \\
						""".format(
						contrast_name = contrast.name,
						analysis_file = analysis_file,
						report_entry = report_entry[0],
						data_file = data_file
						),
										name = contrast.name+"minfi_pos_enr_table"))

					if config.param("differential_methylated_pos", "limma") == 'y':
						data_file = os.path.join(self.REP_DIR, contrast.name + "_limma_pos_enr.txt")
						analysis_file = os.path.join(self.REP_DIR, "data","enrichment_analysis",contrast.name+"_limma_position_enrichment_analysis.csv")
						report_entry = fill_in_entry.format(
						contrast_name = contrast.name,
						source = "limma",
						results_link=os.path.join(report_data, os.path.basename(analysis_file))),
						dmps_file = os.path.join(self.CPG_DIR,contrast.name+"_limma_dmps.csv")
						jobs.append(
									enrichment.pos_enr(
											grset_file,
											dmps_file,
											config.param('position_enrichment_analysis', 'cluster_cpu').split('=')[-1],
											config.param('position_enrichment_analysis', 'LOLA_dir'),
											tuple(config.param('position_enrichment_analysis', 'LOLA_bed_files', type='dirpathlist')),
											config.param('position_enrichment_analysis', 'assembly'),
											tuple(config.param('position_enrichment_analysis', 'collection', type='list')),
											analysis_file,
											contrast.name,"limma"))
						#generate data tables
						data_files.append(data_file)
						jobs.append(
								Job(
									[analysis_file],
									[data_file],
									command="""\
						cat {analysis_file} | sed 's/,/\t/g' | cut -f2-5,12-16 | sed 's/\"//g' > {contrast_name}_limma_temp_pos_enr.tsv && \\
						pos_enr_table_md=`head -7 {contrast_name}_limma_temp_pos_enr.tsv | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print "dbSet|collection|pValueLog|b|c|d|description|";print "-----|-----:|-----:|-----:|-----:|-----:|-----:"}}else{{print $1, sprintf($2), sprintf( $3), sprintf($5), sprintf( $6), sprintf($7), $8}}}}' ;` && \\ 
						echo "{report_entry}" >> {data_file} && \\
						""".format(
									contrast_name = contrast.name,
									analysis_file = analysis_file,
									report_entry = report_entry[0],
									data_file = data_file
									),
									name = contrast.name+"_limma_pos_enr_table"
									))

				#Generate final report
				jobs.append(
						Job(
							data_files,
							[report_file],
							[['position_enrichment_analysis','module_pandoc']],
							command="""\
				mkdir -p {data_dir} && \\
			table=$(cat report/*enr.txt) && \\
			pandoc \\
				{report_template_dir}/{basename_report_file} \\
				--template {report_template_dir}/{basename_report_file} \\
				--variable data_table="$table" \\
				--to markdown > {report_file}
				""".format(
							data_dir = os.path.join(self.REP_DIR,report_data),
							data_files = data_files,
							report_file = report_file,
							report_template_dir=self.report_template_dir,
							basename_report_file=os.path.basename(report_file)
						),
							report_files = [report_file],
							name = 'report_position_enrichment'
							))
				

		return jobs

	def region_enrichment_analysis(self):
		"""
		Similar to position_enrichment_analysis, this step uses the R package LOLA to test for 
		enriched regions identified in the differential_methylated_regions step.
		"""

		report_file = 'report/Epigene.enrichment_analysis.md'
		report_data = 'data/enrichment_analysis'
		fill_in_entry= '| Regions | {contrast_name} | [download link]({results_link}) |'
		jobs = []
		data_files = []
		for contrast in self.contrasts:
			data_file = os.path.join(self.REP_DIR, contrast.name+"_reg_enr.txt")
			data_files.append(data_file)
		for contrast in self.contrasts:
			grset_file = os.path.join(self.GR_DIR,"grset_norm_qc.rds")
			data_file = os.path.join(self.REP_DIR, contrast.name+"_reg_enr.txt")
			analysis_file = os.path.join(self.REP_DIR,"data","enrichment_analysis",contrast.name+"_region_enrichment_analysis.csv")

			report_entry = fill_in_entry.format(
			contrast_name=contrast.name,
			results_link=os.path.join(report_data, os.path.basename(analysis_file)))

			dmrs_file = os.path.join(self.CPG_DIR,contrast.name+"_dmrs.csv")

			command="""\
	mkdir -p report/data/enrichment_analysis && \\
	R --vanilla <<-'EOF'
	suppressPackageStartupMessages(library(BiSeq))
	suppressPackageStartupMessages(library(GenomicRanges))
	suppressPackageStartupMessages(library(data.table))
	suppressPackageStartupMessages(library(doParallel))
	suppressPackageStartupMessages(library(LOLA))
	# suppressPackageStartupMessages(library(matrixStats))
	# suppressPackageStartupMessages(library(minfi))
	#source(file.path(Sys.getenv('R_TOOLS'), 'LOLAsearch.R'))
	source('/hpf/largeprojects/ccmbio/jonBarenboim/mugqic_tools/R-tools/LOLAsearch.R')
	registerDoParallel(cores={cores})
	grset <- readRDS('{grset_file}')
	dmrs <- read.csv('{dmrs_file}', as.is = TRUE)
	universe <- rowRanges(grset)
	userset <- granges(GRanges(dmrs))
	# The universe should contain every region in the dmrs. If not, something probably went wrong in earlier steps
	tryCatch(
		{{ checkUniverseAppropriateness(userset, universe) }},
		warning = function(w) {{ stop("userset is not a subset of universe") }},
		error = function(e) {{ stop("userset is not a subset of universe") }})
	LOLAcoreDB <- suppressWarnings(suppressMessages(buildRegionDB(rootdir='{LOLA_root}', genome='{LOLA_genome}', 
					collection=c{LOLA_collection}, filename=c(), description=c(), any=c())))
	userFileDB <- foreach(dir=c{LOLA_user_dirs}, .combine=mergeRegionDBs) %dopar% {{
		files <- Sys.glob(paste(dir, "/*", sep=""))
		tmpDB <- list()
		tmpDB$regionGRL <- suppressWarnings(suppressMessages(readCollection(files)))
		tmpDB$collectionAnno <- data.table(collectionname=basename(dir), collector=NA, date=NA,
			source=dir, description=paste('User supplied .bed files in directory', dir))
		tmpDB$regionAnno <- data.table(filename=basename(files), cellType=NA, description=NA, tissue=NA, dataSource=NA,
			antibody=NA, treatment=NA, collection=basename(dir), size=sapply(tmpDB$regionGRL, length))
		tmpDB
	}}
	regionDB <- mergeRegionDBs(LOLAcoreDB, userFileDB)
	LOLAresult <- runLOLA(userset, universe, regionDB)
	write.csv(LOLAresult, file='{analysis_file}', row.names=FALSE)
	EOF""".format(
				directory=os.path.dirname(analysis_file),
				entry=report_entry,
				cores=config.param('position_enrichment_analysis', 'cluster_cpu').split('=')[-1],
				grset_file=grset_file,
				dmrs_file=dmrs_file,
				LOLA_root=config.param('region_enrichment_analysis', 'LOLA_dir'),
				LOLA_user_dirs=tuple(config.param('position_enrichment_analysis', 'LOLA_bed_files', type='dirpathlist')),
				LOLA_genome=config.param('region_enrichment_analysis', 'assembly'),
				LOLA_collection=tuple(config.param('region_enrichment_analysis', 'collection', type='list')),
				analysis_file=analysis_file,
				data_dir=os.path.join('report', report_data),
				zip_file=os.path.join('report', report_data, os.path.basename(report_data) + '.zip'),
				report_template_dir=self.report_template_dir,
				basename_report_file=os.path.basename(report_file),
				report_file=report_file)

			job = Job(
					[dmrs_file, grset_file],
					[analysis_file],
					[
						["region_enrichment_analysis", "module_R"],
						["region_enrichment_analysis", "module_pandoc"],
						["region_enrichment_analysis", "module_mugqic_tools"]
					],
					command=command,
					report_files=[report_file],
					name="region_enrichment_analysis."+contrast.name)
			jobs.append(job)
			#generate data tables
			jobs.append(
					Job(
						[analysis_file],
						[data_file],
						command="""\
                cat {analysis_file} | sed 's/,/\t/g' | cut -f2-5,12-16 | sed 's/\"//g' > {contrast_name}_temp_reg_enr.tsv && \\
                reg_enr_table_md=`head -7 {contrast_name}_temp_reg_enr.tsv | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print "dbSet|collection|pValueLog|b|c|d|description|";print "-----|-----:|-----:|-----:|-----:|-----:|-----:"}}else{{print $1, sprintf($2), sprintf( $3), sprintf($5), sprintf( $6), sprintf($7), $8}}}}' ;` && \\ 
        echo "{report_entry}" >> {data_file} && \\
        """.format(
					contrast_name = contrast.name,
					analysis_file = analysis_file,
					report_entry = report_entry,
					data_file = data_file
                                        ),
						name = contrast.name+"_reg_enr_table"
                                        ))



		#Generate final report
		jobs.append(
					Job(
                                        data_files,
                                        [report_file],
                                        [['region_enrichment_analysis','module_pandoc']],
                                        command="""\
                mkdir -p {data_dir} && \\
    table=$(cat report/*enr.txt) && \\
    pandoc \\
        {report_template_dir}/{basename_report_file} \\
        --template {report_template_dir}/{basename_report_file} \\
        --variable data_table="$table" \\
        --to markdown > {report_file}
                """.format(
                                        data_dir = os.path.join(self.REP_DIR,report_data),
                                        data_files = data_files,
                                        report_file = report_file,
                                        report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file)
                                ),
                                        report_files = [report_file],
                                        name = 'region_enrichment_report'
                                        )) 

		return jobs
	
	def apply_models(self):
		"""
		This step applies models provided by the user on each sample and scores them.
		The result is a csv file containing raw and probability scores for each sample-model pair.
		"""
		jobs=[]
		if config.param('apply_models','apply_models') == 'y':
			design_csv = csv.DictReader(open(config.param('apply_models','design_file')), delimiter='\t')
			diseases = []
			cpg_files = []
			models = []
			output_file = os.path.join(self.REP_DIR, 'data/models/scores.csv')
			rgset_file = os.path.join(self.GR_DIR, "grset_norm.rds")
			for line in design_csv:
				diseases.append(line['disease'])
				cpg_files.append(line['cpgs_file'])
				models.append(line['model_file'])
			jobs.append(Job(
					[rgset_file],
					[output_file],
					command="""\
			mkdir -p report/data/models
			R --vanilla <<-'EOF'
			suppressPackageStartupMessages(library(minfi))
			suppressPackageStartupMessages(library(caret))
			rgset <- readRDS("{rgset_file}")
			pdata <- pData(rgset)
			beta <- getBeta(rgset)
			beta_rows <- rownames(beta)
			colnames(beta) <- pdata$Sample_Name
			diseases <- c{diseases}
			cpg_files <- c{cpg_files}
			output <- data.frame(row.names = colnames(beta),Age = pdata$Age, Sex = pdata$Sex, Group = pdata$Sample_Group)
			models <- c{models}
			for (i in 1:length(diseases)){{
				cpgs <- read.csv(cpg_files[i],header=FALSE,row.names=1, as.is=TRUE)
				dim(cpgs)
				rows <- match(rownames(cpgs), beta_rows)
				B <- beta[rows,]
				dim(B)
				model <- load(models[i])
				raw_scores <- predict(fit.model, newdata = t(B), type="raw")
				prob_scores <- predict(fit.model, newdata = t(B), type="prob")
				output[,paste(diseases[i],"_raw",sep="")] <- raw_scores
				output[,paste(diseases[i],"_prob",sep="")] <- prob_scores
			}}
			write.csv(output, '{output_file}')
			EOF
			""".format(
					diseases = tuple(diseases) if len(diseases) > 1 else "('"+diseases[0]+"')",
					cpg_files = tuple(cpg_files) if len(cpg_files) > 1 else "('"+cpg_files[0]+"')",
					models = tuple(models) if len(models) > 1 else "('"+models[0]+"')",
					rgset_file = rgset_file,
					output_file = output_file),
					
					name = "apply_models"))

			#generate report tables
			report_file = 'report/Epigene.model_scores.md'
			jobs.append(Job(
							[output_file],
							[report_file],
							[['apply_models', 'module_pandoc']],
							command="""\
			cat {scores_file} | sed 's/,/\t/g' | cut -f1-6 | sed 's/\"//g' > report/scores.tsv && \\
			scores_table_md=`head -7 report/scores.tsv | LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print "Sample Name|Age|Sex|Group",$5,$6; print "-----|-----:|-----:|-----:|-----:|-----:|"}} else {{print $1,$2,$3,$4,$5, sprintf("%.4f", $6)}}}}' ;` && \\
			pandoc \\
                {report_template_dir}/{basename_report_file} \\
                --template {report_template_dir}/{basename_report_file} \\
                --variable scores_table="$scores_table_md" \\
                --to markdown > {report_file}
				""".format(
							scores_file = output_file,
							report_template_dir = self.report_template_dir,
							basename_report_file = os.path.basename(report_file),
							report_file = report_file
					),
							report_files = [report_file],
							name = 'apply_models_report'))
	
		return jobs		
			
				
			


	def build_model(self):
		"""
		This step uses the R caret package to build a model for the data provided.
		The output is an .RData file with the object fit.model that can later be used to score other samples.
		"""
		jobs = []
		if config.param('build_model','build_model') == 'y':
			rgset_file = os.path.join(self.GR_DIR, 'grset_norm.rds')
			method = config.param('build_model','method') 
			for contrast in self.contrasts:
				if config.param('differential_methylated_pos','minfi') == 'y':
					dmps_file = os.path.join(self.CPG_DIR, contrast.name+'_minfi_dmps.csv')
					jobs.append(model.build_model(rgset_file, dmps_file, contrast, method, 'minfi'))
				if config.param('differential_methylated_pos','limma') == 'y':
					dmps_file = os.path.join(self.CPG_DIR, contrast.name+'_limma_dmps.csv')
					jobs.append(model.build_model(rgset_file, dmps_file, contrast, method,'limma'))
		return jobs
	
	
if __name__ == '__main__':
	Epigene()


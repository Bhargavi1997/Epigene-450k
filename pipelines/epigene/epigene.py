#!/usr/bin/env python


#Epignetic pipeline for Illumina 450k

#Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys

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
	The final output is a list CpG sites and their methylation levels, a list of differentially methylated positions and a list of 
	differentially methylated regions
	"""

	RG_DIR = "rgchannel_set"
	GR_DIR = "gr_set"
	QC_DIR = "quality_control"
	CPG_DIR = "cpg"

	def __init__(self):
		self.argparser.add_argument("-r", "--samples", help="sample file", type=file)
		super(Epigene, self).__init__()
		
	@property
	def samples(self):
		if not hasattr(self, "_samples"):
			if self.args.samples:
				self._samples = parse_illumina450k_sample_file(self.args.samples.name)
			else:
				self.argparser.error("arguement -r/--samples is required!")

		return self._samples 

	def read_idat(self):	
		"""
		Read raw .idat files together with the experiment sample sheet to generate an annotated 
		minfi RGChannelSet object
		"""
		basenames = []
		
		for sample in self.samples:
			basenames.append(sample.basename)

		job = Job(
					[self.args.samples.name],
					[os.path.join(self.RG_DIR,"rgset.rds")],
					[
						["read_idat","module_R"],
						["read_idat","module_mugqic_R_packages"]
					],
					command="""\
	mkdir -p {rgchannel_dir}
	R --no-save --no-restore <<-EOF
	suppressPackageStartupMessages(library(minfi))
	targets <- read.csv("{samplesheet}")
	targets['Basename'] <- c{basenames} 
	rgset <- read.450k.exp(target = targets)
	saveRDS(rgset, file = "{output_file}")

	EOF""".format(
			rgchannel_dir=self.RG_DIR,
			#samplesheet=config.param("read_idat", "samplesheet", type = "filepath"),
			samplesheet=self.args.samples.name,
			basenames=tuple(basenames),
			output_file=os.path.join(self.RG_DIR,"rgset.rds")
			),
			name="read_idat")

		return [job]
	
	def cell_counts(self):
		"""
		Estimate confounding levels between phenotype and cell type composition. The cell type
		composition of each sample is returned. This step should only be included for blood samples 
		"""
		for sample in self.samples:
			if (sample.cell_type).lower() != "blood":
				log.info("Warning: Sample " + sample.name + "is not a blood sample!")

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
	rgset\$Sex <- as.character(rgset\$Sex)
	cellCounts <- estimateCellCounts(rgset)
	write.csv(as.data.frame(cellCounts), file="{output_file}")

	EOF""".format(
			qc_dir=self.QC_DIR,
			rgset_file=os.path.join(self.RG_DIR,"rgset.rds"),
			output_file=os.path.join(self.QC_DIR,"cell_counts.csv")
			),
			name="cell_counts")

		return [job]

	def normalize(self):

		"""
		Normalize the data from a selection of various normalization methods
		"""						
			
		if config.param("normalize","method").lower() == "funnorm":
			job = concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.func(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.abspath(os.path.join(self.GR_DIR,"grset_norm.rds"))
				)
			])
		elif config.param("normalize","method").lower() == "swan":
			job = concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.swan(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			])
		elif config.param("normalize","method").lower() == "noob":
			job = concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.noob(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			])
		elif config.param("normalize","method").lower() == "raw":
			job = concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.raw(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			])
		elif config.param("normalize","method").lower() == "quantile":
			job = concat_jobs([
				Job(command="mkdir -p " + self.GR_DIR),
				norm.quantile(
					os.path.join(self.RG_DIR,"rgset.rds"),
					os.path.join(self.GR_DIR,"grset_norm.rds")
				)
			])
		job.name = "normalize"

		return [job]


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
							qc.find_failures(
								os.path.join(self.RG_DIR,"rgset.rds"),
								os.path.join(self.QC_DIR,"failed_probes.csv"),
								config.param("find_failures","padjust_method"),
								config.param("find_failures","pvalue"),
								config.param("find_failures","percent_failed_threshold"))
							], name="finding_failures")
			)
		if config.param("probe_filter","find_extremes").lower() == "y":
			jobs.append(concat_jobs([
						Job(command="mkdir -p " + self.QC_DIR),
						qc.find_extremes(
							os.path.join(self.RG_DIR,"rgset.rds"),
							os.path.join(self.QC_DIR,"failed_proves.csv"))
						], name="find_extremes")
			)
		if config.param("probe_filter","find_chromosomes").lower() == "y":
			jobs.append(concat_jobs([
						Job(command="mkdir -p " + self.QC_DIR),
						qc.find_chromosomes(
							os.path.join(self.RG_DIR,"rgset.rds"),
							os.path.join(self.QC_DIR,"chr_probes.csv"),
							config.param("find_chromosomes","chromosomes",type = "list"))
						], name="find_chromosomes")
			)
		if config.param("probe_filter","find_snps").lower() == "y":
			jobs.append(concat_jobs([
						Job(command="mkdir -p " + self.QC_DIR),
						qc.find_snps(
							os.path.join(self.RG_DIR,"rgset.rds"),
							os.path.join(self.QC_DIR,"snp_probes.csv"))
						], name="find_chromosomes")
			)
		
		#Remove chosen probes from normalized data		
		jobs.append(concat_jobs([
			Job(command="mkdir -p " + self.QC_DIR),
			qc.merge_remove(
				os.path.join(self.GR_DIR,"grset_norm.rds"),
				os.path.join(self.GR_DIR,"grset_norm_qc.rds"),
				os.path.abspath(self.QC_DIR),
			)
		], name="probe_filter"))

						
		return jobs

	def differential_methylated_pos(self):
		"""
		This step finds CpG sites that are differentially methylated by sample group (i.e control vs case)
		Essentially, this a wraparound the minfi function dmpFinder and uses a F-test to generate p-values
		for each site. The output is a table of differentially methylated CpG sites with associated information
		such as chromsomal position, average methylation difference (from controls), p-value, and whether it is
		confounded by blood composition (for blood samples only) 
		"""
		job = concat_jobs([
			Job(command="mkdir -p " + self.CPG_DIR), 
			dmp.differential_methylated_pos(
				os.path.join(self.GR_DIR,"grset_norm_qc.rds"),
				os.path.join(self.CPG_DIR,"dmps.csv"),
				config.param("differential_methylated_pos","padjust_method"),
				config.param("differential_methylated_pos","pvalue"),
				config.param("differential_methylated_pos","delta_beta_threshold"))
			])	

		job.name="differential_methylated_pos"

		return [job]
	
	def differential_methylated_regions(self):
		"""
		This step finds differentially methalated regions using the bumphunting algorithm
		"""
		job = concat_jobs([
			Job(command="mkdir -p " + self.CPG_DIR),
			dmp.find_dmr(
				os.path.join(self.GR_DIR,"grset_norm_qc.rds"),
				os.path.join(self.CPG_DIR,"dmrs.csv"),
				config.param("differential_methylated_regions","permutations"),
				config.param("differential_methylated_regions","delta_beta_threshold"),
				config.param("differential_methylated_regions","cores"))
			])

		job.name="differential_methylated_regions"

		return [job]
		
	def principal_comp(self):
		"""
		This step generals a principal component analysis plot of the first to PC in the data. Samples are
		differentiated by their group (i.e controls or cases)
		"""
		job = Job(
					[os.path.join(self.GR_DIR,"grset_norm_qc.rds")],
					["pcaplot.pdf"],
					[
						['principal_comp','module_R'],
						['principal_comp','module_mugqic_R_packages']
					],
					command="""\
	R --no-save --no-restore <<-EOF
	suppressPackageStartupMessages(library(minfi))
	library(ggplot2)
	
	gset <- readRDS("{input_rds_file}")
	B <- getBeta(gset)
	pca <- prcomp(t(B), scale = TRUE)
	pdf("{output_file}")
	df <- as.data.frame(pca\$x)
	df["Groups"] <- pData(gset)[,"Sample_Group"]
	qplot(x=PC1, y=PC2, data=df, color=Groups)
	ggsave("pcaplot.pdf", plot=last_plot())
	EOF""".format(
			input_rds_file=os.path.join(self.GR_DIR,"grset_norm_qc.rds"),
			output_file="pcaplot.pdf"
			),
			name="principal_comp")

		return [job]
	
	@property
	def steps(self):
		return [
				self.read_idat,
				self.cell_counts,
				self.normalize,
				self.probe_filter,
				self.differential_methylated_pos,
				self.differential_methylated_regions,
				self.principal_comp
				]	
	
if __name__ == '__main__':
	Epigene()

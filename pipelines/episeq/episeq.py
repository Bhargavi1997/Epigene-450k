#!/usr/bin/env python


#Epignetic pipeline for RRBS sequential data

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
from pipelines import common
import utils

log = logging.getLogger(__name__)

class Episeq(common.Illumina):

	def trim_galore(self):	

		#Trim FASTQ files in readset file using Trim Galore
 		jobs = []
		for readset in self.readsets:
			trim_directory = os.path.join("trimmed", readset.sample.name)
			trimmed_fq = os.path.join(trim_directory, readset.name)
			jobs.append(
				Job(
					[readset.fastq1, readset.fastq2],
					[trimmed_fq + "_R1_val_1.fq.gz", trimmed_fq + "_R2_val_2.fq.gz"],
					#["trim_galore","module_trim_galore"],
					command="""\
	mkdir -p {directory}
	module load trim_galore/0.4.1
	trim_galore --rrbs --non_directional --paired --output_dir {directory} {fastq1} {fastq2}
	""".format(
				directory=trim_directory,
				fastq1=readset.fastq1,
				fastq2=readset.fastq2			
				),
				name="trim_galore." + readset.name)	
			)	
		
		return jobs
	
	
	def bismark_align(self):

		#Align trimmed reads to reference genome using Bismark
		#Multicore and basname not currently supported together
		jobs = []
		for readset in self.readsets:
			trim_prefix = os.path.join("trimmed", readset.sample.name, readset.name)
			align_directory = os.path.join("aligned", readset.sample.name)
			readset_sam = os.path.join(align_directory, readset.name + "_aligned_pe.sam.gz")
			
			job = Job(
					[trim_prefix + "_R1_val_1.fq.gz", trim_prefix + "_R2_val_2.fq.gz"],
					[readset_sam],
					[["bismark_align","module_bowtie2"]],
					command="""\
	mkdir -p {directory}
	module load bismark/0.15
	bismark -q --non_directional --output_dir {directory} --basename {basename} {bs_refgene} -1 {fastq1} -2 {fastq2}
	""".format(
				directory=align_directory,
				fastq1=trim_prefix + "_R1_val_1.fq.gz",
				fastq2=trim_prefix + "_R2_val_2.fq.gz",
				basename=readset.name + '_aligned',
				bs_refgene=config.param("bismark_align", "bs_refgene", type="dirpath")
				),
				name="bismark_align." + readset.name)
	
			jobs.append(job)

		return jobs	

	def merge_sam_files(self):
		
		#Merge SAM files belonging to the same sample 
		jobs = []
		for sample in self.samples:
			align_directory = os.path.join("aligned", sample.name)
			readset_sams = [os.path.join(align_directory, readset.name + "_aligned_pe.sam.gz") for readset in sample.readsets]
			
			#job=picard.merge_sam_files(readset_sams, os.path.join(align_directory, sample.name + "_merged_aligned_pe.sam.gz"), sort_order="unsorted")
			#job.name="picard_merge_sam_files." + sample.nameb

			job = Job(
					readset_sams,
					[os.path.join(align_directory, sample.name + "_merged_aligned_pe.sam.gz")],
					command="""\
			cat {files} > {merged_file}
			""".format(
						files=" ".join(readset_sams),
						merged_file=os.path.join(align_directory, sample.name + "_merged_aligned_pe.sam.gz")
						),
						name="merge_sam_files." + sample.name)

			jobs.append(job)

		return jobs

	def bismark_methylation_caller(self):
		
		#Call methylated positions from all samples using Bismark
		aligned_samples = []
		for sample in self.samples:
			aligned_samples.append(os.path.join("aligned", sample.name, sample.name + "_merged_aligned_pe.sam.gz"))
		
		job = Job(
				aligned_samples,
				[os.path.join("methyl_calls", sample.name + "_merged_aligned_pe.sam.bismark.cov.gz") for sample in self.samples],
				command="""\
		mkdir -p {directory}
		module load bismark/0.15
		bismark_methylation_extractor -p --no_overlap --output {directory} --multicore {cores} --bedGraph {sample}
		""".format(
					directory="methyl_calls",
					sample=" ".join(aligned_samples),
					cores=config.param("bismark_methylation_caller","cores",type="int")/3
					),
					name="bismark_methylation_caller")
				
		return [job]

	def find_dmps(self):
		#Find list of differentially methylated CpGs

		methyl_files = [os.path.join("methyl_calls", sample.name + "_merged_aligned_pe.sam.bismark.cov.gz") for sample in self.samples]
		dmps_file = os.path.join("dmp", "RRBS_dmps.csv")
		job = Job(
				methyl_files,
				[dmps_file],
				[
					["find_dmps","module_R"],
					["find_dmps","module_mugqic_R_packages"]
				],
				command="""\
		mkdir -p {directory}
		R --no-save --no-restore <<-EOF
		suppressPackageStartupMessages(library(bumphunter))
		suppressPackageStartupMessages(library(BiSeq))
		suppressPackageStartupMessages(library(minfi))
		library(doParallel)
		registerDoParallel(cores={cores})
		samples <- c{samples}
		rrbs <- readBismark(samples, colData=DataFrame(group=factor(c(rep("control",3), rep("disease",3))), row.names=c{sample_names}))
		
		rrbs.filtered <- rrbs[apply(totalReads(rrbs), 1, function(x) all(x > 10)),]
		beta <- methLevel(rawToRel(rrbs.filtered))
		chr <- as.character(seqnames(colData(rrbs.filtered))
		dmp <- dmpFinder(beta, pheno=colData(rrbs.rel)[,"group"], type="categorical")
 		dmp["pval"] <- p.adjust(dmp[,"pval"], method = "bonferroni")
 		dmp <- dmp[dmp["pval"] < 0.05,]["pval"]
 		
		write.csv(dmp, file="{dmps_file}", quote=FALSE)

		EOF""".format(
						directory=os.path.dirname(dmps_file),
						samples=tuple(methyl_files),
					  	sample_names=tuple([sample.name for sample in self.samples]),
					 	cores=config.param("find_dmps","cores", type="int"),
						dmps_file=dmps_file
						),
						name="find_dmps")

		return [job]
	
	def find_dmrs(self):
		
		#Find and define DMR regions for each sample using Bioconductor package BiSeq
		methyl_files = [os.path.join("methyl_calls", sample.name + "_merged_aligned_pe.sam.bismark.cov.gz") for sample in self.samples]
		dmrs_file = os.path.join("dmrs", "dmrs.csv")
		job = Job(
				methyl_files,
				[dmrs_file],
				[
					["find_dmrs","module_R"],
					["find_dmrs","module_mugqic_R_packages"]
				],
				command="""\
		mkdir -p {directory}
		R --no-save --no-restore <<-EOF
		suppressPackageStartupMessages(library(bumphunter))
		suppressPackageStartupMessages(library(BiSeq))
		library(doParallel)
		registerDoParallel(cores={cores})
		samples <- c{samples}
		rrbs <- readBismark(samples, colData=DataFrame(group=factor(c(rep("control",3), rep("disease",3))), row.names=c{sample_names}))
		rrbs.filtered <- rrbs[apply(totalReads(rrbs), 1, function(x) !all(x < 10)),]
	
		beta <- methLevel(rawToRel(rrbs.filtered))
		chr <- as.character(seqnames(rowData(rrbs.rel)))
		pos <- start(ranges(rowData(rrbs.rel)))
		pheno <- colData(rrbs.rel)[,"group"]
		designM <- model.matrix(~pheno)

		dmrs <- bumphunterEngine(methLevel(rrbs.rel),
								 chr=chr,
								 pos=pos,
								 design=designM,
								 cutoff=0.8,
								 pickCutoffQ=0.99,
								 null_method=c("permutation","bootstrap"),
								 smooth=FALSE,
								 smoothFunction=locfitByCluster,
								 B=1000,
								 verbose=TRUE,
								 maxGap=500)

		write.csv(dmrs\$table, "{dmrs_file}", quote=FALSE)

		EOF""".format(
						directory=os.path.dirname(dmrs_file),
						samples=tuple(methyl_files),
					  	sample_names=tuple([sample.name for sample in self.samples]),
					 	cores=config.param("find_dmrs","cores", type="int"),
						dmrs_file=dmrs_file
						),
						name="find_dmrs")

		return [job]
	

	@property
	def steps(self):
		return [
			self.trim_galore,
			self.bismark_align,
			#self.picard_sort_sam_files,
			self.merge_sam_files,
			self.bismark_methylation_caller,
			self.find_dmps,
			self.find_dmrs
		]	
	
if __name__ == '__main__':
	Episeq()

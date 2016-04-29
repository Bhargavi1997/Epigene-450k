#!/usr/bin/env python 
import os 
from core.config import *
from core.job import * 

def find_failures(rds_file, output_file, padjust_method, pvalue, threshold ):
	return Job(
		[rds_file],
        [output_file],
        [
        	['qc_find_failures', 'module_R'],
            ['qc_find_failures', 'module_mugqic_R_packages']
        ],
        command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi)) 
rgset <- readRDS("{rds_file}")
detP <- detectionP(rgset)
detP <- apply(detP, 2, function(x) p.adjust(x, method = "{padjust_method}"))
fails <- detP > {pvalue}
failedProbes <- rownames(fails)[rowMeans(fails) > {threshold}]
write(failedProbes, "{output_file}", sep=',')
q(save="no")
  
EOF""".format(
		rds_file=rds_file,
		output_file=output_file,
		padjust_method=padjust_method,
      	pvalue=pvalue,
        threshold=threshold
    ))

def find_extremes(rds_file, output_file):
	return Job(
		[rds_file],
		[output_file],
		[
			['qc_find_extremes', 'module_R'],
			['qc_find_extremes', 'module_mugqic_R_packages']
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
rgset <- readRDS("{rds_file}")
Mset <- preprocessRaw(rgset)
extremeBetas <- apply(getBeta(Mset), 1, function(x) any(x==0 | x==1))
extremeProbes <- names(which(extremeBetas))
write(extremeProbes, "{output_file}", sep=',')

EOF""".format(
		rds_file=rds_file,
		output_file=output_file
	))

def find_chromosomes(rds_file, output_file, chromosomes):
	return Job(
		[rds_file],
		[output_file],
		[
			['qc_find_chromosomes', 'module_R'],
			['qc_find_chromosomes', 'module_mugqic_R_packages']
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
rgset <- readRDS("{rds_file}")
anno <- getAnnotation(rgset)
chrProbes <- rownames(anno[anno[,'chr'] %in% c{chromosomes},])
write(chrProbes, "{output_file}", sep=',')

EOF""".format(
		rds_file=rds_file,
		output_file=output_file,
		chromosomes=tuple(chromosomes)
	))

def find_snps(rds_file, output_file):
	return Job(
		[rds_file],
		[output_file],
		[	
			["qc_find_snps", "module_R"],
			["qc_find_snps", "module_mugqic_R_packages"]
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
rgset <- readRDS("{rds_file}")
snp <- getSnpInfo(rgset)
snpProbes <- rownames(snp[!(is.na(snp[,"CpG_rs"]) & is.na(snp[,"SBE_rs"])),])
write(snpProbes, "{output_file}", sep=',')

EOF""".format(
		rds_file=rds_file,
		output_file=output_file
	))

def merge_remove(input_rds_file, output_rds_file, path):

	return Job(
		[input_rds_file],
		[output_rds_file],
		[
			['qc_merge_remove', 'module_R'],
			['qc_merge_remove', 'module_mugqic_R_packages']
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
grset <- readRDS("{input_rds_file}")

probe_list <- dir(path="{path}", pattern="*.csv")
removedProbes <- lapply(probe_list, scan, what=character(), sep=',')
removedProbes <- unique(unlist(removedProbes))

grset <- grset[!(rownames(grset) %in% removedProbes)]
saveRDS(grset, "{output_rds_file}")

EOF""".format(
		input_rds_file=input_rds_file,
		output_rds_file=output_rds_file,
		path=path
	))


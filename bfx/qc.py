#!/usr/bin/env python 
import os 
from core.config import *
from core.job import * 


def quality_control_plots(input_rds_file, output_file,path,groupname):
	return Job(
		[input_rds_file],
		[output_file],
		[
			['qc_quality_control_plots', 'module_R'],
			['qc_quality_control_plots', 'module_mugqic_R_packages']
		],
		command="""\
mkdir -p report && \\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi)) 
library(beanplot)

rgset <- readRDS("{input_rds_file}")
Mset <- preprocessRaw(rgset)
qc <- getQC(Mset)
Mset <- addQC(Mset, qc=qc)
png("{output_file}")
plotQC(qc)
dev.off()

png(file = "{path}/densityPlot.png", units = "in", width = 8.5, height = 8.5, res = 300)
groups <- pData(rgset)\${group}
par(cex.axis=0.7)
densityPlot(rgset, sampGroups = groups, main = "Pre-normalization density plot", legend = TRUE)
dev.off()

png(file = "{path}/densityBeanPlot.png", units = "in", width = 8.5, height = 8.5, res = 300)
names <- pData(rgset)\$Sample_Name
groups <- pData(rgset)\${group}
par(cex.axis=0.7)
densityBeanPlot(rgset, sampName = names, sampGroups = groups, main = "Pre-normalization density bean plot")
dev.off()

pdf("report/qcReport.pdf")
names <- pData(rgset)\$Sample_ID
groups <- pData(rgset)\${group}
qcReport(rgset, sampName = names, sampGroups = groups, pdf="report/qcReport.pdf")
dev.off()

EOF""".format(
		input_rds_file=input_rds_file,
		output_file=output_file,
		group=groupname,
		path = path
	))

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


def cross_reactive_probes(input_file, output_file):
	return Job(
		[input_file],
		[output_file],
		command="""\
cp {input_file} {output_file}
""".format(
		input_file=input_file,
		output_file=output_file
	))

def find_snps_chen(input_file, output_file, minimal_distance_to_SNPs, minimal_minor_allele_frequency):
	return Job(
		[input_file],
		[output_file],
		command="""\
awk -F, '$2 >= {minimal_distance_to_SNPs} && $3 >= {minimal_minor_allele_frequency}' {input_file} > report/data/quality_control/temp.csv
cat report/data/quality_control/temp.csv | cut -d, -f1 > {output_file}
rm report/data/quality_control/temp.csv
""".format(
		input_file=input_file,
		output_file=output_file,
		minimal_distance_to_SNPs=minimal_distance_to_SNPs,
		minimal_minor_allele_frequency=minimal_minor_allele_frequency
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

def find_snps_minfi(rds_file, output_file):
	return Job(
		[rds_file],
		[output_file],
		[	
			["qc_find_snps_minfi", "module_R"],
			["qc_find_snps_minfi", "module_mugqic_R_packages"]
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
cat {path}/*.csv > {path}/total.csv
sort {path}/total.csv | uniq > {path}/uniq.csv 
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
grset <- readRDS("{input_rds_file}")

probe_list <- read.csv("{path}/uniq.csv", sep=",")
removedProbes <- unlist(probe_list)

grset <- grset[!(rownames(grset) %in% removedProbes)]
saveRDS(grset, "{output_rds_file}")

EOF""".format(
		input_rds_file=input_rds_file,
		output_rds_file=output_rds_file,
		path=path
	))


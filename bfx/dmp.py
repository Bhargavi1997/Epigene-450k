#!/usr/bin/env python 
import os
from core.config import *
from core.job import *

def differential_methylated_pos(input_rds_file, output_file, padjust, pvalue, dnam_threshold):
	return Job(
		[input_rds_file],
		[output_file],
		[
			['differential_methylated_pos','module_R'],
			['differential_methylated_pos','module_mugqic_R_packages']
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
library(FlowSorted.Blood.450k)

grset <- readRDS("{input_rds_file}")
beta <- getBeta(grset)
#dmpFinder does not work if the M value matrix contins zeros or Inf
beta[beta == 0] <- 0.0001
beta[beta == 1] <- 1.0001
M <- log2(beta/(1-beta))
P <- pData(grset)

dmp <- dmpFinder(M, pheno=P[,"Sample_Group"], type="categorical")
dmp["pval"] <- p.adjust(dmp[,"pval"], method = "{padjust}")
dmp <- dmp[dmp["pval"] < {pvalue},]["pval"]

#Find which samples are controls and which or not in order to calculate methylation differences
control <- rownames(P[tolower(P[,'Sample_Group']) == 'control',])
disease <- rownames(P[tolower(P[,'Sample_Group']) != 'control',])

avgDiseaseBeta = rowMeans(beta[,disease])
avgControlBeta = rowMeans(beta[,control])
avgDeltaBeta = avgDiseaseBeta - avgControlBeta
chrs <- getAnnotation(grset)[,'chr']
genes <- getAnnotation(grset)[,'UCSC_RefGene_Name']

result <- data.frame(avgDiseaseBeta, avgControlBeta, avgDeltaBeta, chrs, genes)
rownames(result) <- rownames(beta)
result <- merge(result, dmp, by=0)
names(result) <- c("CpG Site", "Avg Case Beta", "Avg Control Beta", "Avg Delta Beta", "Chromosome", "Gene", "P-value")
result <- result[abs(result["Avg Delta Beta"]) > {dnam_threshold},]

blood <- FlowSorted.Blood.450k.compTable
blood["p.value"] <- p.adjust(blood[,"p.value"], method="fdr")
result['Cell Affected'] <- result[,'CpG Site'] %in% rownames(blood)[blood['p.value'] < {pvalue}]

write.csv(result, "{output_file}")

EOF""".format(
		input_rds_file=input_rds_file,
		output_file=output_file,
		padjust=padjust,
		pvalue=pvalue,
		dnam_threshold=dnam_threshold
	))

def find_dmr(input_rds_file, output_file, permutations, cutoff, cores):
	return Job(
		[input_rds_file],
		[output_file],
		[
			["find_dmr","module_R"],
			["find_dmr","module_mugqic_R_packages"]
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(doParallel))
registerDoParallel(cores={cores})

grset <- readRDS("{input_rds_file}")
pheno <- pData(grset)[,'Sample_Group']
designMatrix <- model.matrix(~ pheno)

dmrs <- bumphunter(grset, design = designMatrix, cutoff={cutoff}, B={permutations}, type="Beta")

write.csv(dmrs\$table, "{output_file}")


EOF""".format(
		input_rds_file=input_rds_file,
		output_file=output_file,
		permutations=permutations,
		cutoff=cutoff,
		cores=cores
	))

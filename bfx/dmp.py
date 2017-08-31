#!/usr/bin/env python 
import os
from core.config import *
from core.job import *

def minfi_differential_methylated_pos(contrast,grset_file,input_rds_file, output_file, padjust, pvalue, dnam_threshold):
	contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
        sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
	
	return Job(
		[input_rds_file],
		[output_file],
		[
			['minfi_differential_methylated_pos','module_R'],
			['minfi_differential_methylated_pos','module_mugqic_R_packages']
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))

grset <- readRDS("{grset_file}")

beta <- getBeta(grset)
#dmpFinder does not work if the M value matrix contins zeros or Inf
beta[beta == 0] <- 0.0001
beta[beta == 1] <- 1.0001

M <- getM(readRDS("{input_rds_file}"))
colnames(M) <- pData(grset)\$Sample_Name
P <- DataFrame(group=factor(c{group}), row.names=c{sample_names})
samprows <- rownames(P)
dmp <- dmpFinder(M[,c(samprows)], pheno=P\$group, type="categorical")
dmp["adj.pval"] <- p.adjust(dmp[,"pval"], method = "{padjust}")
dmp <- dmp[dmp["adj.pval"] < {pvalue},][c("pval","adj.pval")]

controls <- c({controls})
cases <- c({cases})

colnames(beta) = pData(grset)\$Sample_Name

avgControlBeta = rowMeans(beta[,controls])
avgCaseBeta = rowMeans(beta[,cases])
avgDeltaBeta = avgCaseBeta - avgControlBeta
chrs <- getAnnotation(grset)[,'chr']
genes <- getAnnotation(grset)[,'UCSC_RefGene_Name']
result <- data.frame(avgCaseBeta, avgControlBeta, avgDeltaBeta, chrs, genes)
result\$genes = strsplit(as.character(result\$genes), ";")
result\$genes = as.character(lapply(result\$genes, function(x) {{paste(unique(as.character(x)), collapse=";")}}))
rownames(result) <- rownames(beta)
result <- merge(result, dmp, by=0)
names(result) <- c("CpG Site", "Avg Case Beta", "Avg Control Beta", "Avg Delta Beta", "Chromosome", "Gene", "P-value","Adj.P-value")
result <- result[abs(result["Avg Delta Beta"]) > {delta_beta_threshold},]

write.csv(result, "{dmps_file}", row.names=FALSE)

EOF""".format(
		group=tuple(sample_group),
                sample_names=tuple([sample.name for sample in contrast_samples]),
                controls=', '.join(["'" + sample.name + "'" for sample in contrast.controls]),
                cases=', '.join(["'" + sample.name + "'" for sample in contrast.treatments]),
		grset_file=grset_file,
		input_rds_file=input_rds_file,
		dmps_file=output_file,
		padjust=padjust,
		pvalue=pvalue,
		delta_beta_threshold=dnam_threshold
	))

def limma_differential_methylated_pos(contrast, grset_file,input_rds_file, output_file, adjust_method, pval, delta_beta_threshold, design_matrix_columns):
	
	contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
        sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
	
	return Job(
		[input_rds_file],
		[output_file],
		[
			['minfi_differential_methylated_pos','module_R'],
                        ['minfi_differential_methylated_pos','module_mugqic_R_packages']
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(minfi))
M <- getM(readRDS("{input_rds_file}"))
grset <- readRDS("{grset_file}")
beta <- getBeta(grset)
P <- DataFrame(group=factor(c{group}), row.names=c{sample_names})
samprows <- rownames(P)
pdata <- pData(grset)
colnames(grset) <- pdata\$Sample_Name
Age <- as.factor(pdata\$Age)
Sex <- as.factor(pdata\$Sex)
design <- model.matrix(~P\$group {design_matrix_columns}, data=pdata)
design
colnames(M) <- colnames(grset)
fit <- lmFit(M[,c(samprows)], design)
fit2 <- eBayes(fit)
probe <- topTable(fit2,coef=2,adjust="{adjust_method}",num=Inf)

ann <- getAnnotation(grset)
probe\$genes <- ann[match(rownames(probe), ann\$Name),]\$UCSC_RefGene_Name
probe\$chromosome <- as.character(ann[match(rownames(probe), ann\$Name),]\$chr)
probe\$genes = strsplit(as.character(probe\$genes), ";")
probe\$genes = as.character(lapply(probe\$genes, function(x) {{paste(unique(as.character(x)), collapse=";")}}))
controls <- c({controls})
cases <- c({cases})

colnames(beta) = pData(grset)\$Sample_Name
avgControlBeta = rowMeans(beta[,controls])
avgCaseBeta = rowMeans(beta[,cases])
avgDeltaBeta = data.frame(avgCaseBeta - avgControlBeta)
avgDeltaBeta\$Name = rownames(beta)


probe\$avgDeltaBeta <- avgDeltaBeta[match(rownames(probe), avgDeltaBeta\$Name),1]
probe <- probe[which(probe\$adj.P.Val <= {pval}),]

probe <- probe[which(abs(probe\$avgDeltaBeta) > {delta_beta_threshold}),]
#na.omit(probe)
probe <- cbind("CpG Site" = rownames(probe), probe)
write.csv(probe, "{output_file}", row.names = FALSE)
EOF""".format(
		group=tuple(sample_group),
                sample_names=tuple([sample.name for sample in contrast_samples]),
                controls=', '.join(["'" + sample.name + "'" for sample in contrast.controls]),
                cases=', '.join(["'" + sample.name + "'" for sample in contrast.treatments]),
		grset_file=grset_file,
		design_matrix_columns=design_matrix_columns,
		input_rds_file = input_rds_file,
		output_file = output_file,
		adjust_method = adjust_method,
		pval = pval,
		delta_beta_threshold = delta_beta_threshold,
		contrast_name = contrast.name
))
		

def find_dmr(contrast, input_rds_file, output_file, nullMethod, permutations, cutoff,pvalue, LCutoff, cores):
	contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
	sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
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
colnames(grset) <- pData(grset)\$Sample_Name
P <- DataFrame(group=factor(c{group}), row.names=c{sample_names})
samprows <- rownames(P)
designMatrix <- model.matrix(~ P\$group)

dmrs <- bumphunter(grset[,c(samprows)], design = designMatrix, cutoff={cutoff}, nullMethod="{nullMethod}", B={permutations}, type="Beta")
dmrs <- dmrs\$table
dmrs <- dmrs[which(dmrs\$p.value <= {pvalue}),]
dmrs <- dmrs[which(dmrs\$L >= {LCutoff}),]
ann <- getAnnotation(grset)
pos <- mapply(function(x,y) {{c(x:y)}}, dmrs\$start, dmrs\$end)
dmrs\$genes <- mapply(function(x,y) {{ann[which(ann\$pos %in% data.frame(x)[,1] & ann\$chr == y),]\$UCSC_RefGene_Name}}, pos,dmrs\$chr)
dmrs\$genes <- sapply(dmrs\$genes, unlist)
dmrs\$genes <- sapply(dmrs\$genes, function(x){{paste(x[x!=""],sep=";",collapse=";")}})
dmrs\$genes <- sapply(dmrs\$genes, function(x){{strsplit(x, split=";")}})
dmrs\$genes <- as.character(sapply(dmrs\$genes, function(x){{paste(unique(x),collapse=";")}}))
write.csv(dmrs , "{output_file}")
EOF""".format(
		group=tuple(sample_group),
                sample_names=tuple([sample.name for sample in contrast_samples]),
		input_rds_file=input_rds_file,
		output_file=output_file,
		nullMethod=nullMethod,
		permutations=permutations,
		cutoff=cutoff,
		LCutoff = LCutoff,
		pvalue = pvalue,
		cores=cores
	))


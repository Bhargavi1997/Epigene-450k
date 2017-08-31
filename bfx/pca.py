#!/usr/bin/env python 
import os
from core.config import *
from core.job import *


def pca_plot(input_rds_file,input_file, output_file, contrast, jobname):
	contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
        sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
	return Job(
		[input_file, input_rds_file],
		[output_file + "-1.png"],
		[
			['principal_comp','module_R'],
			['principal_comp','module_mugqic_R_packages']
		],
		command="""\
	R --no-save --no-restore <<-EOF
	suppressPackageStartupMessages(library(minfi))
	library(ggplot2)
	extension <- tools::file_ext("{input_file}")
	if(extension == 'rds'){{
		input <- readRDS("{input_rds_file}")
	}}else{{
		input <- read.csv("{input_file}")
	}}
	gset <- readRDS("{input_rds_file}")
	P <- DataFrame(group=factor(c{group}), row.names=c{sample_names})
	samprows <- rownames(P)
	pdata <- pData(gset)
	colnames(gset) <- pdata\$Sample_Name
	
	beta <- getBeta(gset)
	#select appropriate columns depending on samprows and rows (cpg sites) depending on input file
	if(extension == 'rds'){{
		 B <- beta[ , c(samprows)]
	}}else {{ B <- beta[input\$CpG.Site, c(samprows)]}}
	B[B == 0] <- 0.0001
	B[B == 1] <- 1.0001
	x <- scale(t(B))
	x[is.nan(x)] <- 0
	pca <- prcomp(x)
	pdf("{output_file}.pdf")
        df <- as.data.frame(pca\$x)
        #df["Groups"] <- pData(gset)[,"Sample_Group"]
        df["Groups"] <- P\$group
	qplot(x=PC1, y=PC2, data=df, color=Groups)
        ggsave("{output_file}.pdf", plot=last_plot())
	EOF
        pdftoppm -rx 300 -ry 300 -png {output_file}.pdf {output_file}
        """.format(
			group=tuple(sample_group),
               	 	sample_names=tuple([sample.name for sample in contrast_samples]),
                	controls=', '.join(["'" + sample.name + "'" for sample in contrast.controls]),
                	cases=', '.join(["'" + sample.name + "'" for sample in contrast.treatments]),
                        input_rds_file=input_rds_file,
			input_file=input_file,
                        output_file=output_file),
		name = jobname)
                    

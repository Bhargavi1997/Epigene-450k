#!/usr/bin/env python 
import os
from core.config import *
from core.job import *


def build_model(rgset_file, dmps_file, contrast, method, source='minfi'):
	contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
        sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]

	return Job(
			[dmps_file,rgset_file],
			[contrast.name+'_'+source+'_model.RData'],                
                        command="""\                    
            mkdir -p report/data/models
	    R --vanilla <<-'EOF'
        suppressPackageStartupMessages(library(minfi))
            suppressPackageStartupMessages(library(caret))
            rgset <- readRDS("{rgset_file}")
            pdata <- pData(rgset) 
            beta <- getBeta(rgset)
	    colnames(beta) <- pdata$Sample_Name
            dmps <- read.csv("{dmps_file}")
            target_cpgs <- dmps$CpG.Site
            beta <- beta[match(target_cpgs, rownames(beta)),c{sample_names}]
            fit.model <- train(t(beta), y=as.factor(c{group}),method='{method}')
            save(fit.model, file='report/data/models/{file_name}')
EOF""".format(
			    group=tuple(sample_group),
			    sample_names=tuple([sample.name for sample in contrast_samples]),
                            rgset_file=rgset_file,
                            dmps_file=dmps_file,
			    method=method,
			    file_name = contrast.name+'_'+source + '_model.RData'
                            ),
                                name='build_model.'+contrast.name+'_'+source
                            )

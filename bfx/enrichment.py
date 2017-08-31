#!/usr/bin/env python 
import os
from core.config import *
from core.job import *

def pos_enr(grset_file, dmps_file, cores, LOLA_root, LOLA_user_dirs, LOLA_genome, LOLA_collection, analysis_file, contrast_name="user",source="user"):
	command="""\
    mkdir -p report/data/enrichment_analysis && \\
    R --vanilla <<-'EOF'
    #suppressPackageStartupMessages(library(BiSeq))
    library(BiSeq)
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(doParallel))
    # suppressPackageStartupMessages(library(matrixStats))
    suppressPackageStartupMessages(library(minfi))
    suppressPackageStartupMessages(library(LOLA))
    source('/hpf/largeprojects/ccmbio/jonBarenboim/mugqic_tools/R-tools/LOLAsearch.R')
    registerDoParallel(cores={cores})
    grset <- readRDS('{grset_file}')
    tmp <- read.csv('{dmps_file}')
    dmps <- tmp
    anno <- getAnnotation(grset)
    
    dmps$start <- anno[match(dmps$CpG.Site,anno$Name), 'pos']
    dmps$end <- dmps$start

    universe <- rowRanges(grset)
    userset <- granges(GRanges(dmps))
    # The universe should contain every region in the dmps. If not, something probably went wrong in earlier steps
    tryCatch(
        {{ checkUniverseAppropriateness(userset, universe) }},
        warning = function(w) {{ stop("userset is not a subset of universe") }},
        error = function(e) {{ stop("userset is not a subset of universe") }})
    # search LOLAcore for region sets
    LOLAcoreDB <- suppressWarnings(suppressMessages(buildRegionDB(rootdir='{LOLA_root}', genome='{LOLA_genome}',collection=c{LOLA_collection}, filename=c(), description=c(), any=c())))
    # Read region sets supplied by the user and coerce into the format expectd by LOLA
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
					cores=cores,
					grset_file=grset_file,
					dmps_file=dmps_file,
					LOLA_root=LOLA_root,
					LOLA_user_dirs=LOLA_user_dirs,
					LOLA_genome=LOLA_genome,
					LOLA_collection=LOLA_collection,
					analysis_file=analysis_file)

	job = Job(
	[dmps_file, grset_file],
	[analysis_file],
	[
		["position_enrichment_analysis", "module_R"],
		["position_enrichment_analysis", "module_pandoc"],
		["position_enrichment_analysis", "module_mugqic_tools"]
	],
	command=command,
	name="position_enrichment_analysis."+contrast_name+"_"+source)
	return job


#!/usr/bin/env python 
import os 
from core.config import *
from core.job import * 

#All normalization methods currently available in minfi
#SWAN, functional, raw, noob, quantile, illumina 
#Each method takes in a minfi RGChannel object and returns a Genomic ratio object

def swan(input_rds_file, output_rds_file):
	return Job(
		[input_rds_file],
		[output_rds_file],
		[
			["norm_swan", "module_R"],
			["norm_swan", "module_mugqic_R_packages"]
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
rgset <- readRDS("{input_rds_file}")
set.seed(33)
grset <-ratioConvert(mapToGenome(preprocessSWAN(rgset)))
saveRDS(grset, "{output_rds_file}")

EOF""".format(
		input_rds_file=input_rds_file,
		output_rds_file=output_rds_file
	))

def func(input_rds_file, output_rds_file):
	return Job(
		[input_rds_file],
		[output_rds_file],
		[
			["norm_func", "module_R"],
			["norm_func", "module_mugqic_R_packages"]
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
rgset <- readRDS("{input_rds_file}")
grset <- preprocessFunnorm(rgset)
saveRDS(grset, "{output_rds_file}")

EOF""".format(
		input_rds_file=input_rds_file,
		output_rds_file=output_rds_file
	))

def raw(input_rds_file, output_rds_file):
	return Job(
		[input_rds_file],
		[output_rds_file],
		[
			["norm_raw", "module_R"],
			["norm_raw", "module_mugqic_R_packages"]
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
rgset <- readRDS("{input_rds_file}")
grset <- ratioConvert(mapToGenome(preprocessRaw(rgset))
saveRDS(grset, "{output_rds_file}")

EOF""".format(
		input_rds_file=input_rds_file,
		output_rds_file=output_rds_file
	))

def noob(input_rds_file, output_rds_file):
	return Job(
		[input_rds_file],
		[output_rds_file],
		[
			["norm_noob", "module_R"],
			["norm_noob", "module_mugqic_R_packages"]
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
rgset <- readRDS("{input_rds_file}")
grset <- ratioConvert(mapToGenome(preprocessNoob(rgset)))
saveRDS(grset, "{output_rds_file}")

EOF""".format(
		input_rds_file=input_rds_file,
		output_rds_file=output_rds_file
	))

def quantile(input_rds_file, output_rds_file):
	return Job(
		[input_rds_file],
		[output_rds_file],
		[
			["norm_quantile", "module_R"],
			["norm_quantile", "module_mugqic_R_packages"]
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
rgset <- readRDS("{input_rds_file}")
grset <- preprocessQuantile(rgset)
saveRDS(grset, "{output_rds_file}")

EOF""".format(
		input_rds_file=input_rds_file,
		output_rds_file=output_rds_file
	))


def illumina(input_rds_file, output_rds_file):
	return Job(
		[input_rds_file],
		[output_rds_file],
		[
			["norm_illumina", "module_R"],
			["norm_illumina", "module_mugqic_R_packages"]
		],
		command="""\
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(minfi))
rgset <- readRDS("{input_rds_file}")
grset <- ratioConvert(mapToGenome(preprocessIllumina(rgset))
saveRDS(grset, "{output_rds_file}")

EOF""".format(
		input_rds_file=input_rds_file,
		output_rds_file=output_rds_file
	))

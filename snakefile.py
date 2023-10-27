import os
import pandas as pd
import tempfile
#import pprint
import ntpath
from Bio import SeqIO
import gzip
from collections import Counter
import itertools
import logging


configfile: "/mnt/data2/wuzengding/05.pipeline_dev/pipeline_SAKit_v1.0/config/config.yml"
configfile: config["SampleYml"]

#print(config)
SampleID = config["SampleID"]
outpath = config["OUTPATH"]
rule all:
    input:
        ## Isoform  outputs
        sqanti3_classification_final_results = os.path.join(outpath,"results/05.Isoform_Novel/{0}.isoform_classification_final_results.tsv".format(SampleID)),
        sqanti3_classification_filtered_out = os.path.join(outpath,"results/05.Isoform_Novel/{0}.isoform_classification_filtered_out_isoforms_results.tsv".format(SampleID)),

        ## Fusion outputs
        fusion_class_fusionhub = os.path.join(outpath,"results/06.Fusion_calling/{0}.fusion_classification_final_results_fusionhub.tsv".format(SampleID)),
        
        ## SVSNV outputs
        pbsv_vcf_annot_tsv = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.snpeff.tsv".format(SampleID)),
        snpeff_res = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.single_variant.snpeff.vcf".format(SampleID))

include: "workflow/rules/01.ccs.smk.py"
include: "workflow/rules/02.demultiplex.smk.py"
include: "workflow/rules/03.refine.smk.py"
include: "workflow/rules/04.cluster.smk.py"
include: "workflow/rules/05.align.smk.py"
include: "workflow/rules/06.1.Isoform_calling.smk.py"
include: "workflow/rules/06.2.Isoform_filter_bycounts.smk.py"
include: "workflow/rules/06.3.Isoform_filter_subset.smk.py"
include: "workflow/rules/12.Align_NGS.smk.py"
include: "workflow/rules/07.1.Isoform_quantfs.smk.py"
include: "workflow/rules/07.2.Isoform_classification.smk.py"
include: "workflow/rules/07.3.Isoform_classification_filter.smk.py"
include: "workflow/rules/08.Isoform_summary.smk.py"
include: "workflow/rules/09.1.Fusion_calling.smk.py"
include: "workflow/rules/09.2.Fusion_classification.smk.py"
include: "workflow/rules/10.Fusion_annotation.py"
include: "workflow/rules/11.SV_calling.smk.py"
include: "workflow/rules/13.SnvIndelCalling.smk.py"

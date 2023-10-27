rule Isoform_calling:
    input:
        sorted_sam = rules.Align.output.samsort_tx,
        hq_transcripts = rules.Cluster.output.hq_transcripts
    
    output:
        cupcake_collapsed_gff = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.gff".format(sample=SampleID)),
        cupcake_collapsed_rep_fa = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.rep.fa".format(sample=SampleID)),
        cupcake_collapsed_group = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.group.txt".format(sample=SampleID)),
        cupcake_ignored_ids = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.ignored_ids.txt".format(sample=SampleID))
    
    params:
        collapse = config["ScriptTools"]["collapse_isoforms"],
        short5merge = config["collapse"]["parameter"]["short5"],
        mincov = config["collapse"]["parameter"]["mincov"],
        miniden = config["collapse"]["parameter"]["miniden"],
        prefix = os.path.join(outpath,"results/04.Isoform_Calling/{sample}".format(sample=SampleID)),
        fq = ""


    log:
        os.path.join(outpath,"log/collapse.log")

    threads:
        1
    
    shell:
        """
            {config[Env][python3]} {params.collapse} \
            --input {input.hq_transcripts} \
            {params.fq} \
            -s  {input.sorted_sam}  \
            -o {params.prefix} {params.short5merge}  \
            -c {params.mincov} \
            -i {params.miniden} > {log} 2>&1
        """
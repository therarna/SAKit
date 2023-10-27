cutoff = config["isoform_filter_bycounts"]["parameter"]["cutoff"]
rule Isoform_filter_bycounts:
    input:
        cluster_report = rules.Cluster.output.cluster_report,
        collapsed_gff = rules.Isoform_calling.output.cupcake_collapsed_gff,
    
    output:
        cupcake_collapsed_read_stat = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.read_stat.txt".format(sample=SampleID)),
        cupcake_collapsed_abundance = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.abundance.txt".format(sample=SampleID)),
        cupcake_filter_gff = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.min_fl_{cutoff}.gff" \
                                  .format(sample=SampleID, cutoff = cutoff)),
        cupcake_filter_abundance = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.min_fl_{cutoff}.abundance.txt" \
                                  .format(sample=SampleID, cutoff = cutoff)),
        cupcake_filter_fa = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed.min_fl_{cutoff}.rep.fa" \
                                  .format(sample=SampleID, cutoff = cutoff))
                                  
    params:
        get_abun_script = config["ScriptTools"]["get_abundance"],
        get_abun_prefix = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed".format(sample=SampleID)),
        filter_by_count = config["ScriptTools"]["filter_by_count"],
        filt_bycount_prefix = os.path.join(outpath,"results/04.Isoform_Calling/{sample}.collapsed".format(sample=SampleID)),
        filt_bycount_cutoff = cutoff,
        
    log:
        os.path.join(outpath,"log/isoform_filter_bycounts.log")

    threads:
        1
    
    shell:
        """
           {config[Env][python3]} {params.get_abun_script} {params.get_abun_prefix} {input.cluster_report} > {log} 2>&1
           
           {config[Env][python3]} {params.filter_by_count}  --min_count {params.filt_bycount_cutoff} \
           --dun_use_group_count {params.filt_bycount_prefix} \
            > {log} 2>&1
        """
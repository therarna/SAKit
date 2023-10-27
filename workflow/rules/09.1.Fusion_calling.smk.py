
rule Fusion_calling:
    input:
        hq_transcripts = rules.Cluster.output.hq_transcripts,
        sorted_sam = rules.Align.output.samsort_tx,
        cluster_report = rules.Cluster.output.cluster_report

    output:
        fusion_gff = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion.gff".format(sample=SampleID)),
        fusion_rep = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion.rep.fa".format(sample=SampleID)),
        fusion_rename_rep = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion.rep.rename.fa".format(sample=SampleID)),
        fusion_abundance = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion.abundance.txt".format(sample=SampleID)),
        fusion_finder_stat= os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion.read_stat.txt".format(sample=SampleID)),

    params:
        min_locus_coverage = 0.05,
        min_total_coverage = 0.99,
        min_dist_between_loci = 10000,
        prefix_fusion_finder_out = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion".format(sample=SampleID))

    threads:
        1
    
    log:
        os.path.join(outpath,"log/Fusion_calling.log")

    shell:
        """
            {config[Env][python3]} {config[ScriptTools][fusion_finder]} \
                --min_locus_coverage {params.min_locus_coverage} \
                --min_total_coverage {params.min_total_coverage} \
                --min_dist_between_loci {params.min_dist_between_loci} \
                --input {input.hq_transcripts} \
                -s {input.sorted_sam} \
                -o {params.prefix_fusion_finder_out} \
                --cluster_report_csv {input.cluster_report} >{log} 2>&1
            
            sed  1,8d {output.fusion_abundance} > temp
            mv temp {output.fusion_abundance}
            
            sed 's/|.*//g' {output.fusion_rep} > {output.fusion_rename_rep}
        """
    


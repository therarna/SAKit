rule Isoform_summary:
    input:
        sqanti3_output_filtered = rules.Isoform_classification_filter.output[4],
        cupcake_collapsed_stat = rules.Isoform_filter_bycounts.output[0],
        sqanti3_output_unfiltered = rules.Isoform_classification.output[0]

    output:
        os.path.join(outpath,"results/05.Isoform_Novel/{0}.isoform_classification_final_results.tsv".format(SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{0}.isoform_classification_filtered_out_isoforms_results.tsv".format(SampleID)),
        os.path.join(outpath,"results/05.Isoform_Novel/{0}.isoform_disease_associate_gene.tsv".format(SampleID))
    params:
        pythonpath = config["Env"]["python_lib"],
        transcript2gene = config["Database"]["transcript2gene"],
        outdir = os.path.join(outpath,"results/05.Isoform_Novel"),
        disgenet = config["Database"]["disgenet"]
        

    threads:
        1

    log:
        os.path.join(outpath,"log/Isoform_summary.log")

    threads:
        1
    
    shell:
        """
           {config[Env][python3]} {config[InhouseScript][iosform_attibution]} \
           --collapsed_stat {input.cupcake_collapsed_stat} \
           --isoform_class {input.sqanti3_output_unfiltered} \
           --isoform_class_filter {input.sqanti3_output_filtered} \
           --tr2gene {params.transcript2gene} \
           --output1 {output[0]} \
           --output2 {output[1]} \
           > {log} 2>&1 
           
           {config[Env][python3]} {config[InhouseScript][gene_annotation_disease]} \
           --infile {output[0]} \
           --digenet {params.disgenet} \
           --output {output[2]} \
           >> {log} 2>&1 
        """

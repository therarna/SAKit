from workflow.scripts.iosform_attibution import attribute_ccs_read_counts

rule Fusion_annotation:
    input:
        fusion_class = rules.Fusion_classification.output.fusion_class,
        fusion_finder_stat= rules.Fusion_calling.output.fusion_finder_stat
    output:
        fusion_annot = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion.annotated.txt".format(sample=SampleID)),
        fusion_annot_fltrd = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion.annotated_ignored.txt".format(sample=SampleID)),
        fusion_classification_final_results = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion_classification_final_results.tsv".format(sample=SampleID)),
        fusion_class_fusionhub = os.path.join(outpath,"results/06.Fusion_calling/{sample}.fusion_classification_final_results_fusionhub.tsv".format(sample=SampleID))
        
    params:
        genome = config["Reference"]["genome"],
        genome_annotation = config["Reference"]["annotgtf"],
        transcript2gene = config["Database"]["transcript2gene"],
        min_fl_count = 2,
        fusion_finder_prefix = rules.Fusion_calling.params.prefix_fusion_finder_out,
        fusionhub_file = config["Database"]["fusionhubdb"]
        
    threads:
        1
    log:
        os.path.join(outpath,"log/Fusion_annotation.log")
    run:
        shell(\
            """
                {config[Env][python3]} {config[ScriptTools][fusion_collate_info]} \
                     --min_fl_count {params.min_fl_count} \
                     --genome {params.genome} \
                     {params.fusion_finder_prefix} \
                     {input.fusion_class} {params.genome_annotation} > {log} 2>&1
            """)
        
        fusion_results = attribute_ccs_read_counts( \
                                        output.fusion_annot, 
                                        input.fusion_finder_stat, 
                                        event_type = "fusion" \
                                        )
        try:
            fusion_results.to_csv( \
                            output.fusion_classification_final_results, 
                            sep ="\t", 
                            index =  False \
                            )
        except Exception as e:
            raise e
          
        shell(\
        """
           {config[Env][python3]} {config[InhouseScript][fusion_annotation]} \
           --fusion_class_final_result {output.fusion_classification_final_results} \
           --output_file {output.fusion_class_fusionhub} \
           --fusion_database {params.fusionhub_file} \
           --tr2gene {params.transcript2gene} \
           --thread_num {threads} >> {log} 2>&1
        """)
           
            
        
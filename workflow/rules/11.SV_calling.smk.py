
rule SV_calling:
    input:
        hq_mapped_hg38_bam = rules.Align.output.bamsort_tx

    output:
        pbsv_svsig = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.svsig.gz".format(SampleID)),
        pbsv_vcf = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.var.vcf".format(SampleID)),
        pbsv_vcf_annot = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.snpeff.vcf".format(SampleID)),
        report = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.snpeff.html".format(SampleID)),
        pbsv_vcf_annot_tsv = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant.snpeff.tsv".format(SampleID)),
    
    params:
        genome = config["Reference"]["genome"],
        pbsv = config["SoftwareTools"]["pbsv"],
        min_ref_span= config["PBsvParam"]["min_ref_span"],
        call_min_read_perc_one_sample = config["PBsvParam"]["call_min_read_perc_one_sample"],
        ccs_flag = "--ccs",
        snpeff = config["SoftwareTools"]["snpeff"],
        snpsift = config["SoftwareTools"]["snpsift"],
        species = config["Reference"]["species"],
        sep = "';'",
        empty = "'.'",
        prefix = os.path.join(outpath,"results/08.SvInDelSnvCalling/{0}.structure_variant".format(SampleID))



    threads:
        16

    log:
         os.path.join(outpath,"log/SV_calling.log")

    run:
        if params.species == 'hs':
            species = "hg38"
        elif  params.species == 'mm':
            species = "mm10"

        shell("""
            {params.pbsv} discover \
                -m {params.min_ref_span} \
                {input.hq_mapped_hg38_bam} {output.pbsv_svsig} > {log} 2>&1

            {params.pbsv} call  \
            -P {params.call_min_read_perc_one_sample} \
            {params.ccs_flag} \
            -j {threads} \
            {params.genome}  {output.pbsv_svsig}  {output.pbsv_vcf} >> {log} 2>&1

            /usr/bin/java -Xmx20g -jar {params.snpeff}/snpEff.jar \
            -c {params.snpeff}/snpEff.config \
            -noLog -v -stats {output.report} {species} \
            -fastaProt {params.prefix}.prot.fa  \
            {output.pbsv_vcf} > {output.pbsv_vcf_annot} 2> {log}
            
            /usr/bin/java -Xmx20g -jar {params.snpsift}/SnpSift.jar extractFields \
                 -s {params.sep} -e {params.empty} {output.pbsv_vcf_annot}  CHROM POS END SVLEN SVTYPE FILTER ANN[*].GENE ANN[*].IMPACT ANN[*].EFFECT ANN[*].RANK > {output.pbsv_vcf_annot_tsv} 2> {log}
            """)
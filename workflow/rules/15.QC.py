rule QC:
    input:
        bamfile = rules.Align.output.flnc_bamsort_tx
    output:
        qcsum = os.path.join(outpath,"results/00.QC/{0}.runqc.txt".format(SampleID))
    log:
        log_geneBody_coverage = os.path.join(outpath,"log/{0}.QC.geneBody_coverage.log".format(SampleID)),
        log_FPKM_count_genome = os.path.join(outpath,"log/{0}.QC.log_FPKM_count_genome.log".format(SampleID)),
        log_FPKM_count_wgs = os.path.join(outpath,"log/{0}.QC.log_FPKM_count_wgs.log".format(SampleID)),
        log_junction_saturation = os.path.join(outpath,"log/{0}.QC.log_junction_saturation.log".format(SampleID)),
        log_tin = os.path.join(outpath,"log/{0}.QC.log_tin.log".format(SampleID)),
        log_CollectRnaSeqMetrics = os.path.join(outpath,"log/{0}.QC.log_CollectRnaSeqMetrics.log".format(SampleID))
        
    params:
        pytool = config["Env"]["python3"],
        qctool = config["SoftwareTools"]["qctool"],
        picard = config["SoftwareTools"]["picard"],
        rseqc = config["SoftwareTools"]["rseqc"],
        bedfile = config["Reference"]["bedfile"],
        rRNA_bed = config["Reference"]["rRNA_bed"],
        rRNA_interval = config["Reference"]["rRNA_interval"],
        refFlat = config["Reference"]["refFlat"],
        outpath = os.path.join(outpath,"results/00.QC"),
        java = config["Env"]["java"]
        
    threads: 20
    
    shell:
        """ 
        
        {params.pytool} {params.rseqc}/geneBody_coverage.py -i {input.bamfile} \
        -r {params.bedfile} -o {params.outpath}/{SampleID} > {log.log_geneBody_coverage} 2>&1
        
        {params.pytool} {params.rseqc}/bam_stat.py -i {input.bamfile} > \
        {params.outpath}/{SampleID}.bamstat.txt 2>&1 
        
        {params.pytool} {params.rseqc}/infer_experiment.py -i {input.bamfile} \
        -r {params.bedfile} >{params.outpath}/{SampleID}.experiment.txt 2>&1 
        
        {params.pytool} {params.rseqc}/FPKM_count.py -i {input.bamfile} \
        -r {params.bedfile} -o {params.outpath}/{SampleID}.genome > {log.log_FPKM_count_genome} 2>&1
    
        {params.pytool} {params.rseqc}/FPKM_count.py -i {input.bamfile} \
        -r {params.bedfile} -o {params.outpath}/{SampleID}.wgs > {log.log_FPKM_count_wgs} 2>&1
         
        {params.pytool} {params.rseqc}/junction_saturation.py -i {input.bamfile} \
        -r {params.bedfile} -o {params.outpath}/{SampleID}.rRNA > {log.log_junction_saturation} 2>&1
        
        {params.pytool} {params.rseqc}/read_distribution.py -i {input.bamfile} \
        -r {params.bedfile} > {params.outpath}/{SampleID}.read_distribution.txt 2>&1
        
        cd {params.outpath}
        {params.pytool} {params.rseqc}/tin.py -i {input.bamfile} \
        -r {params.bedfile} > {log.log_tin} 2>&1
        
        {params.java} -jar {params.picard} CollectRnaSeqMetrics --INPUT {input.bamfile} \
        --OUTPUT {params.outpath}/{SampleID}.RnaSeqMetrix.tsv \
        --VALIDATION_STRINGENCY  LENIENT --REF_FLAT {params.refFlat}  \
        --RIBOSOMAL_INTERVALS {params.rRNA_interval} --STRAND_SPECIFICITY NONE > {log.log_CollectRnaSeqMetrics} 2>&1
        """

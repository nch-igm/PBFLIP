
rule run_sv_caller:
    input:
        hq_mapped_hg38_bam = rules.run_sam_to_bam.output[0]

    output:
        pbsv_svsig = rules.all.input[38],
        pbsv_vcf =  rules.all.input[39],
        pbsv_vcf_annot = rules.all.input[40],
        report = rules.all.input[41],
        pbsv_vcf_annot_tsv = rules.all.input[42]

    
    params:
        genome = config["REFERENCES"]["genome"],
        min_ref_span= config["PBSVCALLERPARAM"]["min_ref_span"],
        call_min_read_perc_one_sample = config["PBSVCALLERPARAM"]["call_min_read_perc_one_sample"],
        ccs_flag = "--ccs",
        snpeff_java_mem = "-Xmx8g",
        snpeff = config["SNPEFF"]["default"],
        snpsift = config["SNPEFF"]["snpsift"],
        species = config["REFERENCES"]["species"],
        sep = "';'",
        empty = "'.'"



    threads:
        16

    log:
        "log/run_sv_caller.log"

    run:
        if params.species == 'hs':
            species = "hg38"
        elif  params.species == 'mm':
            species = "mm10"

        shell("""
            pbsv discover \
                -m {params.min_ref_span} \
                {input.hq_mapped_hg38_bam} {output.pbsv_svsig} > {log} 2>&1

            pbsv call  \
            -P {params.call_min_read_perc_one_sample} \
            {params.ccs_flag} \
            -j {threads} \
            {params.genome}  {output.pbsv_svsig}  {output.pbsv_vcf} >> {log} 2>&1

            java {params.snpeff_java_mem} \
                 -jar {params.snpeff} \
                 -v -stats {output.report} {species} {output.pbsv_vcf} > {output.pbsv_vcf_annot} 2> {log}
            
             java {params.snpeff_java_mem} \
                 -jar {params.snpsift} extractFields \
                 -s {params.sep} -e {params.empty} {output.pbsv_vcf_annot}  CHROM POS END SVLEN SVTYPE FILTER ANN[*].GENE ANN[*].IMPACT ANN[*].EFFECT ANN[*].RANK > {output.pbsv_vcf_annot_tsv} 2> {log}
            

        """)

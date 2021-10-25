
rule run_sam_to_bam:
    input:
        rules.run_samsort.output
    
    output:
        rules.all.input[32],
        hq_isoforms_mapped_bai = "results/FinalResults/{0}_hq_isoforms_mapped_hg38_sorted.bam.bai".format(CASENAME)
    
    log:
        "log/run_sam_to_bam.log"

    threads:
        18
    
    shell:
        """
            samtools view  \
                -@ {threads} \
                -S -b -h \
                -o {output[0]} {input} > {log}
            samtools index {output} >> {log}
        """
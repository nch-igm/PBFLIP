import ntpath

rule run_short_reads_genome_alignment:
    output:
        rules.all.input[16]
    
    log:
        "log/run_short_reads_genome_alignment.log"

    threads:
        1
    
    shell:
        """           
            touch {output}
        """

rule run_transcriptome_indexing:
    input:
        rules.all.input[14]
    output:
        rules.all.input[17]
    
    log:
        "log/run_transcriptome_indexing.log"

    threads:
        1
    
    shell:
        """
            touch {output}
        """


rule run_transcriptome_alignment:
    output:
        rules.all.input[18]
    
    log:
        "log/run_transcriptome_alignment.log"

    threads:
        1
    
    shell:
        """
           touch {output}
        """

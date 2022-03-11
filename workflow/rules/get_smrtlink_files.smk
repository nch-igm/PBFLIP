
rule run_get_smrtlink_clusterreport:
    input:
        cluster_report = config["SMRTLINKFILES"]["cluster_report"]

    output:
        cluster_report = "results/smrt_link_data/cluster_report.csv"


    log:
        "log/get_smrtlink_files_cluster_report.log"

    threads:
        1

    run:
        if(os.path.isfile(input[0])):
            shell( 
                """
                cp {input.cluster_report} {output.cluster_report} > {log}
                """
            )
        else:
            print ( "No such a file called unpolished.cluster_report.csv")
            exit


rule run_get_smrtlink_hq_tx:
    input:
        hq_transcripts = config["SMRTLINKFILES"]["hq_transcripts"]

    output:
        hq_transcripts = "results/smrt_link_data/hq_isoforms.fasta",

    log:
        "log/get_smrtlink_files_hq_tx.log"

    threads:
        1

    run:
        if(os.path.isfile(input[0])):
            shell( 
                    """
                    cat {input.hq_transcripts} | sed 's/.*_/>/' > {output.hq_transcripts}
		    echo 'hq_transcripts header changed' >> {log}
                    """
                )
            
        else:
            print ( "No such a file called hq_transcripts.fasta")
            exit


rule run_get_smrtlink_flnc:
    input:
        flnc = config["SMRTLINKFILES"]["flnc"]

    output:
        flnc_bam = "results/smrt_link_data/flnc.bam",
        ccs_fq = "results/smrt_link_data/ccs.fastq"


    log:
        "log/get_smrtlink_files_flnc.log"

    params:
        prefix_output =  "results/smrt_link_data/ccs",
        dont_compress = "-u"

    threads:
        1

    run:
        if(os.path.isfile(input[0])):
            shell( 
                    """
                    cp {input.flnc} {output.flnc_bam} > {log}
                    {config[PBBAM][index]} \
                        {output.flnc_bam} >> {log}
                    {config[PBBAM][bam2fastq]} \
                        {params.dont_compress} \
                        -o {params.prefix_output} \
                        {output.flnc_bam} >> {log}
                    """
                )

        else:
            print ( "No such a file called ccs.fastq.zip")
            exit

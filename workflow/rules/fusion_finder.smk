
rule run_fusion_finder:
    input:
        hq_transcripts = rules.run_get_smrtlink_hq_tx.output.hq_transcripts,
        sorted_sam = rules.run_samsort.output[0],
        cluster_report = rules.run_get_smrtlink_clusterreport.output.cluster_report

    output:
        fusion_gff = rules.all.input[47],
        fusion_rep =  rules.all.input[48],
        fusion_abundance = rules.all.input[49],
        fusion_finder_stat= "results/FusionProducts/lq_isoforms.fasta.fusion.read_stat.txt",

    params:
        min_locus_coverage = 0.05,
        min_total_coverage = 0.99,
        min_dist_between_loci = 10000,
        prefix_fusion_finder_out = "results/FusionProducts/lq_isoforms.fasta.fusion"

    threads:
        1
    
    log:
        "log/run_fusion_finder.log"

    shell:
        """
            {config[ISOSEQSCRIPTS][fusion_finder]} \
                --min_locus_coverage {params.min_locus_coverage} \
                --min_total_coverage {params.min_total_coverage} \
                --min_dist_between_loci {params.min_dist_between_loci} \
                --input {input.hq_transcripts} \
                -s {input.sorted_sam} \
                -o {params.prefix_fusion_finder_out} \
                --cluster_report_csv {input.cluster_report} >{log} 2>&1
            
            sed  1,8d {output.fusion_abundance} > temp
            mv temp {output.fusion_abundance}

        """
    


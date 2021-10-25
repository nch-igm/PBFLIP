rule run_fusion_collect_info:
    input:
        fusion_product_class = rules.run_fusion_classifier.output.fusion_product_class
    output:
        fusion_annot = rules.all.input[52],
        fusion_annot_fltrd = "results/FusionProducts/lq_isoforms.fasta.fusion.annotated_ignored.txt"
        
    params:
        genome = config["REFERENCES"]["genome"],
        genome_annotation = config["REFERENCES"]["annotation"],
        min_fl_count = 2,
        fusion_finder_prefix = rules.run_fusion_finder.params.prefix_fusion_finder_out

    threads:
        1
    log:
        "log/run_fusion_collect_info.log"
    shell:
        """
            {config[ISOSEQSCRIPTS][fusion_collate_info]} \
                 --min_fl_count {params.min_fl_count} \
                 {params.fusion_finder_prefix}  \
                 {input.fusion_product_class} {params.genome_annotation} \
                 --genome {params.genome} > {log} 2>&1
        """

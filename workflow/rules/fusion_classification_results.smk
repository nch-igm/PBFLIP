include: "sqanti3_classification_results.smk"


rule run_fusion_classification_results:
    input:
        fusion_collate_info_annot = rules.run_fusion_collect_info.output.fusion_annot,
        fusion_finder_stat= "results/FusionProducts/lq_isoforms.fasta.fusion.read_stat.txt"

    output:
        fusion_classification_final_results = rules.all.input[53]

    threads:
        1
    log:
        "log/run_fusion_classification_results.log"
    run:
        fusion_product_results = associate_ccs_read_counts( \
                                        input.fusion_collate_info_annot, 
                                        input.fusion_finder_stat, 
                                        event_type = "fusion" \
                                        )

        try:
            fusion_product_results.to_csv( \
                            output.fusion_classification_final_results, 
                            sep ="\t", 
                            index =  False \
                            )
        except Exception as e:
            raise e


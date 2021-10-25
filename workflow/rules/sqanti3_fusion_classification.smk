rule run_fusion_classifier:
    input:
        fusion_product_gff = rules.run_fusion_finder.output.fusion_gff

    output:
        fusion_product_class = rules.all.input[50],
        fusion_product_junc = rules.all.input[51],
        refAnnotation_lq_isoforms = "results/FusionClassification/refAnnotation_lq_isoforms.fasta.fusion.genePred",
        fusion_gtf_corrected = "results/FusionClassification/lq_isoforms.fasta.fusion_corrected.gtf",
        fusion_rep_corrected = "results/FusionClassification/lq_isoforms.fasta.fusion_corrected.fasta"

    params:
        genome = config["REFERENCES"]["genome"],
        genome_annotation = config["REFERENCES"]["annotation"],
        classification_dir = "results/FusionClassification/",
        classification_prefix = "lq_isoforms.fasta.fusion"

    threads:
        1
    log:
        "log/run_fusion_classifier.log"

    shell:
        """
            {config[ISOSEQSCRIPTS][sqanti_qc]} \
                --is_fusion \
                -d {params.classification_dir} \
                -o {params.classification_prefix} \
                {input.fusion_product_gff} \
                {params.genome_annotation} \
                {params.genome}   > {log} 2>&1
        """


rule run_collapse:
    input:
        sorted_sam = rules.run_samsort.output,
        hq_tx = rules.run_get_smrtlink_hq_tx.output.hq_transcripts
    
    output:
        rules.all.input[5],
        rules.all.input[6],
        rules.all.input[7],
        cupcake_ignored_ids = "results/Collapsed_Isoforms/cupcake.ignored_ids.txt"
    
    params:
        pythonpath = config["LIBPATHS"]["python_lib"],
        collapse = config["ISOSEQSCRIPTS"]["collapse_isoforms"],
        collapse_param = config["COLLAPSEPARAM"]["default"],
        prefix = "results/Collapsed_Isoforms/cupcake",
        fq = ""


    log:
        "log/run_collapse.log"

    threads:
        1
    
    shell:
        """
            PYTHONPATH={params.pythonpath}
            {params.collapse} \
            --input {input.hq_tx} \
            {params.fq} \
            -s  {input.sorted_sam}  \
            -o {params.prefix} {params.collapse_param} > {log} 2>&1
        """
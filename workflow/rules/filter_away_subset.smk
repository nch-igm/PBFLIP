
rule run_filter_away_subset:
    input:
        rules.run_filter_by_counts.output
    
    output:
        rules.all.input[12],
        rules.all.input[13],
        rules.all.input[14],
        rules.all.input[15]
    
    params:
        pythonpath = config["LIBPATHS"]["python_lib"],
        filter_away_subset = config["ISOSEQSCRIPTS"]["filter_away_subset"],
        prefix = "results/Collapsed_Isoforms/cupcake.collapsed.min_fl_{filter_count_cutoff}".format(filter_count_cutoff = cutoff ),
        sed = "'s/|.*//g'"

    log:
        "log/run_filter_away_subset.log"

    threads:
        1
    
    shell:
        """
            PYTHONPATH={params.pythonpath}
            {params.filter_away_subset} {params.prefix}> {log} 2>&1
            sed {params.sed} {output[1]}> {output[2]}
        """

rule run_filter_by_counts:
    input:
        rules.run_get_abundance.output
    
    output:
        rules.all.input[9],
        rules.all.input[10],
        rules.all.input[11]
    
    params:
        pythonpath = config["LIBPATHS"]["python_lib"],
        filter_by_count = config["ISOSEQSCRIPTS"]["filter_by_count"],
        filter_by_count_parameters = config["FILTERBYCOUNTS"]["default"],
        prefix = "results/Collapsed_Isoforms/cupcake.collapsed"


    log:
        "log/run_filter_by_counts.log"

    threads:
        1
    
    shell:
        """
            PYTHONPATH={params.pythonpath}
           {params.filter_by_count}  {params.filter_by_count_parameters} {params.prefix} > {log} 2>&1
        """
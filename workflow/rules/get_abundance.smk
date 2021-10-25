
rule run_get_abundance:
    input:
        cluster_report = rules.run_get_smrtlink_clusterreport.output,
        collapsed_gff = rules.run_collapse.output[0]
    
    output:
        rules.all.input[31],
        rules.all.input[8]
    
    params:
        get_abundance = config["ISOSEQSCRIPTS"]["get_abundance"],
        prefix = "results/Collapsed_Isoforms/cupcake.collapsed"


    log:
        "log/run_get_abundance.log"

    threads:
        1
    
    shell:
        """
           {params.get_abundance} {params.prefix} {input.cluster_report} > {log} 2>&1
        """
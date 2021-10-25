rule run_hq_tx_mapping:
    input:
        hq_tx = rules.run_get_smrtlink_hq_tx.output.hq_transcripts

    output:
        sam_tx = "results/Minimap2_Mapping/hq_isoforms_mapped_hg38.sam"

    params:
        genome = config["REFERENCES"]["genome"],
        ax = "splice", 
        secondary = "no",  
        o = "-O6,24", 
        b4 = "-B4",
        uf = "-uf",
        read_groups = "'@RG\\tID:HQ\\tSM:{sample}'".format(sample = CASENAME),
        hard_clip_off = "-Y"

    log:
        "log/run_hq_tx_mapping.log"

    threads:
        18

    shell:
        """
            {config[MAPPERS][mapper_default]} \
                -t {threads} \
                -ax {params.ax} \
                --secondary={params.secondary} \
                -R {params.read_groups} \
                {params.hard_clip_off} \
                {params.o} \
                {params.b4} \
                {params.uf} \
                {params.genome} \
                {input} > {output} 2> {log}

        """

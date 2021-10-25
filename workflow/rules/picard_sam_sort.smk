
rule run_samsort:
    input:
        rules.run_hq_tx_mapping.output.sam_tx
    
    output:
        "results/Minimap2_Mapping/hq_isoforms_mapped_hg38.sorted.sam"
    
    params:
        picard = config["PICARD"]["default"],
        mem = "-Xmx40g",
        parallel_threads = "-XX:ParallelGCThreads=1",
        temp_dir =  "-Djava.io.tmpdir=/igm/temp",
        sort_order = "coordinate"

    log:
        "log/run_samsort.log"

    threads:
        1
    
    shell:
        """
            java {params.mem} {params.parallel_threads}  {params.temp_dir} \
                -jar  {params.picard} \
                SortSam \
                I={input} \
                O={output} \
                SORT_ORDER={params.sort_order}  2> {log}
        """
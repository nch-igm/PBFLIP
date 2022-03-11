rule run_classified_tx_mapping:
    input:
        squanti3_output_filtered  = rules.run_sqanti3_filter.output[5]

    output:
        squanti3_output_filtered_bam_tx = rules.all.input[33],
        collapsed_filtered_hq_bam_lite_bai = rules.all.input[105]

    params:
        genome = config["REFERENCES"]["genome"],
        ax = "splice", 
        secondary = "no",  
        o = "-O6,24", 
        b4 = "-B4",
        uf = "-uf",
        temp_out_sam = "results/FinalResults/{0}_sqanti3_classification.filtered_lite_mapped_hg38_sorted.sam".format(CASENAME),
        picard = config["PICARD"]["default"],
        mem = "-Xmx40g",
        parallel_threads = "-XX:ParallelGCThreads=1",
        temp_dir =  "-Djava.io.tmpdir=/data/tmp",
        sort_order = "coordinate"

    log:
        "log/run_classified_tx_mapping.log"

    threads:
        18

    shell:
        """
            {config[MAPPERS][mapper_default]} \
                -t {threads} \
                -ax {params.ax} \
                --secondary={params.secondary} \
                {params.o} \
                {params.b4} \
                {params.uf} \
                {params.genome} \
                {input.squanti3_output_filtered} > {params.temp_out_sam} 2> {log}
            
            java {params.mem} {params.parallel_threads}  {params.temp_dir} \
                -jar  {params.picard} \
                SortSam \
                I={params.temp_out_sam} \
                O={output.squanti3_output_filtered_bam_tx} \
                SORT_ORDER={params.sort_order}  2>> {log}
            
            samtools index {output.squanti3_output_filtered_bam_tx} 2>> {log}

            rm {params.temp_out_sam}
        """

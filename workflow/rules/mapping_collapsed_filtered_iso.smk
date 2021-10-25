rule run_collapsed_filtered_hq_mapping:
    input:
        hq_tx = rules.run_filter_away_subset.output[2]

    output:
        collapsed_filtered_hq_bam_tx = rules.all.input[68]

    params:
        genome = config["REFERENCES"]["genome"],
        ax = "splice", 
        secondary = "no",  
        o = "-O6,24", 
        b4 = "-B4",
        uf = "-uf",
        temp_out_sam = "results/FinalResults/{0}_cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered_mapped_hg38_sorted.sam".format(CASENAME, filter_count_cutoff = cutoff),
        picard = config["PICARD"]["default"],
        mem = "-Xmx40g",
        parallel_threads = "-XX:ParallelGCThreads=1",
        temp_dir =  "-Djava.io.tmpdir=/data/tmp",
        sort_order = "coordinate"

    log:
        "log/run_collapsed_filtered_hq_mapping.log"

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
                {input.hq_tx} > {params.temp_out_sam} 2> {log}
            
            java {params.mem} {params.parallel_threads}  {params.temp_dir} \
                -jar  {params.picard} \
                SortSam \
                I={params.temp_out_sam} \
                O={output.collapsed_filtered_hq_bam_tx} \
                SORT_ORDER={params.sort_order}  2>> {log}
            
            samtools index {output.collapsed_filtered_hq_bam_tx} 2>> {log}

            rm {params.temp_out_sam}
        """

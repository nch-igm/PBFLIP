rule run_sqanti3_classification:
    input:
        collapsed_isoform_gtf = rules.run_filter_away_subset.output[0],
        #orf_input = rules.run_filter_away_subset.output[2],
        fl_count_file = rules.run_filter_away_subset.output[3]

    
    output:
        rules.all.input[19],
        rules.all.input[20],
        rules.all.input[21],
        rules.all.input[22]
    
    params:
        pythonpath = config["LIBPATHS"]["python_lib"],
        aligner_choice = "minimap2",
        chunks = 5,
        isoannot = "--isoAnnotLite",
        outprefix = "sqanti3",
        output_dir = os.path.join(workflow.basedir, "results/IsoformClassification/"),
        isoannotlitegff3 = config["REFERENCES"]["isoannotlitegff3"],
        genome = config["REFERENCES"]["genome"],
        genome_annotation = config["REFERENCES"]["annotation"],
        orf_input_path = os.path.join(workflow.basedir, rules.run_filter_away_subset.output[2])



    log:
        "log/run_sqanti3_classification.log"

    threads:
        18
    
    run:
        if not config["ILLUMINASHORTREADS"]["ill_fastq_R1"]:
            shell(\
                """
                    PYTHONPATH={params.pythonpath}
                    {config[ISOSEQSCRIPTS][sqanti_qc]} \
                        --dir {params.output_dir} \
                        -t {threads} \
                        --aligner_choice {params.aligner_choice} \
                        --chunks {params.chunks}  \
                        -o {params.outprefix} \
                        {params.isoannot} \
                        --gff3 {params.isoannotlitegff3}  \
                        --fl_count {input.fl_count_file}  \
                        --orf_input {params.orf_input_path} \
                        {input.collapsed_isoform_gtf} {params.genome_annotation}  {params.genome} >{log} 2>&1
                """
                )
        else:
            junction_shortreads = rules.run_short_reads_genome_alignment.params.output_prefix + "SJ.out.tab"
            expression = rules.run_transcriptome_alignment.params.quant_file
            shell(\
                """
                    PYTHONPATH={params.pythonpath}
                     {config[ISOSEQSCRIPTS][sqanti_qc]} \
                        --dir {params.output_dir} \
                        -t {threads} \
                        --aligner_choice {params.aligner_choice} \
                        --chunks {params.chunks}  \
                        -o {params.outprefix} \
                        {params.isoannot} \
                        --gff3 {params.isoannotlitegff3}  \
                        --fl_count {input.fl_count_file}  \
                        --orf_input {params.orf_input_path} \
                        --expression {expression} \
                        --coverage {junction_shortreads} \
                        {input.collapsed_isoform_gtf} {params.genome_annotation}  {params.genome}  >{log} 2>&1
                """
                 )
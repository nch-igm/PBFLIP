rule run_sqanti3_filter:
    input:
        sqanti3_classification = rules.run_sqanti3_classification.output[0],
        sqanti3_classification_corrected_fa = rules.run_sqanti3_classification.output[2],
        sqanti3_classification_corrected_gtf = rules.run_sqanti3_classification.output[3]

    
    output:
        rules.all.input[23],
        rules.all.input[24],
        rules.all.input[25],
        rules.all.input[26],
        rules.all.input[27],
        rules.all.input[28]
    
    params:
        pythonpath = config["LIBPATHS"]["python_lib"],
        intrapriming = 0.6,
        runAlength = 6,
        max_dist_to_known_end = 50,
        min_cov = 3,
        input_dir = os.path.join(workflow.basedir, "results/IsoformClassification/")


    log:
        "log/run_sqanti3_filter.log"

    threads:
        1
    
    run:
        shell(\
            """
                PYTHONPATH={params.pythonpath}
                {config[ISOSEQSCRIPTS][sqant_filter]}  \
                --intrapriming {params.intrapriming} \
                --runAlength {params.runAlength} \
                --max_dist_to_known_end {params.max_dist_to_known_end} \
                --min_cov {params.min_cov} \
                {input.sqanti3_classification} \
                {input.sqanti3_classification_corrected_fa} \
                {input.sqanti3_classification_corrected_gtf} 2>{log}
        """)

        outputs = [output[0], output[1], output[2], output[3], output[4],  output[5]]
        for f in outputs:
            file_name = ntpath.basename( f )
            input_file_path = os.path.join(params.input_dir, file_name )
            shell("mv {input_file_path} {f} >> {log}" ) 
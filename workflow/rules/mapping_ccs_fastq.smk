
def filter_css_reads( cupcake_read_stat):

    read_cluster_report_pd =  pd.read_csv( cupcake_read_stat, sep ='\t')
    read_cluster_report_pd_filtered_FL = read_cluster_report_pd[ read_cluster_report_pd[ 'is_fl']== 'Y']

    fl_read_name_list = set(read_cluster_report_pd_filtered_FL['id'])


    fd, path = tempfile.mkstemp()
    try:
        with os.fdopen(fd, 'w') as tmp:
            # do stuff with temp file
            tmp.writelines(line + '\n' for line in fl_read_name_list)
    
    except:
        raise
    return path
            



rule run_mapping_css_fastq:
    input:
        cupcake_collapsed_stat = rules.run_get_abundance.output[0],
        ccs_fastq_file = rules.run_get_smrtlink_flnc.output.ccs_fq

    
    output:
        ccs_fl_fastq = rules.all.input[34],
        ccs_fl_bam = rules.all.input[35]

    params:
        pythonpath = config["LIBPATHS"]["python_lib"],
        genome = config["REFERENCES"]["genome"],
        ax = "splice/splice:hq", 
        temp_out_sam = "results/FinalResults/{0}_SQII_ccs_fl.sam".format(CASENAME),
        picard = config["PICARD"]["default"],
        mem = "-Xmx40g",
        parallel_threads = "-XX:ParallelGCThreads=1",
        temp_dir =  "-Djava.io.tmpdir=/igm/temp/",
        sort_order = "coordinate"

    log:
        "log/run_mapping_css_fastq.log"

    threads:
        18
    
    run:
        file_name = filter_css_reads( input[0] )
        shell(\
                """
                    filterbyname.sh  in={input.ccs_fastq_file}  \
                                       names={file_name} \
                                       out={output.ccs_fl_fastq} \
                                       overwrite=true \
                                       include=true \
                                       qin=33 >{log} 2>&1

                """    
        )
        shell(\
                """
                    {config[MAPPERS][mapper_default]} \
                        -t {threads} \
                        -ax {params.ax} \
                        {params.genome} \
                        {output.ccs_fl_fastq} > {params.temp_out_sam} 2>> {log}
                """
        )

        shell(\
                """
                    java {params.mem} {params.parallel_threads}  {params.temp_dir} \
                        -jar  {params.picard} \
                        SortSam \
                        I={params.temp_out_sam} \
                        O={output.ccs_fl_bam} \
                        SORT_ORDER={params.sort_order}  2>> {log}
                    
                    samtools index {output.ccs_fl_bam} 2>> {log}
                    rm {params.temp_out_sam}
                """
        )

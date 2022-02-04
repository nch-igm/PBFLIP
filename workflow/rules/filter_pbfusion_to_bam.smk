def filter_fusion_reads( final_pbfusion_fa, fusion_class_cor_fasta, fusion_final_results):
    if os.stat(fusion_final_results).st_size != 0:
        fusion_final_results_pd =  pd.read_csv( fusion_final_results, sep ='\t')

        pbfusion_name_list  = fusion_final_results_pd ["UniqueID"]

        with open( final_pbfusion_fa, "w" ) as output_handle:
            for seq_record in SeqIO.parse( fusion_class_cor_fasta, "fasta"):
                for pb in pbfusion_name_list:
                    if pb.lower() in seq_record.description.lower():
                        #print(seq_record.format("fasta"))
                        sequences =seq_record
                        SeqIO.write(sequences, output_handle, "fasta")
                        break
        return 0
    else:
        return 1


rule run_filter_pbfusion_to_bam:
    input:
        fusion_rep = rules.run_fusion_finder.output.fusion_rep,
        fusion_final_results = rules.run_fusion_events.output.final_fusion_class_fusionhub

    output:
        final_pbfusion_fa = rules.all.input[55],
        final_pbfusion_mapped_bam = rules.all.input[56],
        final_pbfusion_mapped_bam_bai = "results/FinalResults/{0}_final_pbfusion_mapped_hg38_sorted.bam.bai".format(CASENAME)

    params:
        genome = config["REFERENCES"]["genome"],
        ax = "splice", 
        secondary = "no",  
        o = "-O6,24", 
        b4 = "-B4",
        uf = "-uf",
        temp_out_sam = "results/FinalResults/final_pbfusion.sam",
        picard = config["PICARD"]["default"],
        mem = "-Xmx40g",
        parallel_threads = "-XX:ParallelGCThreads=1",
        temp_dir =  "-Djava.io.tmpdir=/data/tmp",
        sort_order = "coordinate"
    threads:
        18
    log:
        "log/run_filter_pbfusion_to_bam.log"
    run:
        filter_status = filter_fusion_reads( \
                                output.final_pbfusion_fa, 
                                input.fusion_rep, 
                                input.fusion_final_results \
                            )
        if filter_status == 0:
            shell(\
                    """
                        {config[MAPPERS][mapper_default]} \
                            -t {threads} \
                            -ax {params.ax} \
                            --secondary={params.secondary} \
                            {params.o} \
                            {params.b4} \
                            {params.uf} \
                            {params.genome} \
                            {output.final_pbfusion_fa} > {params.temp_out_sam} 2> {log}
                
                        java {params.mem} {params.parallel_threads}  {params.temp_dir} \
                            -jar  {params.picard} \
                            SortSam \
                            I={params.temp_out_sam} \
                            O={output.final_pbfusion_mapped_bam} \
                            SORT_ORDER={params.sort_order}  2>> {log}
                
                        samtools index {output.final_pbfusion_mapped_bam} 2>> {log}

                        rm {params.temp_out_sam}
                    """  
                )
        else:
            shell(\
                    """
                        touch {output.final_pbfusion_fa}

                        touch {output.final_pbfusion_mapped_bam}
                        
                        touch {output.final_pbfusion_mapped_bam_bai}
                    """
                )

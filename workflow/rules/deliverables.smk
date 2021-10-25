rule run_deliverables:
    input:
        collapsed_gff = rules.run_collapse.output[0],
        collapsed_rep_fa = rules.run_collapse.output[1],
        collapsed_group = rules.run_collapse.output[2],
        collapsed_abundance = rules.run_get_abundance.output[1],
        cupcake_ignored_ids = "results/Collapsed_Isoforms/cupcake.ignored_ids.txt",
        collapsed_gff_fl = rules.run_filter_by_counts.output[0],
        collapsed_rep_fa_fl =  rules.run_filter_by_counts.output[2],
        collapsed_abundance_fl =  rules.run_filter_by_counts.output[1],
        collapsed_gff_fl_away = rules.run_filter_away_subset.output[0],
        collapsed_rep_fa_fl_away =  rules.run_filter_away_subset.output[1],
        collapsed_abundance_fl_away =  rules.run_filter_away_subset.output[3],
        #collapsed_filtered_hq_bam_tx  = rules.run_collapsed_filtered_hq_mapping.output.collapsed_filtered_hq_bam_tx
        sqanti3_class = rules.run_sqanti3_classification.output[0],
        sqanti3_junc = rules.run_sqanti3_classification.output[1],
        sqanti3_corrected_fa = rules.run_sqanti3_classification.output[2],
        sqanti3_corrected_gtf = rules.run_sqanti3_classification.output[3],
        sqanti3_class_flt = rules.run_sqanti3_filter.output[4],
        sqanti3_junc_flt = rules.run_sqanti3_filter.output[3],
        sqanti3_corrected_fa_flt = rules.run_sqanti3_filter.output[5],
        sqanti3_corrected_gtf_flt = rules.run_sqanti3_filter.output[2],
        sqanti3_corrected_reason_flt = rules.run_sqanti3_filter.output[0],
        sqanti3_corrected_report_flt = rules.run_sqanti3_filter.output[1],
        pbsv_svsig = rules.run_sv_caller.output[0],
        pbsv_vcf =  rules.run_sv_caller.output[1],
        pbsv_vcf_annot = rules.run_sv_caller.output[2],
        toi_list = rules.run_toi_list.output[0],
        pbsv_result = rules.run_toi_list.output[1],
        squanti3_classification_final_results = rules.run_sqanti3_classification_results.output[0],
        squanti3_classification_filtered_out = rules.run_sqanti3_classification_results.output[1],
        dis_gene_asso = rules.run_dis_gene_asso.output[0],
        # Fusion Pipeline Outputs
        fusion_gff = rules.run_fusion_finder.output[0],
        fusion_rep = rules.run_fusion_finder.output[1],
        fusion_abundance = rules.run_fusion_finder.output[2],
        hq_isoforms_mapped = rules.run_sam_to_bam.output[0],
        hq_isoforms_mapped_bai = "results/FinalResults/{0}_hq_isoforms_mapped_hg38_sorted.bam.bai".format(CASENAME),
        fusion_product_class = rules.run_fusion_classifier.output[0],
        fusion_rep_corrected = "results/FusionClassification/lq_isoforms.fasta.fusion_corrected.fasta",
        fusion_gtf_corrected = "results/FusionClassification/lq_isoforms.fasta.fusion_corrected.gtf",
        refAnnotation_lq_isoforms = "results/FusionClassification/refAnnotation_lq_isoforms.fasta.fusion.genePred",
        final_pbfusion_mapped_bam =  rules.run_filter_pbfusion_to_bam.output[1],
        final_pbfusion_mapped_bam_bai = "results/FinalResults/{0}_final_pbfusion_mapped_hg38_sorted.bam.bai".format(CASENAME),
        fusion_annot =  rules.run_fusion_collect_info.output[0],
        fusion_annot_fltrd = "results/FusionProducts/lq_isoforms.fasta.fusion.annotated_ignored.txt",
        fusion_classification_final_results =  rules.run_fusion_events.output[0]



    output:
        collapsed_gff = rules.all.input[57],
        collapsed_rep_fa = rules.all.input[58],
        collapsed_group =  rules.all.input[59],
        collapsed_abundance =  rules.all.input[60],
        cupcake_ignored_ids = rules.all.input[61],
        collapsed_gff_fl = rules.all.input[62],
        collapsed_rep_fa_fl = rules.all.input[63],
        collapsed_abundance_fl =  rules.all.input[64],
        collapsed_gff_fl_away = rules.all.input[65],
        collapsed_rep_fa_fl_away =  rules.all.input[66],
        collapsed_abundance_fl_away =  rules.all.input[67],
        sqanti3_class = rules.all.input[69],
        sqanti3_junc = rules.all.input[70],
        sqanti3_corrected_fa = rules.all.input[71],
        sqanti3_corrected_gtf = rules.all.input[72],
        sqanti3_class_flt = rules.all.input[73],
        sqanti3_junc_flt = rules.all.input[74],
        sqanti3_corrected_fa_flt = rules.all.input[75],
        sqanti3_corrected_gtf_flt = rules.all.input[76],
        sqanti3_corrected_reason_flt = rules.all.input[77],
        sqanti3_corrected_report_flt = rules.all.input[78],
        pbsv_svsig =  rules.all.input[79],
        pbsv_vcf =   rules.all.input[80],
        pbsv_vcf_annot =  rules.all.input[81],
        toi_list = rules.all.input[82],
        pbsv_result = rules.all.input[83],
        squanti3_classification_final_results = rules.all.input[84],
        squanti3_classification_filtered_out = rules.all.input[85],
        dis_gene_asso = rules.all.input[86],
        fusion_gff =  rules.all.input[87],
        fusion_rep =  rules.all.input[88],
        fusion_abundance =  rules.all.input[89],
        hq_isoforms_mapped = rules.all.input[90],
        hq_isoforms_mapped_bai = rules.all.input[91],
        fusion_product_class = rules.all.input[92],
        fusion_rep_corrected = rules.all.input[93],
        fusion_gtf_corrected = rules.all.input[94],
        refAnnotation_lq_isoforms = rules.all.input[95],
        final_pbfusion_mapped_bam =  rules.all.input[96],
        final_pbfusion_mapped_bam_bai = rules.all.input[97],
        fusion_annot =  rules.all.input[98],
        fusion_annot_fltrd = rules.all.input[99],
        fusion_classification_final_results = rules.all.input[100]
    
    params:
        junction_shortreads_out = "results/Deliverables/Isoform/Isoform_Step2.1/{0}SJ.out.tab".format(CASENAME),
        expression_out = "results/Deliverables/Isoform/Isoform_Step2.1/cupcake.collapsed.filtered.rep.renamed.fasta.salmon"


    threads:
        1

    log:
        "log/run_deliverables.log"

    run:
        for no in range(len(output)):
            orig_file = input[no]
            dest_file = output[no]
            shell(\
                    """
                    ln -s $PWD/{orig_file} {dest_file}
                    """ \
            )
        
        # Check if short reads are processed, then link STAR Junc and Salmon results to Isoform_Step2.1
        step2_1_dir = os.path.dirname(params.junction_shortreads_out)
        shell("mkdir {step2_1_dir}")
        if  config["ILLUMINASHORTREADS"]["ill_fastq_R1"]:
            junction_shortreads = rules.run_short_reads_genome_alignment.params.output_prefix + "SJ.out.tab"
            expression = rules.run_transcriptome_alignment.params.output_dir
            shell(\
                    """  
                        ln -s $PWD/{junction_shortreads} $PWD/{params.junction_shortreads_out} 
                        ln -s $PWD/{expression} $PWD/{params.expression_out} 
                    """ \
                )
        else:
            shell(\
                    """
                        touch {params.junction_shortreads_out} 
                        mkdir {params.expression_out} 
                    """ \
                )



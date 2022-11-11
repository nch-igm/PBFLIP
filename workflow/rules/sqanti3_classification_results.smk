
TX2G = config["TX2G"]



def associate_ccs_read_counts( squanti3_output, cupcake_collapsed_stat, event_type = "isoform"):
    #case_parameters = read_case_config( wd )
    ccs_counter = []
    if event_type == 'isoform':
        sep_squanti  = '\t'
        sep_stat = "\t"
        left_on = 'isoform'
    elif event_type == 'fusion':
        sep_squanti  = ','
        sep_stat = "\t"
        left_on = 'UniqueID'

    if os.path.exists( squanti3_output ) and os.path.exists (cupcake_collapsed_stat):
        squant3_class_output_pd_df = pd.read_csv(squanti3_output, sep = sep_squanti)
        cupcake_collapsed_stat_pd_df = pd.read_csv(cupcake_collapsed_stat, sep = sep_stat)
    else:
        raise Exception( "File not found {0} or {1}".format( squanti3_output, cupcake_collapsed_stat))
    
    # Select full length reads only 
    is_fl = cupcake_collapsed_stat_pd_df["is_fl"] == "Y"
    cupcake_collapsed_stat_pd_df_fl = cupcake_collapsed_stat_pd_df[is_fl]
    cupcake_collapsed_stat_pd_df_fl_dedup = cupcake_collapsed_stat_pd_df_fl.groupby( ['pbid'])['id'].apply(', '.join).reset_index()
    result = pd.merge(squant3_class_output_pd_df, cupcake_collapsed_stat_pd_df_fl_dedup, left_on = left_on, right_on = 'pbid')
    for ccs in result['id'].values:
        ccs_counter.append(len(ccs.split(',')))


    result["CCS_Counts"]= ccs_counter
    return result



def geneid_to_gene_name():
    # Reading transcipt ids ensemble ids and gene names file
    #case_parameters = read_case_config()
    tr2g = TX2G
    if os.path.exists(  tr2g ):
        tr2g_df = pd.read_table( tr2g)
    else:
        raise Exception(" No file found: {0}". format( tr2g ))

    g2genename = tr2g_df[["gene", "gene_name"]]
    g2genename  = g2genename.drop_duplicates()

    return g2genename


def associate_gene_names_isoforms( isoform_detection_results_input):
    gid2g={}
    geneid_to_gene_name_df =  geneid_to_gene_name( )
    for gene, gene_name in zip(geneid_to_gene_name_df['gene'], geneid_to_gene_name_df['gene_name']):
        gid2g[gene_name]=gene
    
    associated_gene_id = []

    for geneid in isoform_detection_results_input['associated_gene']:
        if type( geneid ) == str:
            if  not geneid.startswith("novelGene"):
                geneid_split = geneid.split('_')
                if len(geneid_split)==1:
                    if geneid in gid2g.keys():
                        associated_gene_id.append(gid2g[geneid])
                    else:
                        associated_gene_id.append(geneid)
                else:
                    g1 = geneid_split[0]
                    g2 = geneid_split[1]
                    if (g1 in gid2g.keys()) and (g2 in gid2g.keys()):
                        gene_pair = gid2g[g1] +'_' + gid2g[g2]
                    elif (g1 +"_"+ g2 in gid2g.keys()):
                        gene_pair = gid2g[g1 +"_"+ g2]

                    associated_gene_id.append(gene_pair )
            else:
                associated_gene_id.append(geneid )
        else:
                 associated_gene_id.append("")
    
    isoform_detection_results_input['associated_gene_id'] = associated_gene_id 

    return isoform_detection_results_input

def filtered_out_lite_classification_isoforms(unfiltered_classification_file_pd , filtered_classification_file_pd ):
    #unfiltered_classification_file_pd = pd.read_csv(unfiltered_class_file_path, sep ="\t")
    #filtered_classification_file_pd = pd.read_csv(filtered_class_file_path, sep ="\t")
    filtered_out_isoforms = unfiltered_classification_file_pd[~unfiltered_classification_file_pd.isoform.isin(filtered_classification_file_pd.isoform)]
    return filtered_out_isoforms 


rule run_sqanti3_classification_results:
    input:
        squanti3_output_filtered = rules.run_sqanti3_filter.output[4],
        cupcake_collapsed_stat = rules.run_get_abundance.output[0],
        squanti3_output_unfiltered = rules.run_sqanti3_classification.output[0]

    
    output:
        rules.all.input[29],
        rules.all.input[30]

    params:
        pythonpath = config["LIBPATHS"]["python_lib"]

    threads:
        1

    log:
        "log/run_sqanti3_classification_results.log"

    threads:
        1
    
    run:

        isoform_detection_results = associate_ccs_read_counts( \
                                         input.squanti3_output_filtered, 
                                         input.cupcake_collapsed_stat, 
                                         event_type = "isoform" \
                                         )
        isoform_detection_results_unflitered = associate_ccs_read_counts( \
                                                     input.squanti3_output_unfiltered , 
                                                     input.cupcake_collapsed_stat, 
                                                     event_type = "isoform" \
                                                     )
        isoform_detection_results_gene_names_added = associate_gene_names_isoforms( \
                                                                isoform_detection_results )

        isoform_detection_results_unflitered_gene_names_added = associate_gene_names_isoforms( \
                                                                            isoform_detection_results_unflitered \
                                                                             )

        isoform_detection_filtered_out_results =  filtered_out_lite_classification_isoforms( \
                                                            isoform_detection_results_unflitered_gene_names_added , isoform_detection_results_gene_names_added \
                                                            )
        
        isoform_detection_results_gene_names_added.to_csv( \
                                                                      output[0], 
                                                                      sep ="\t", 
                                                                      index =  False \
                                                                    )
        isoform_detection_filtered_out_results.to_csv( \
                                                                output[1], 
                                                                sep ="\t", 
                                                                index =  False \
                                                            )

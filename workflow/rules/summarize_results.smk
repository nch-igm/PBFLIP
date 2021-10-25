
def split_list( list_to_split, data_type ):
    if any(isinstance(i, list) for i in list_to_split):
        concat_list = list(itertools.chain.from_iterable(list_to_split ))
    else:
        concat_list = list_to_split

    split_list = []
    if data_type == "pbsv":
        for item in concat_list:
            split_list.extend((item.split(",")))

    elif data_type == "sqanti3":
        for item in concat_list:
            sub_class, junct = item.split(",")
            if junct == "canonical":
               split_list.append( "canonical"+":"+sub_class)
               split_list.append( junct )
            elif junct == "non_canonical":
                split_list.append( "non_canonical"+":"+sub_class)
                split_list.append( junct )
            elif junct == "":
                split_list.append( "no_classification"+":" + sub_class)
                split_list.append( "no_classification" )

            else:
                split_list.append( sub_class)
    
    return Counter( split_list)   




def sqanti3_output_filter( squanti_output):
    gene_dict = {}
    for index in range( len(squanti_output)):
        if index == 0:
            names = squanti_output[ index ].split("\t")
        else:
            #gene = squanti_output[ index ].split("\t")[49].replace("\n", "")
            gene = squanti_output[ index ].split("\t")[6].replace("\n", "")
            if gene in gene_dict.keys():
                gene_dict[ gene ].append( squanti_output[ index ].split("\t")[5] +","+squanti_output[ index ].split("\t")[16])
            else:
                gene_dict [ gene ] = [squanti_output[ index ].split("\t")[5] +"," +squanti_output[ index ].split("\t")[16]]
    return ( gene_dict)

"""
The following function is taken from here:
https://towardsdatascience.com/reordering-pandas-dataframe-columns-thumbs-down-on-standard-solutions-1ff0bc2941d5
"""
def movecol(df, cols_to_move=[], ref_col=1, place='After', table_type = 'sqanti3'):
    
    cols = df.columns.tolist()
    if place == 'After':
        seg1 = cols[:ref_col + 1]
        seg2 = cols_to_move
    if place == 'Before' and table_type == 'sqanti3':
        seg1 = cols[:ref_col]
        seg2 = cols_to_move + [cols[ref_col]]
    elif place == 'Before' and table_type == 'pbsv':
        seg1 = cols[:ref_col]
        seg2 = cols_to_move        
    
    seg1 = [i for i in seg1 if i not in seg2]
    seg3 = [i for i in cols if i not in seg1 + seg2]
    
    return(df[seg1 + seg2 + seg3])

def Merge(dict1, dict2):
    return(dict2.update(dict1))

def run_results_summarization ( sample_name, final_results_file, final_results_filtered_out_file, pbsv_results_file, summary_result_files_sqanti3, summary_result_files_pbsv, log):
    
    #logging.basicConfig(filename=log, level=logging.DEBUG)
    with open(final_results_file, "r") as rf1:
        data1 = rf1.readlines()
    with open(final_results_filtered_out_file, "r") as rf2:
        data2 = rf2.readlines()

    data1 += data2


    gene_dict = sqanti3_output_filter( data1 )
    #gene_dict2 = sqanti3_output_filter( data2 )
    #logging.debug(pprint.pformat(gene_dict))

    final_result_gene_dict = {}
    SVTYPE = []
    IMPACT = []
    with open ( pbsv_results_file, "r") as pbsv:
        header= pbsv.readline()
        for body in pbsv:
            split_line = body.split ("\t")
            svtype = split_line [4]
            gene = split_line [6].split(";")[0]
            impact = split_line [7].split(";")[0]
            effect  = split_line [8].split(";")[0]
            #print ( svtype, gene, impact, effect)
            temp_dict ={}
            if gene in gene_dict.keys():
                #print ( gene, svtype, impact, effect)
                squanti_class = gene_dict [ gene]
                if  gene in final_result_gene_dict.keys():
                    final_result_gene_dict[ gene ].append(final_result_gene_dict[ gene ] [0]["pbsv"].append([svtype, impact, effect] ))
                    SVTYPE.append( svtype)
                    IMPACT.append( impact )
                else:

                    temp_dict = {"pbsv":[[svtype, impact, effect]], \
                                        "sqanti3":squanti_class
                                        }
                    final_result_gene_dict[ gene ] = [ temp_dict  ]
                    SVTYPE.append( svtype)
                    IMPACT.append( impact )

            else:
                #logging.debug( "No match" + gene)
                pass



    final_combined_data_dict  = {}
    for key in final_result_gene_dict.keys():
        #print (key)
        #logging.debug(pprint.pformat((split_list(  final_result_gene_dict [ key] [0]['pbsv'], "pbsv"))))
        pbsv =  dict(split_list(  final_result_gene_dict [ key] [0]['pbsv'], "pbsv"))
        #logging.debug(pprint.pformat(pprint.pprint ( split_list(  final_result_gene_dict [ key] [0]['sqanti3'], "sqanti3")))
        squant3 = dict(split_list(  final_result_gene_dict [ key] [0]['sqanti3'], "sqanti3"))
        #final = Merge( pbsv, squant3)
        #print ( squant3)
        final_combined_data_dict [ key ] = {'squant3':squant3,  'pbsv': pbsv}

    #logging.debug("final_combined_data_dict" )
    #logging.debug(pprint.pformat( final_combined_data_dict ))
    #pprint.pprint ( final_combined_data_dict )



    l2_sqanti3 = list()

    for name, dcts in final_combined_data_dict.items():
        sqant3_dict = dict(dcts ['squant3'])
        sqant3_dict["Gene"] = name
        l2_sqanti3.append( sqant3_dict )

    df_sqanti3= pd.DataFrame(l2_sqanti3)
    df_sqanti3.rename(columns = {"" : "No Classification"}, inplace = True)
    df_sqanti3.fillna(0, inplace=True, downcast='infer')
    df_sqanti3_final = movecol(df_sqanti3, 
                cols_to_move=['Gene','canonical', 'non_canonical', "no_classification"], 
                ref_col=0,
                place='Before')


    # Writing final sumarized output
    df_sqanti3_final.to_csv(summary_result_files_sqanti3 , sep = "\t", index=False)

    l2_pbsv = list()
    SVTYPE = list( set( SVTYPE ) )
    IMPACT = list ( set ( IMPACT ))
    for name, dcts in final_combined_data_dict.items():
        pbsv_dict = dict(dcts ['pbsv'])
        pbsv_dict["Gene"] = name
        l2_pbsv.append( pbsv_dict )

    df_pbsv= pd.DataFrame(l2_pbsv)
    df_pbsv.fillna(0, inplace=True, downcast='infer')
    df_pbsv_final = movecol(df_pbsv, 
                cols_to_move= ['Gene'] + SVTYPE +IMPACT , 
                ref_col=0,
                place='Before', table_type = "pbsv")



    #df_pbsv_final.to_csv("test_pbsv.csv", index=False)


    annotation = []
    for i, row in df_pbsv_final.iterrows():
        concat_str = ""
        for cl in df_pbsv_final.columns:
            if cl not in  ['Gene'] + SVTYPE +IMPACT:
                val = df_pbsv_final.at[i,cl]
                if val == 0:
                    pass
                else:
                    concat_str = concat_str+ cl+":"+ str( val ) +"|"
        #logging.debug (concat_str  )
                    
        annotation.append( concat_str )

    df_pbsv_final["Annotation"] = annotation

    df_pbsv_final_filtered = df_pbsv_final[ ['Gene'] + SVTYPE +IMPACT + ["Annotation"]]

    df_pbsv_final_filtered.to_csv(summary_result_files_pbsv , sep = '\t', index=False)




rule run_summarize_results:
    input:    
        final_results_file = rules.run_sqanti3_classification_results.output[0],
        final_results_filtered_out_file = rules.run_sqanti3_classification_results.output[1],
        pbsv_vcf_annot_tsv  = rules.run_sv_caller.output.pbsv_vcf_annot_tsv

    output:
        summary_result_files_sqanti3 =rules.all.input[45],
        summary_result_files_pbsv = rules.all.input[46]


    
    params:
        sample_name = config["CASENAME"]


    threads:
        1

    log:
        "log/run_summarize_results.log"
        

    run:
        
        run_results_summarization ( \
                                            params.sample_name, 
                                            input.final_results_file, 
                                            input.final_results_filtered_out_file, 
                                            input.pbsv_vcf_annot_tsv, 
                                            output.summary_result_files_sqanti3, 
                                            output.summary_result_files_pbsv,
                                            log[0] \
                                        )
        
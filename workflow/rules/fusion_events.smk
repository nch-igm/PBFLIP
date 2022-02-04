include: "sqanti3_classification_results.smk"

FALSEPOSDB = ["Bodymap2", 	"HPA", 	"Non_Tumor_Cells", 	"Babiceanu_Dataset",  "Banned_Dataset",  "GTEx"]

def count_no_db( db_list ):
    item_length = []
    for item in db_list:
        item_length.append(len(item.split(",") ))
    
    return item_length

def check_in_fusionhub_db( fusionhub_global_summary, fusion_events, gene, false_pos_db, false_pos_db_perc ):
    no_databases_false_pos = 0
    false_pos_db_perc_value = 0
    found_status = fusion_events.str.contains(gene, regex=False)
    #print( gene)
    item_contain_df = fusionhub_global_summary[ found_status ]
    databases_names_col = item_contain_df.columns[item_contain_df.isin(['+']).any()]
    
    for db in databases_names_col:
        if db in FALSEPOSDB:
            no_databases_false_pos +=1


    no_db = len( databases_names_col )
    
    if no_databases_false_pos != 0:
        false_pos_db_perc_value = round((no_databases_false_pos/no_db) * 100, 2)

    databases_names_col = ",".join(databases_names_col )
    return databases_names_col , no_db, no_databases_false_pos,  false_pos_db_perc_value


def read_fusion_table( fusion_candid_file):
    df = pd.DataFrame()
    
    with open(fusion_candid_file, 'r') as f:
        header = f.readline().split()

        for line in f:
            split_items = line.strip().split('\t')
            main_col = split_items[0:24]
            #ccs_col = ",".join(split_items[18:])
            #main_col.append(ccs_col)
            df = pd.concat( [df, pd.DataFrame([tuple(main_col)])], ignore_index=True )
    #df = df.iloc[1:]
    #df.columns = ["ID.1" ,   "CHR1" ,"STR1" ,"LENGTH1",  "EXONS1",   "CLASS1",   "GENE1" ,"GRCh38_pos1" , "ID.2" ,"CHR2" ,"STR2" ,"LENGTH2",  "EXONS2",    "CLASS2",    "GENE2" ,"GRCh38_pos2" , "COUNT", "CCS"]
    if df.empty:
        return df
    df.columns = header
    df.drop('pbid',axis='columns', inplace=True)
    df.drop('Sequence',axis='columns', inplace=True)
    df=df.rename(columns={"id": "CCS"})
    return df

def assign_db_status(fusionhub_global_summary, fusion_events, fusion_db_list, no_fusion_db, both_found, pair, false_pos_db, false_pos_db_perc):

    db_string, no_db, no_databases_false_pos, false_pos_db_perc_value = check_in_fusionhub_db( fusionhub_global_summary, fusion_events,  pair, false_pos_db, false_pos_db_perc )
    if db_string and no_db:
        fusion_db_list.append( db_string )
        no_fusion_db.append( no_db )
        both_found.append("Pair Found")
        false_pos_db.append( no_databases_false_pos )
        false_pos_db_perc.append( false_pos_db_perc_value )
    else:
        db_string, no_db, no_databases_false_pos, false_pos_db_perc_value = check_in_fusionhub_db( fusionhub_global_summary, fusion_events,  pair.split("--") [0], false_pos_db , false_pos_db_perc)
        if db_string and no_db:
                #print( db_string, no_db)
                fusion_db_list.append( db_string )
                no_fusion_db.append( no_db )
                both_found.append("Gene1 Found")
                false_pos_db.append( no_databases_false_pos )
                false_pos_db_perc.append( false_pos_db_perc_value )
        else:
            db_string, no_db, no_databases_false_pos, false_pos_db_perc_value = check_in_fusionhub_db( fusionhub_global_summary, fusion_events,  pair.split("--") [1], false_pos_db, false_pos_db_perc)
            if db_string and no_db:
                fusion_db_list.append( db_string )
                no_fusion_db.append( no_db )
                both_found.append("Gene2 Found")
                false_pos_db.append( no_databases_false_pos )
                false_pos_db_perc.append( false_pos_db_perc_value )
            else:
                fusion_db_list.append( "" )
                no_fusion_db.append( 0 )
                both_found.append("Not Found")
                false_pos_db.append( no_databases_false_pos )
                false_pos_db_perc.append( false_pos_db_perc_value )
    return fusion_db_list,  no_fusion_db, both_found, false_pos_db, false_pos_db_perc
                
 

rule run_fusion_events:
    input:
        fusion_candidates_file = rules.run_fusion_classification_results.output[0]
    output:
        final_fusion_class_fusionhub = rules.all.input[54]
    params:
        fusionhub_file = config["FUSIONHUBDB"]
    threads:
        1
    log:
        "log/run_fusion_events.log"
    run:
        geneid_2_gene_name =  geneid_to_gene_name( )
        read_fusion_candidates  = read_fusion_table ( \
                                        input.fusion_candidates_file\
                                    )
        if not read_fusion_candidates.empty:
            fusion_name_left_gene = read_fusion_candidates.LeftGeneName.str.split( \
                                        "_", expand=False \
                                        )
            fusion_name_right_gene = read_fusion_candidates.RightGeneName.str.split( \
                                        "_", expand=False \
                                        )
            if os.path.exists( params.fusionhub_file):
                fusionhub_global_summary =pd.read_table(params.fusionhub_file)
            else:
                read_fusion_candidates.to_csv( \
                                                    read_fusion_candidates , 
                                                    sep ="\t", 
                                                    index =  False \
                                                )
                raise Exception( " Fusion Hub database file doesnt exist...")
            
            gene1 = []
            gene2 = []
            gene3 = []
            gene4 = []
            fusion_db_list = []
            no_fusion_db = []
            both_found = []
            false_pos_db = []
            false_pos_db_perc = []
            gid2g={}
            fusion_events = fusionhub_global_summary [ 'Fusion_gene']

            for gene, gene_name in zip(geneid_2_gene_name['gene'], geneid_2_gene_name['gene_name']):
                gid2g[gene]=gene_name

            g1_c = 0
            g2_c = 0
            for lg, rg  in zip(fusion_name_left_gene, fusion_name_right_gene):
                if len(lg) == 1 and len(rg) == 1:
                    g1 = lg[0].replace("PAR", "").replace("Y", "")
                    gene1.append(gid2g[g1])
                    g2 = rg[0].replace("PAR", "").replace("Y", "")
                    gene2.append(gid2g[g2])
                    pair = gid2g[g1] + '--' + gid2g[g2]

                    # Check the databases
                    #db_string, no_db = check_in_fusionhub_db( fusionhub_global_summary, fusion_events,  pair )
                    fusion_db_list,  no_fusion_db, both_found, no_databases_false_pos, false_pos_db_perc = assign_db_status(fusionhub_global_summary, fusion_events, fusion_db_list, no_fusion_db, both_found, pair,false_pos_db, false_pos_db_perc )

                elif len(lg) > 1 and len(rg) == 1:
                    g1tmp = []
                    for splitgene in lg:
                        if  splitgene not in ["PAR", "Y"]:
                            g1tmp.append(gid2g [splitgene ])
                    gene1.append('_'.join(g1tmp))
                    g2 = rg[0].replace("PAR", "").replace("Y", "")
                    gene2.append(gid2g[g2])
                    pair = '_'.join(g1tmp) + '--' + gid2g[g2]

                    fusion_db_list,  no_fusion_db, both_found, no_databases_false_pos, false_pos_db_perc = assign_db_status(fusionhub_global_summary, fusion_events, fusion_db_list, no_fusion_db, both_found, pair, false_pos_db, false_pos_db_perc)
                elif len(rg) > 1 and  len(lg) == 1:
                    g1tmp = []
                    for splitgene in rg:
                        if  splitgene not in ["PAR", "Y"]:
                            g1tmp.append(gid2g [splitgene ])
                    gene2.append('_'.join(g1tmp))
                    g1 = lg[0].replace("PAR", "").replace("Y", "")
                    gene1.append(gid2g[g1])
                    pair = gid2g[g1] +  '--'  + '_'.join(g1tmp) 

                    fusion_db_list,  no_fusion_db, both_found, no_databases_false_pos, false_pos_db_perc = assign_db_status(fusionhub_global_summary, fusion_events, fusion_db_list, no_fusion_db, both_found, pair, false_pos_db, false_pos_db_perc)
                elif len(lg) > 1 and len(rg) >= 1:
                    g1tmp = []
                    for splitgene in lg:
                        if  splitgene not in ["PAR", "Y"]:
                            g1tmp.append(gid2g [splitgene ])
                    gene1.append('_'.join(g1tmp))
                    g2tmp = []
                    for splitgene in rg:
                        if  splitgene not in ["PAR", "Y"]:
                            g2tmp.append(gid2g [splitgene ])
                    gene2.append('_'.join(g2tmp))
                    pair = '_'.join(g1tmp) +  '--'  + '_'.join(g2tmp) 

                    fusion_db_list,  no_fusion_db, both_found, no_databases_false_pos, false_pos_db_perc = assign_db_status(fusionhub_global_summary, fusion_events, fusion_db_list, no_fusion_db, both_found, pair, false_pos_db, false_pos_db_perc)
            
            read_fusion_candidates [ 'LeftGeneName'] = gene1
            #read_fusion_candidates [ 'GENE12'] = gene2
            read_fusion_candidates [ 'RightGeneName'] = gene2
            FusionName2 = []
            read_fusion_candidates [ 'FusionName2'] = [i +"--" + j for i, j in zip(gene1, gene2)]
            read_fusion_candidates['Found Status'] = both_found
            read_fusion_candidates['No_DB'] =  no_fusion_db
            read_fusion_candidates['No_None_Disease_DB'] =  no_databases_false_pos
            read_fusion_candidates['None_Disease_DB%'] =  false_pos_db_perc
            read_fusion_candidates['DB'] =  fusion_db_list

            # Writing final table to same FinalResults

            read_fusion_candidates.to_csv(output.final_fusion_class_fusionhub , sep ="\t", index =  False )
        else:
            # Writing a empty dataframe to final results file
            with open( output.final_fusion_class_fusionhub, "w") as wf:
                pass
            #pd.DataFrame({}).to_csv(output.final_fusion_class_fusionhub , sep ="\t", index =  False )

    

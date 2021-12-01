def get_no_txs_support ( db_item):
    return len(db_item )

def final_results_headers( wf ):
    wf.write("CHR_POS" +"\t" + "igv_url"+ "\t" +"HQ_Isoform_ID" +"\t" + "Event_Type" +"\t" + "No_HQ_Isoform_Support"+"\t" + "isoform"+"\t"+"chrom"+"\t"+"strand"+"\t"+"length"+"\t"+"exons"+"\t"+"structural_category"+"\t"+"associated_gene"+"\t"+"associated_transcript"+"\t"+"ref_length"+"\t"+"ref_exons"+"\t"+"diff_to_TSS"+"\t"+"diff_to_TTS"+"\t"+"diff_to_gene_TSS"+"\t"+"diff_to_gene_TTS"+"\t"+"subcategory"+"\t"+"RTS_stage"+"\t"+"all_canonical"+"\t"+"min_sample_cov"+"\t"+"min_cov"+"\t"+"min_cov_pos"+"\t"+"sd_cov"+"\t"+"FL"+"\t"+"n_indels"+"\t"+"n_indels_junc"+"\t"+"bite"+"\t"+"iso_exp"+"\t"+"gene_exp"+"\t"+"ratio_exp"+"\t"+"FSM_class"+"\t"+"coding"+"\t"+"ORF_length"+"\t"+"CDS_length"+"\t"+"CDS_start"+"\t"+"CDS_end"+"\t"+"CDS_genomic_start"+"\t"+"CDS_genomic_end"+"\t"+"predicted_NMD"+"\t"+"perc_A_downstream_TTS"+"\t"+"seq_A_downstream_TTS"+"\t"+"dist_to_cage_peak"+"\t"+"within_cage_peak"+"\t"+"dist_to_polya_site"+"\t"+"within_polya_site"+"\t"+"polyA_motif"+"\t"+"polyA_dist"+"\t"+"ORF_seq"+"\t"+"pbid"+"\t"+"id"+"\t"+"CCS_Counts"+"\t"+"associated_gene_name\n" )

def separate_items (db_item ):
    unnamedSample_HQ_l = []
    event_type_l = []
    pb_id_l = []
    for item in db_item:
        unnamedSample_HQ, event_type,  pb_id = item.split(",")
        unnamedSample_HQ_l.append( unnamedSample_HQ)
        event_type_l.append( event_type)
        pb_id_l.append( pb_id )
    return  ",".join(unnamedSample_HQ_l) + "\t" + ",".join( list(set(event_type_l))) + "\t" +  ",".join(pb_id_l)


def get_iso_id ( rci_db_ids_holder , rci_db, id):
    pb_id_holder  = []
    if id in rci_db_ids_holder:
            pb_id_holder.append(rci_db[id])
    return pb_id_holder 


def toi_list_generator( sample_name, collapsed_group, svsig_file, final_results_file, final_results_filtered_out_file, toi_list, pbsv_results, igvurl, log ):

    with open( collapsed_group, "r") as rci:
        rci_db = {}
        for line in rci:
            split_line = line.split()
            key = split_line[0]
            iso_header =  split_line[1].split()
            for hd in  iso_header:
                if hd in rci_db.keys():
                    rci_db[ hd ].append(key)
                else:
                    rci_db[ hd ] = [key]

    rci_db_ids_holder = []
    for id_merge in rci_db.keys ():
        id_split = [id_split for id_split in id_merge.split(',')]
        rci_db_ids_holder.append( "".join(id_split ) )
    list_df = []
    event_dict = {}
    with gzip.open( svsig_file, "r") as svsig:
        for line in svsig:
            line = line.decode('utf-8')
            #if not line.startswith("#") and not line.startswith("c") and not line.startswith("s"):
            if not line.startswith("#") and not line.startswith("c"):
                split_line = line.split()
                #print (split_line)
                list_df.append(split_line)
                
                if split_line[0] == 'b':
                    
                    key = split_line[2] +":"+split_line[3]
                    items = split_line[9] + "," +split_line[0]
                    iso_id = get_iso_id ( rci_db_ids_holder, rci_db,split_line[9])

                    if iso_id:
                        items = items + ","+",".join(iso_id[0])
                    else:
                        items = items + ","+"No PBID"
                else:
                    key = split_line[2] +":"+split_line[3]
                    items = split_line[4] + "," +split_line[0]
                    iso_id = get_iso_id ( rci_db_ids_holder, rci_db,split_line[4])
                    if iso_id:
                        items = items +"," +",".join(iso_id[0])
                    else:
                        items = items + ","+"No PBID"


                if key in event_dict.keys():
                    event_dict [ key ].append(items)
                else:
                    event_dict [ key ] = [items]


    with open( pbsv_results, "w") as wf:
        wf.write("CHR_POS" +"\t" + "HQ_Isoform_ID" +"\t" + "Event_Type" +"\t" + "PBIsoform_ID" +"\t" + "No_HQ_Isoform_Support" "\n")
        for key in event_dict.keys():
            item = event_dict[key]
            line = key +"\t" + separate_items (item ) + "\t" + str(get_no_txs_support ( item))
            #print ( line)
            wf.write( line +"\n")

    with open(final_results_file, "r") as rf1:
        data1 = rf1.readlines()
    with open(final_results_filtered_out_file, "r") as rf2:
        data2 = rf2.readlines()

    data1 += data2 

    temp_dict = {} # pbids in event_dict
    temp_2 = [] # pbids in data1

    for data_item in data1:
        pbid = data_item.split("\t")[0]
        temp_2.append(pbid)
        temp_dict[pbid] = data_item

    with open(toi_list, "w") as f:
        final_results_headers( f )
        for key in event_dict.keys():
            items = event_dict[key]
            line_key = key
            for item in items:
                unnamedSample_HQ, event_type,  pb_id = item.split(",")
                #temp_1.append(pb_id)
                #print(line_key, unnamedSample_HQ, event_type,pb_id )
 
                if pb_id in temp_2:
                    data_item = temp_dict[ pb_id ]
                    #print ( line_key, unnamedSample_HQ, event_type, pb_id, data_item)
                    igv_url = '=HYPERLINK("' + IGVURL + '={location}")'.format( location = line_key )
                    wr_line = "\t".join( [line_key, igv_url,  unnamedSample_HQ, event_type, str(get_no_txs_support ( items))]) +"\t"+data_item
                    f.write( wr_line )
    print ("Done")

rule run_toi_list:
    input:
        collapsed_group     = rules.run_collapse.output[2],
        pbsv_svsig           = rules.run_sv_caller.output[0],
        final_results_file   = rules.run_sqanti3_classification_results.output[0],
        final_results_filtered_out_file = rules.run_sqanti3_classification_results.output[1]


    output:
        toi_list = rules.all.input[43],
        pbsv_result = rules.all.input[44]

    
    params:
        igv_url = 'http://localhost:60151/goto?locus='


    threads:
        1

    log:
        "log/run_toi_list.log"

    run:
        toi_list_generator( \
                                CASENAME, 
                                input.collapsed_group, 
                                input.pbsv_svsig, 
                                input.final_results_file, 
                                input.final_results_filtered_out_file,
                                output[0],
                                output[1],
                                params.igv_url,
                                log[0] \
                            )
        

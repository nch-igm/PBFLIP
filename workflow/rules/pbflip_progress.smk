class collectTranscriptsInfo:
    def __init__(self, fastx):
        def merge_header_space(identifier):
            identifier = identifier.replace(" ", "_")
            return identifier
        if (fastx.endswith(".fastq") or fastx.endswith(".fq")) :
            self.record_dict = SeqIO.index(fastx, "fastq", key_function=merge_header_space)
        elif (fastx.endswith(".fasta") or fastx.endswith(".fa")):
            self.record_dict = SeqIO.index(fastx, "fasta")
        else:
            print("Only fastq or fasta are supported")
            exit()
    


    def return_number_records( self):
        return len(self.record_dict)
    
    def return_records_ids( self):
        return self.record_dict.keys()

    def return_records( self):
        return self.record_dict.items()

def pb_progress( collapsed_obj_fa, collapsed_filtered_obj_fa, classification_filtered_obj_fa,  collapsed_5prim_filtered_obj_fa, collapsed_obj_info_collector):
    collapsed_obj =  collectTranscriptsInfo(collapsed_obj_fa)
    collapsed_filtered_obj =  collectTranscriptsInfo(collapsed_filtered_obj_fa)
    collapsed_5prim_filtered_obj =  collectTranscriptsInfo(collapsed_5prim_filtered_obj_fa)
    classification_filtered_obj =  collectTranscriptsInfo(classification_filtered_obj_fa)


    collapsed_obj_info_collector_ds = set(collapsed_obj_info_collector)
    collapsed_obj_ids = set(collapsed_obj.return_records_ids())
    collapsed_filtered_obj_ids = set(collapsed_filtered_obj.return_records_ids())
    collapsed_5prim_filtered_obj_ids = set(collapsed_5prim_filtered_obj.return_records_ids())
    classification_filtered_obj_ids = set(classification_filtered_obj.return_records_ids())

    
    collapsed_dict = {}

    def match_collapsed_ids( dict, set_get_compared, flag):
        for fid in set_get_compared:
            try:
                id, collapsed_iso =  fid.split("|",1)
            except:
                id = fid
            if id in dict.keys():
                dict[id].append(collapsed_iso)
            elif flag == 1:
                dict[id] = [collapsed_iso]
            else:
                dict[id].append(collapsed_iso)
    
        return dict

    def match_ids( dict, set_get_compared):
        for id in dict.keys():
            if id in set_get_compared:
                #print("Passed")
                dict[id].append("Passed")
            else:
                #print( "NtP")
                dict[id].append("Not_Passed")
                #dict[id].append("Not_Passed")
    
        return dict

    def split_header( set_with_headers):
        split_headers = []
        for fid in set_with_headers:
            try:
                id, collapsed_iso =  fid.split("|",1)
                split_headers.append(id)
            except:
                id = fid
                split_headers.append(id)
        
        return split_headers 

    collapsed_dict = match_collapsed_ids( collapsed_dict, collapsed_obj_ids, 1)
    collapsed_dict = match_collapsed_ids( collapsed_dict, collapsed_obj_info_collector_ds, 2)

    split_headers = split_header( collapsed_filtered_obj_ids)
    collapsed_dict = match_ids( collapsed_dict, split_headers)
    split_headers = split_header( collapsed_5prim_filtered_obj_ids)
    collapsed_dict = match_ids( collapsed_dict, split_headers)
    split_headers = split_header(classification_filtered_obj_ids)
    collapsed_dict = match_ids( collapsed_dict,split_headers )

    collapsed_dict_df = pd.DataFrame.from_dict(collapsed_dict, orient='index',  columns=[ 'Collapsed_Alignment_Status', 'Collapsed_Status', 'Collapsed_Abundance_Filter', 'Collapsed_Abundance_5primFilter', 'Classification_Filter'])
    collapsed_dict_df.index.name = "PBID"
    
    return collapsed_dict_df

rule run_pbflip_progress:
    input:
        collapsed_group = "results/Collapsed_Isoforms/cupcake.collapsed.group.txt",
        collapsed_obj_fa = rules.run_collapse.output[1],
        collapsed_filtered_obj_fa = "results/Collapsed_Isoforms/cupcake.collapsed.min_fl_2.rep.fa",
        collapsed_5prim_filtered_obj_fa = rules.run_filter_away_subset.output[1],
        classification_filtered_obj_fa = rules.run_sqanti3_filter.output[5]
    output:
        rules.all.input[37]

    threads:
        1

    log:
        "log/run_pbflip_progress.log"
    run:
        collapsed_obj_info_collector = []

        with open(input.collapsed_group, "r") as fr:
            for line in fr:
                collapsed_obj_info_collector.append(line.replace("\t", "|").replace(",", "|").replace("\n", ""))
        
        collapsed_dict_df = pb_progress( \
                                                input.collapsed_obj_fa, 
                                                input.collapsed_filtered_obj_fa, 
                                                input.classification_filtered_obj_fa,  
                                                input.collapsed_5prim_filtered_obj_fa, 
                                                collapsed_obj_info_collector \
                                            )
        collapsed_dict_df.to_csv( output[0], sep ="\t")



DISGENET = config["DISGENET"]

def dis_gene_asso(isoform_final_result):
    """
        Gets DisGeNET file and associates the infor with sqanti3_final_results.tsv
    """

    dis_gene_network_df = pd.read_csv( \
                                             DISGENET, 
                                             sep = '\t'
                                        )
    isoform_final_result_df = pd.read_csv( \
                                                 isoform_final_result, 
                                                 sep = '\t'
                                                )

    dis_gene_network_df_dis_name = dis_gene_network_df.groupby( \
                                                                ["geneSymbol"])["diseaseName"].apply(', '.join).reset_index()
    dis_gene_network_df_dis_type = dis_gene_network_df.groupby(["geneSymbol"])["diseaseType"].apply(', '.join).reset_index()

    result = pd.merge(dis_gene_network_df_dis_name, dis_gene_network_df_dis_type, on='geneSymbol')

    intersect_result = set(result['geneSymbol']).intersection(set(isoform_final_result_df['associated_gene']))

    gene_dis_ass_results = result[result['geneSymbol'].isin(intersect_result)]

    return gene_dis_ass_results

    

rule run_dis_gene_asso:
    input:
        rules.run_sqanti3_classification_results.output[0]
    output:
        rules.all.input[36]

    threads:
        1

    log:
        "log/run_dis_gene_asso.log"
    run:
        gene_dis_ass_results = dis_gene_asso(input[0])
        gene_dis_ass_results.to_csv( output[0], sep ="\t", index =  False )

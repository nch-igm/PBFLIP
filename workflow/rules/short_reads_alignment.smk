import ntpath
import pandas as pd

def clean_salmon_output( quant_file):
    read_quant_sf = pd.read_table(quant_file , sep ="\t")
    col_name = ["target_id", "length",  "eff_length",  "est_counts",  "tpm"]
    temp = read_quant_sf [["Name", "Length", 	"EffectiveLength", 	"NumReads", "TPM"]]
    temp.columns = col_name 
    temp.to_csv(quant_file, sep ="\t")

rule run_short_reads_genome_alignment:
    input:
        r1 = config["ILLUMINASHORTREADS"]["ill_fastq_R1"],
        r2 = config["ILLUMINASHORTREADS"]["ill_fastq_R2"]

    output:
        rules.all.input[16]
    
    params:
        star_aligner_parameters = "--genomeLoad NoSharedMemory \
                                        --outSAMtype BAM SortedByCoordinate \
                                        --limitBAMsortRAM 100000000000 \
                                        --outSAMunmapped Within KeepPairs \
                                        --twopassMode Basic \
                                        --readFilesCommand zcat",
        star_index = config["GENOMEINDEX"]["star_index"],
        output_prefix = "results/ShortReadSupport/{0}".format(CASENAME)

    log:
        "log/run_short_reads_genome_alignment.log"

    threads:
        18
    
    shell:
        """
            {config[MAPPERS][short_reads_mapper]} {params.star_aligner_parameters} \
                --genomeDir {params.star_index} \
                --runThreadN {threads} \
                --outFileNamePrefix {params.output_prefix} \
                --readFilesIn {input.r1} {input.r2} > {log} 2>&1
            
            touch {output}
        """

rule run_transcriptome_indexing:
    input:
        rules.all.input[14]
    output:
        rules.all.input[17]
    
    params:
        indexoutput = "results/ShortReadSupport/SalmonIndex",
        kmer = 31

    log:
        "log/run_transcriptome_indexing.log"

    threads:
        18
    
    shell:
        """
            salmon index \
            -t {input} \
            -i {params.indexoutput}  \
            -k {params.kmer} \
            -p {threads} > {log} 2>&1
            
            touch {output}
        """


rule run_transcriptome_alignment:
    input:
        index = rules.run_transcriptome_indexing.output[0],
        r1 = rules.run_short_reads_genome_alignment.input.r1,
        r2 = rules.run_short_reads_genome_alignment.input.r2

    output:
        rules.all.input[18]
    
    params:
        output_dir = "results/ShortReadSupport/cupcake.collapsed.filtered.rep.renamed.fasta.salmon",
        indexoutput = "results/ShortReadSupport/SalmonIndex",
        numBootstraps = 50,
        quant_file = "results/ShortReadSupport/cupcake.collapsed.filtered.rep.renamed.fasta.salmon/quant.sf"

    log:
        "log/run_transcriptome_indexing.log"

    threads:
        18
    
    run:
        shell(\
            """
                salmon --no-version-check quant \
                        --validateMappings \
                        -i {params.indexoutput} \
                        --numBootstraps {params.numBootstraps}  \
                        -p {threads} \
                        -l A  \
                        -1 {input.r1} \
                        -2 {input.r2}  \
                        -o  {params.output_dir} > {log} 2>&1
                
                touch {output}
            """)
        
        clean_salmon_output(params.quant_file)



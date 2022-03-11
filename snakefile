import os
import ntpath
import pandas as pd
import tempfile
from Bio import SeqIO
import gzip
import pprint
from collections import Counter
import itertools
import logging


configfile: "config/case.yml"

smrtlink_version = config["SMRTLINKFILES"]["version"]
cutoff = config["FILTERBYCOUNTS"]["default"].split()[1]
CASENAME = config["CASENAME"]

if smrtlink_version < 10:
    print(" Smrtlink version 10 or above needed! ")
    exit

rule all:
    input:
        "results/smrt_link_data/cluster_report.csv",
        "results/smrt_link_data/hq_isoforms.fasta",
        "results/smrt_link_data/flnc.bam",
        "results/Minimap2_Mapping/hq_isoforms_mapped_hg38.sam",
        "results/Minimap2_Mapping/hq_isoforms_mapped_hg38.sorted.sam",
        "results/Collapsed_Isoforms/cupcake.collapsed.gff",
        "results/Collapsed_Isoforms/cupcake.collapsed.rep.fa",
        "results/Collapsed_Isoforms/cupcake.collapsed.group.txt",
        "results/Collapsed_Isoforms/cupcake.collapsed.abundance.txt",
        "results/Collapsed_Isoforms/cupcake.collapsed.min_fl_{filter_count_cutoff}.gff".format(filter_count_cutoff = cutoff),
        "results/Collapsed_Isoforms/cupcake.collapsed.min_fl_{filter_count_cutoff}.abundance.txt".format(filter_count_cutoff = cutoff),
        "results/Collapsed_Isoforms/cupcake.collapsed.min_fl_{filter_count_cutoff}.rep.fa".format(filter_count_cutoff = cutoff),
        "results/Collapsed_Isoforms/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.gff".format(filter_count_cutoff = cutoff),
        "results/Collapsed_Isoforms/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.rep.fa".format(filter_count_cutoff = cutoff),
        "results/Collapsed_Isoforms/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.rep.renamed.fasta".format(filter_count_cutoff = cutoff),
        "results/Collapsed_Isoforms/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.abundance.txt".format(filter_count_cutoff = cutoff),
        "results/ShortReadSupport/short_reads_genome_alignment.status",
        "results/ShortReadSupport/transcriptome_alignment_indexing.status",
        "results/ShortReadSupport/transcriptome_alignment.status",
        "results/IsoformClassification/sqanti3_classification.txt",
        "results/IsoformClassification/sqanti3_junctions.txt",
        "results/IsoformClassification/sqanti3_corrected.fasta",
        "results/IsoformClassification/sqanti3_corrected.gtf",
        "results/IsoformClassificationFiltered/sqanti3_classification.filtered_lite_reasons.txt",
        "results/IsoformClassificationFiltered/sqanti3_classification.filtered_lite_SQANTI3_report.html",
        "results/IsoformClassificationFiltered/sqanti3_classification.filtered_lite.gtf", 
        "results/IsoformClassificationFiltered/sqanti3_classification.filtered_lite_junctions.txt", 
        "results/IsoformClassificationFiltered/sqanti3_classification.filtered_lite_classification.txt",
        "results/IsoformClassificationFiltered/sqanti3_classification.filtered_lite.fasta",
        "results/FinalResults/{0}_squanti3_classification_final_results.tsv".format(CASENAME),
        "results/FinalResults/{0}_squanti3_classification_filtered_out_isoforms_results.tsv".format(CASENAME),
        "results/Collapsed_Isoforms/cupcake.collapsed.read_stat.txt",
        "results/FinalResults/{0}_hq_isoforms_mapped_hg38_sorted.bam".format(CASENAME),
        "results/FinalResults/{0}_sqanti3_classification.filtered_lite_mapped_hg38_sorted.bam".format(CASENAME),
        "results/FinalResults/{0}_ccs_fl.fastq".format(CASENAME),
        "results/FinalResults/{0}_ccs_fl.bam".format(CASENAME),
        "results/FinalResults/{0}_gene_dis_ass_results.tsv".format(CASENAME),
        "results/FinalResults/{0}_pb_isoform_filter_status.tsv".format(CASENAME),
        "results/SVCaller_Output/{0}_hq_isoforms.svsig.gz".format(CASENAME),
        "results/SVCaller_Output/{0}_hq_isoforms.var.vcf".format(CASENAME),
        "results/SVCaller_Output/{0}_hq_isoforms.annot.var.vcf".format(CASENAME),
        "results/SVCaller_Output/{0}.html".format(CASENAME),
        "results/SVCaller_Output/{0}_hq_isoforms.annot.var.tsv".format(CASENAME),
        "results/FinalResults/{0}_toi_list.tsv".format(CASENAME),
        "results/FinalResults/{0}_pbsv_results.tsv".format(CASENAME),
        "results/FinalResults/{0}_squanti3_classification_summary.tsv".format(CASENAME),
        "results/FinalResults/{0}_pbsv_annotation_summary.tsv".format(CASENAME),
        "results/FusionProducts/lq_isoforms.fasta.fusion.gff",
        "results/FusionProducts/lq_isoforms.fasta.fusion.rep.fa",
        "results/FusionProducts/lq_isoforms.fasta.fusion.abundance.txt",
        "results/FusionClassification/lq_isoforms.fasta.fusion_classification.txt",
        "results/FusionClassification/lq_isoforms.fasta.fusion_junctions.txt",
        "results/FusionProducts/lq_isoforms.fasta.fusion.annotated.txt",
        "results/FinalResults/{0}_fusion_classification_final_results.tsv".format(CASENAME),
        "results/FinalResults/{0}_fusion_classification_final_results_fusionhub.tsv".format(CASENAME),
        "results/FinalResults/final_pbfusion.fasta",
        "results/FinalResults/{0}_final_pbfusion_mapped_hg38_sorted.bam".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.gff",
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.rep.fa",
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.group.txt",
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.abundance.txt",
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.ignored_ids.txt",
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.gff".format(filter_count_cutoff = cutoff),
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.rep.fa".format(filter_count_cutoff = cutoff),
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.abundance.txt".format(filter_count_cutoff = cutoff),
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.gff".format(filter_count_cutoff = cutoff),
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.rep.fa".format(filter_count_cutoff = cutoff),
        "results/Deliverables/Isoform/Isoform_Step1/cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered.abundance.txt".format(filter_count_cutoff = cutoff),
        "results/FinalResults/{0}_cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered_mapped_hg38_sorted.bam".format(CASENAME, filter_count_cutoff = cutoff),
        #"results/Deliverables/Isoform/Isoform_Step1/{0}_cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered_mapped_hg38_sorted.bam".format(CASENAME,filter_count_cutoff = cutoff)
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.txt",
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_junctions.txt",
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_corrected.fasta",
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_corrected.gtf",
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite_classification.txt",
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite_junctions.txt",
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite.fasta",
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite.gtf",
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite_reasons.txt",
        "results/Deliverables/Isoform/Isoform_Step2/sqanti3_classification.filtered_lite_SQANTI3_report.html",
        "results/Deliverables/Isoform/Isoform_Step3/{0}_hq_isoforms.svsig.gz".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step3/{0}_hq_isoforms.var.vcf".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step3/{0}_hq_isoforms.annot.var.vcf".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step4/{0}_toi_list.tsv".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step4/{0}_pbsv_results.tsv".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step4/{0}_sqanti3_canonical_results.tsv".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step4/{0}_sqanti3_noncanonical_results.tsv".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step5/gene_dis_ass_results.tsv",
        "results/Deliverables/Fusion/Fusion_Step1/lq_isoforms.fasta.fusion.gff",
        "results/Deliverables/Fusion/Fusion_Step1/lq_isoforms.fasta.fusion.rep.fa",
        "results/Deliverables/Fusion/Fusion_Step1/lq_isoforms.fasta.fusion.abundance.txt",
        "results/Deliverables/Fusion/Fusion_Step2/{0}_hq_isoforms_mapped_hg38_sorted.bam".format(CASENAME),
        "results/Deliverables/Fusion/Fusion_Step2/{0}_hq_isoforms_mapped_hg38_sorted.bam.bai".format(CASENAME),
        "results/Deliverables/Fusion/Fusion_Step2/lq_isoforms.fasta.fusion_classification.txt",
        "results/Deliverables/Fusion/Fusion_Step2/lq_isoforms.fasta.fusion_corrected.fasta",
        "results/Deliverables/Fusion/Fusion_Step2/lq_isoforms.fasta.fusion_corrected.gtf",
        "results/Deliverables/Fusion/Fusion_Step2/refAnnotation_lq_isoforms.fasta.fusion.genePred",
        "results/Deliverables/Fusion/Fusion_Step3/{0}_final_pbfusion_mapped_hg38_sorted.bam".format(CASENAME),
        "results/Deliverables/Fusion/Fusion_Step3/{0}_final_pbfusion_mapped_hg38_sorted.bam.bai".format(CASENAME),
        "results/Deliverables/Fusion/Fusion_Step3/lq_isoforms.fasta.fusion.annotated.txt",
        "results/Deliverables/Fusion/Fusion_Step3/lq_isoforms.fasta.fusion.annotated_ignored.txt",
        "results/Deliverables/Fusion/Fusion_Step4/{0}_fusion_classification_final_results_fusionhub.tsv".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step1/{0}_cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered_mapped_hg38_sorted.bam".format(CASENAME, filter_count_cutoff = cutoff),
        "results/Deliverables/Isoform/Isoform_Step2/{0}_sqanti3_classification.filtered_lite_mapped_hg38_sorted.bam".format(CASENAME),
        "results/FinalResults/{0}_cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered_mapped_hg38_sorted.bam.bai".format(CASENAME,filter_count_cutoff = cutoff),
        "results/Deliverables/Isoform/Isoform_Step1/{0}_cupcake.collapsed.min_fl_{filter_count_cutoff}.filtered_mapped_hg38_sorted.bam.bai".format(CASENAME,filter_count_cutoff = cutoff),
        "results/FinalResults/{0}_sqanti3_classification.filtered_lite_mapped_hg38_sorted.bam.bai".format(CASENAME),
        "results/Deliverables/Isoform/Isoform_Step2/{0}_sqanti3_classification.filtered_lite_mapped_hg38_sorted.bam.bai".format(CASENAME)


include: "workflow/rules/get_smrtlink_files.smk"
include: "workflow/rules/mapping_isoseq.smk"
include: "workflow/rules/picard_sam_sort.smk"
include: "workflow/rules/collapse_isoforms.smk"
include: "workflow/rules/get_abundance.smk"
include: "workflow/rules/filter_by_counts.smk"
include: "workflow/rules/filter_away_subset.smk"
if  config["ILLUMINASHORTREADS"]["ill_fastq_R1"]:
    include: "workflow/rules/short_reads_alignment.smk"
else:
    include: "workflow/rules/short_reads_not_provided.smk"
include: "workflow/rules/sqanti3_classification.smk"
include: "workflow/rules/sqanti3_filter.smk"
include: "workflow/rules/sqanti3_classification_results.smk"
include: "workflow/rules/sam_to_bam.smk"
include: "workflow/rules/mapping_classified_isoforms.smk"
include: "workflow/rules/mapping_ccs_fastq.smk"
include: "workflow/rules/dis_gene_network.smk"
include: "workflow/rules/pbflip_progress.smk"
include: "workflow/rules/sv_caller.smk"
include: "workflow/rules/toi_list.smk"
include: "workflow/rules/summarize_results.smk"
include: "workflow/rules/fusion_finder.smk"
include: "workflow/rules/sqanti3_fusion_classification.smk"
include: "workflow/rules/fusion_collate_info.smk"
include: "workflow/rules/fusion_classification_results.smk"
include: "workflow/rules/fusion_events.smk"
include: "workflow/rules/filter_pbfusion_to_bam.smk"
include: "workflow/rules/deliverables.smk"
include: "workflow/rules/mapping_collapsed_filtered_iso.smk"

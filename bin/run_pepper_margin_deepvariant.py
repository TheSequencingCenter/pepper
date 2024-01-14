import argparse
import sys
from datetime import datetime
import subprocess
import time

__version__ = "0.8.0"


class ModelPaths(object):
    # PEPPER SNP
    PEPPER_MODEL_VARIANT_CALL_ONT_R9_GUPPY5_SUP = "/opt/pepper_models/PEPPER_VARIANT_ONT_R941_GUPPY5_SUP_V8.pkl"
    PEPPER_MODEL_VARIANT_CALL_ONT_R9_GUPPY4_HAC = "/opt/pepper_models/PEPPER_VARIANT_ONT_R941_GUPPY4_HAC_V8.pkl"
    PEPPER_MODEL_VARIANT_CALL_ONT_R10_Q20 = "/opt/pepper_models/PEPPER_VARIANT_ONT_R10_Q20_V8.pkl"
    PEPPER_MODEL_VARIANT_CALL_HIFI = "/opt/pepper_models/PEPPER_VARIANT_HIFI_V8.pkl"

    # DEEPVARIANT
    DEEPVARIANT_MODEL_VARIANT_CALL_ONT_R9_GUPPY5_SUP_SNP = "/opt/dv_models/ont_deepvariant_vc/dv_ont_r9_guppy5_sup_vc_model_snp"
    DEEPVARIANT_MODEL_VARIANT_CALL_ONT_R9_GUPPY5_SUP_INDEL = "/opt/dv_models/ont_deepvariant_vc/dv_ont_r9_guppy5_sup_vc_model_indel"
    DEEPVARIANT_MODEL_VARIANT_CALL_ONT_R10_q20_SNP = "/opt/dv_models/ont_deepvariant_vc/dv_ont_r10_q20_vc_model_snp"
    DEEPVARIANT_MODEL_VARIANT_CALL_ONT_R10_q20_INDEL = "/opt/dv_models/ont_deepvariant_vc/dv_ont_r10_q20_vc_model_indel"
    DEEPVARIANT_MODEL_HIFI = "/opt/dv_models/ont_deepvariant_vc/dv_hifi_model"

    # MARGIN
    MARGIN_HAPLOTAG_ONT = "/opt/margin_dir/params/phase/allParams.haplotag.ont-r94g507.snp.json"
    MARGIN_HAPLOTAG_HIFI = "/opt/margin_dir/params/phase/allParams.haplotag.pb-hifi.snp.json"
    MARGIN_PHASE_MODEL_ONT = "/opt/margin_dir/params/phase/allParams.phase_vcf.ont.json"
    MARGIN_PHASE_MODEL_HIFI = "/opt/margin_dir/params/phase/allParams.phase_vcf.pb-hifi.json"


# STEP 1: PEPPER SNP
def get_pepper_command(options):
    options.output_dir_intermediate = options.output_dir + "/intermediate_files"
    preprocess_command = ""
    preprocess_command = preprocess_command + "mkdir -p " + options.output_dir + "; \n"
    preprocess_command = preprocess_command + "mkdir -p " + options.output_dir + "/logs; \n"
    preprocess_command = preprocess_command + "mkdir -p " + options.output_dir_intermediate + "; \n"
    preprocess_command = preprocess_command + "cp " + options.pepper_model + " " + options.output_dir_intermediate
    options.pepper_model = options.output_dir_intermediate + "/" + options.pepper_model.rstrip().split('/')[-1]

    # build the pepper command
    pepper_command = "time pepper_variant call_variant "
    pepper_command += "-b " + options.bam + " "
    pepper_command += "-f " + options.fasta + " "
    pepper_command += "-t " + str(options.threads) + " "
    pepper_command += "-m " + str(options.pepper_model) + " "
    pepper_command += "-o " + str(options.output_dir) + "/pepper/ "

    if options.pepper_downsample_rate is not None:
        pepper_command += "--downsample_rate " + str(options.pepper_downsample_rate) + " "
    if options.pepper_region_size is not None:
        pepper_command += "--region_size " + str(options.pepper_region_size) + " "
    if options.pepper_include_supplementary:
        pepper_command += "--include_supplementary " + " "
    if options.pepper_min_mapq is not None:
        pepper_command += "--min_mapq " + str(options.pepper_min_mapq) + " "
    if options.pepper_min_snp_baseq is not None:
        pepper_command += "--min_snp_baseq " + str(options.pepper_min_snp_baseq) + " "
    if options.pepper_min_indel_baseq is not None:
        pepper_command += "--min_indel_baseq " + str(options.pepper_min_indel_baseq) + " "
    if options.pepper_snp_frequency is not None:
        pepper_command += "--snp_frequency " + str(options.pepper_snp_frequency) + " "
    if options.pepper_insert_frequency is not None:
        pepper_command += "--insert_frequency " + str(options.pepper_insert_frequency) + " "
    if options.pepper_delete_frequency is not None:
        pepper_command += "--delete_frequency " + str(options.pepper_delete_frequency) + " "
    if options.pepper_min_coverage_threshold is not None:
        pepper_command += "--min_coverage_threshold " + str(options.pepper_min_coverage_threshold) + " "
    if options.pepper_candidate_support_threshold is not None:
        pepper_command += "--candidate_support_threshold " + str(options.pepper_candidate_support_threshold) + " "
    if options.pepper_snp_candidate_frequency_threshold is not None:
        pepper_command += "--snp_candidate_frequency_threshold " + str(options.pepper_snp_candidate_frequency_threshold) + " "
    if options.pepper_indel_candidate_frequency_threshold is not None:
        pepper_command += "--indel_candidate_frequency_threshold " + str(options.pepper_indel_candidate_frequency_threshold) + " "
    if options.pepper_skip_indels:
        pepper_command += "--skip_indels " + " "
    if options.pepper_allowed_multiallelics is not None:
        pepper_command += "--allowed_multiallelics " + str(options.pepper_allowed_multiallelics) + " "
    if options.pepper_snp_p_value is not None:
        pepper_command += "--snp_p_value " + str(options.pepper_snp_p_value) + " "
    if options.pepper_insert_p_value is not None:
        pepper_command += "--insert_p_value " + str(options.pepper_insert_p_value) + " "
    if options.pepper_delete_p_value is not None:
        pepper_command += "--delete_p_value " + str(options.pepper_delete_p_value) + " "

    if options.pepper_snp_p_value_in_lc is not None:
        pepper_command += "--snp_p_value_in_lc " + str(options.pepper_snp_p_value_in_lc) + " "
    if options.pepper_insert_p_value_in_lc is not None:
        pepper_command += "--insert_p_value_in_lc " + str(options.pepper_insert_p_value_in_lc) + " "
    if options.pepper_delete_p_value_in_lc is not None:
        pepper_command += "--delete_p_value_in_lc " + str(options.pepper_delete_p_value_in_lc) + " "

    if options.pepper_snp_q_cutoff is not None:
        pepper_command += "--snp_q_cutoff " + str(options.pepper_snp_q_cutoff) + " "
    if options.pepper_indel_q_cutoff is not None:
        pepper_command += "--indel_q_cutoff " + str(options.pepper_indel_q_cutoff) + " "

    if options.pepper_snp_q_cutoff_in_lc is not None:
        pepper_command += "--snp_q_cutoff_in_lc " + str(options.pepper_snp_q_cutoff_in_lc) + " "
    if options.pepper_indel_q_cutoff_in_lc is not None:
        pepper_command += "--indel_q_cutoff_in_lc " + str(options.pepper_indel_q_cutoff_in_lc) + " "

    if options.pepper_report_snp_above_freq is not None:
        pepper_command += "--report_snp_above_freq " + str(options.pepper_report_snp_above_freq) + " "
    if options.pepper_report_indel_above_freq is not None:
        pepper_command += "--report_indel_above_freq " + str(options.pepper_report_indel_above_freq) + " "
    if options.pepper_quantized is False:
        pepper_command += "--no_quantized " + " "
    if options.pepper_quantized is True:
        pepper_command += "--quantized " + " "

    if options.region is not None:
        pepper_command += "-r " + str(options.region) + " "
    if options.sample_name is not None:
        pepper_command += "-s " + str(options.sample_name) + " "

    if options.gpu:
        pepper_command += "-g "

    if options.ont_r9_guppy5_sup:
        pepper_command += "--ont_r9_guppy5_sup"
    elif options.ont_r10_q20:
        pepper_command += "--ont_r10_q20"
    elif options.hifi:
        pepper_command += "--hifi"
    elif options.clr:
        pepper_command += "--clr"

    pepper_command += " 2>&1 | tee " + options.output_dir + "/logs/1_pepper.log"

    post_process_command = ""
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_FULL.vcf.gz " + options.output_dir_intermediate + "/; \n"
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_FULL.vcf.gz.tbi " + options.output_dir_intermediate + "/; \n"
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz " + options.output_dir_intermediate + "/; \n"
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz.tbi " + options.output_dir_intermediate + "/; \n"
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz " + options.output_dir_intermediate + "/; \n"
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz.tbi " + options.output_dir_intermediate + "/; \n"
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs.vcf.gz " + options.output_dir_intermediate + "/; \n"
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs.vcf.gz.tbi " + options.output_dir_intermediate + "/; \n"
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_INDEL.vcf.gz " + options.output_dir_intermediate + "/; \n"
    post_process_command = post_process_command + "mv " + options.output_dir + "/pepper/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_INDEL.vcf.gz.tbi " + options.output_dir_intermediate + "/; \n"

    # remove pepper directory
    post_process_command = post_process_command + "rm -rf " + options.output_dir + "/pepper/; \n"

    pepper_output_vcf_variant_calling = options.output_dir_intermediate + "/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz"
    pepper_output_vcf_variant_calling_snp = options.output_dir_intermediate + "/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs.vcf.gz"
    pepper_output_vcf_variant_calling_indel = options.output_dir_intermediate + "/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_INDEL.vcf.gz"
    pepper_output_full = options.output_dir_intermediate + "/PEPPER_VARIANT_FULL.vcf.gz"
    pepper_output_pepper = options.output_dir_intermediate + "/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz"

    # print contigs found
    post_process_command = post_process_command + "echo \"CONTIGS FOUND IN PEPPER VCF:\"; \n"
    post_process_command = post_process_command + "zcat " + pepper_output_full + " | grep -v '#' | cut -f1 | uniq"

    return preprocess_command, pepper_command, post_process_command, pepper_output_vcf_variant_calling, pepper_output_vcf_variant_calling_snp, pepper_output_vcf_variant_calling_indel, pepper_output_full, pepper_output_pepper


# STEP 2: MARGIN HAP
def get_margin_haplotag_command(options):
    margin_command = "time margin phase "
    margin_command = margin_command + options.bam + " "
    margin_command = margin_command + options.fasta + " "
    margin_command = margin_command + options.pepper_full_vcf + " "
    margin_command = margin_command + options.margin_haplotag_model + " "
    margin_command = margin_command + "-t " + str(options.threads) + " "
    if options.region:
        margin_command += "-r " + options.region + " "

    margin_command = margin_command + "-V "
    margin_command = margin_command + "-o " + options.output_dir_intermediate + "/PHASED.PEPPER_MARGIN 2>&1 | tee " + options.output_dir + "/logs/2_margin_haplotag.log;\n"

    margin_command = margin_command + "samtools index -@" + str(options.threads) + " " + options.output_dir_intermediate + "/PHASED.PEPPER_MARGIN.haplotagged.bam"

    margin_output_bam = options.output_dir_intermediate + "/PHASED.PEPPER_MARGIN.haplotagged.bam"

    return margin_command, margin_output_bam


def get_make_example_args(options):
    extra_args = []
    if options.dv_alt_aligned_pileup:
        extra_args.append("alt_aligned_pileup=" + str(options.dv_alt_aligned_pileup))
    if options.dv_realign_reads:
        extra_args.append("realign_reads=" + str(options.dv_realign_reads))
    if options.dv_partition_size:
        extra_args.append("partition_size=" + str(options.dv_partition_size))
    if options.dv_min_mapping_quality:
        extra_args.append("min_mapping_quality=" + str(options.dv_min_mapping_quality))
    if options.dv_min_base_quality:
        extra_args.append("min_base_quality=" + str(options.dv_min_base_quality))
    if options.dv_vsc_min_fraction_snps:
        extra_args.append("vsc_min_fraction_snps=" + str(options.dv_vsc_min_fraction_snps))
    if options.dv_vsc_min_fraction_indels:
        extra_args.append("vsc_min_fraction_indels=" + str(options.dv_vsc_min_fraction_indels))
    if options.dv_pileup_image_width:
        extra_args.append("pileup_image_width=" + str(options.dv_pileup_image_width))
    if options.dv_sort_by_haplotypes:
        extra_args.append("sort_by_haplotypes=" + str(options.dv_sort_by_haplotypes))
    if options.dv_parse_sam_aux_fields:
        extra_args.append("parse_sam_aux_fields=" + str(options.dv_parse_sam_aux_fields))
    if options.dv_add_hp_channel:
        extra_args.append("add_hp_channel=" + str(options.dv_add_hp_channel))
    if options.variant_caller:
        extra_args.append("variant_caller=" + str(options.variant_caller))
    if options.proposed_variants:
        extra_args.append("proposed_variants=" + str(options.proposed_variants))

    return ','.join(extra_args)


# STEP 3: DeepVariant
def get_deepvariant_command(options, deepvariant_output_filename, log_name):
    preprocess_command = "mkdir -p " + options.output_dir + "/dv_intermediate_outputs; \n"
    preprocess_command = preprocess_command + "echo \"STARTING DEEPVARIANT\"; \n"

    make_examples_extra_args = get_make_example_args(options)

    # build deepvariant command
    deepvariant_command = "time /opt/deepvariant/bin/run_deepvariant "
    if options.dv_model and not options.hifi:
        deepvariant_command = deepvariant_command + "--model_type=WGS --customized_model=" + options.dv_model + " "
    elif options.hifi and options.dv_model:
        deepvariant_command = deepvariant_command + "--model_type=PACBIO --customized_model=" + options.dv_model + " "
    else:
        raise ValueError('No DeepVariant model provided')

    deepvariant_command = deepvariant_command + "--ref=" + options.fasta + " "
    deepvariant_command = deepvariant_command + "--reads=" + options.margin_output + " "
    deepvariant_command = deepvariant_command + "--output_vcf=" + options.output_dir_intermediate + "/" + deepvariant_output_filename + ".vcf.gz "

    if options.gvcf:
        deepvariant_command = deepvariant_command + "--output_gvcf=" + options.output_dir + "/" + deepvariant_output_filename + ".g.vcf.gz "
    if options.region:
        deepvariant_command = deepvariant_command + "--regions=\"" + options.region + "\" "
    if options.sample_name:
        deepvariant_command = deepvariant_command + "--sample_name=\"" + options.sample_name + "\" "
    if options.dv_use_hp_information:
        deepvariant_command = deepvariant_command + "--use_hp_information "

    deepvariant_command = deepvariant_command + "--intermediate_results_dir=" + options.output_dir + "/dv_intermediate_outputs/ "
    deepvariant_command = deepvariant_command + "--num_shards=" + str(options.threads) + " "

    deepvariant_command = deepvariant_command + "--make_examples_extra_args=\"" + make_examples_extra_args + "\" "

    if options.dv_use_multiallelic_mode:
        deepvariant_command = deepvariant_command + "--postprocess_variants_extra_args=\"use_multiallelic_model=True\" 2>&1 | tee " + options.output_dir + "/logs/" + str(log_name) + ".log"
    else:
        deepvariant_command = deepvariant_command + " 2>&1 | tee " + options.output_dir + "/logs/4_DeepVariant.log"

    deepvariant_output = options.output_dir_intermediate + "/" + deepvariant_output_filename + ".vcf.gz"

    return deepvariant_command, deepvariant_output


# STEP 4: Merge VCF command
def get_merge_vcf_command(options, snp_indel_separate):
    pepper_command = "time pepper_variant merge_variants "
    pepper_command += "-vp " + options.pepper_full_vcf + " "
    if not snp_indel_separate:
        pepper_command += "-vd " + options.deepvariant_output + " "
    else:
        pepper_command += "-vds " + options.deepvariant_output_snp + " "
        pepper_command += "-vdi " + options.deepvariant_output_indel + " "
    pepper_command += "-o " + str(options.output_dir)
    pepper_command += " 2>&1 | tee " + options.output_dir + "/logs/6_merge_vcf.log"

    output_vcf = str(options.output_dir) + "/" + options.final_vcf_name + '.vcf.gz'
    # index the full vcf
    post_process_command = "mv " + str(options.output_dir) + "/" + "PEPPER_MARGIN_DEEPVARIANT_OUTPUT.vcf.gz" + " " + str(output_vcf) + ";\n"
    post_process_command = post_process_command + "mv " + str(options.output_dir) + "/" + "PEPPER_MARGIN_DEEPVARIANT_OUTPUT.vcf.gz.tbi" + " " + str(output_vcf) + ".tbi;\n"
    # generate report
    post_process_command = post_process_command + "/opt/deepvariant/bin/vcf_stats_report --input_vcf " + output_vcf + " --outfile_base " + str(options.output_dir) + "/" + options.final_vcf_name + "; \n"

    return pepper_command, post_process_command, output_vcf


# STEP 5: Margin Phase
def get_margin_phase_command(options):
    margin_command = "time margin phase "
    margin_command = margin_command + options.bam + " "
    margin_command = margin_command + options.fasta + " "
    margin_command = margin_command + options.final_output_vcf + " "
    margin_command = margin_command + options.margin_phase_model + " "
    margin_command = margin_command + "-t " + str(options.threads) + " "
    if options.region:
        margin_command += "-r " + options.region + " "

    if options.skip_final_phased_bam:
        margin_command = margin_command + "-M "
    margin_command = margin_command + "-o " + options.output_dir + "/" + options.final_output_vcf_phased + " 2>&1 | tee " + options.output_dir + "/logs/7_margin_phase.log; \n"

    margin_command = margin_command + "bgzip " + options.output_dir + "/" + options.final_output_vcf_phased + ".phased.vcf; \n"
    margin_command = margin_command + "tabix -p vcf " + options.output_dir + "/" + options.final_output_vcf_phased + ".phased.vcf.gz"

    return margin_command


# STEP 6: Post-process variant call
def get_post_process_variant_call(options):
    post_processing_command = ""
    post_processing_command = post_processing_command + "rm -rf " + options.output_dir + "/dv_intermediate_outputs;\n"

    if not options.only_haplotag:
        if not options.keep_intermediate_bam_files:
            post_processing_command = post_processing_command + "rm -rf " + options.margin_output + ";\n"
            post_processing_command = post_processing_command + "rm -rf " + options.margin_output + ".bai;\n"

        if options.gvcf:
            post_processing_command = post_processing_command + "mv " + options.output_dir + "/" + options.deepvariant_output_filename + ".g.vcf.gz " + " " + options.output_dir + "/" + options.final_vcf_name + ".g.vcf.gz;\n"
            post_processing_command = post_processing_command + "mv " + options.output_dir + "/" + options.deepvariant_output_filename + ".g.vcf.gz.tbi " + " " + options.output_dir + "/" + options.final_vcf_name + ".g.vcf.gz.tbi;\n"

    return post_processing_command


def run_variant_calling(options):
    start_time = time.time()
    # SETUP MODELS AND PARAMETERS
    if options.ont_r9_guppy5_sup:
        # dv chunk options
        if not options.dv_partition_size:
            options.dv_partition_size = 10000
        # model options
        if not options.pepper_model:
            options.pepper_model = ModelPaths.PEPPER_MODEL_VARIANT_CALL_ONT_R9_GUPPY5_SUP
        if not options.margin_haplotag_model:
            options.margin_haplotag_model = ModelPaths.MARGIN_HAPLOTAG_ONT
        # setup deepvariant models
        if not options.dv_model:
            if not options.dv_model_snp:
                options.dv_model_snp = ModelPaths.DEEPVARIANT_MODEL_VARIANT_CALL_ONT_R9_GUPPY5_SUP_SNP
            if not options.dv_model_indel:
                options.dv_model_indel = ModelPaths.DEEPVARIANT_MODEL_VARIANT_CALL_ONT_R9_GUPPY5_SUP_INDEL

        if not options.margin_phase_model:
            options.margin_phase_model = ModelPaths.MARGIN_PHASE_MODEL_ONT

        # DeepVariant parameters
        if not options.dv_alt_aligned_pileup:
            if not options.dv_alt_aligned_pileup_snp:
                options.dv_alt_aligned_pileup_snp = "none"
            if not options.dv_alt_aligned_pileup_indel:
                options.dv_alt_aligned_pileup_indel = "rows"

        if not options.dv_realign_reads:
            options.dv_realign_reads = "false"
        if options.dv_min_mapping_quality is None:
            options.dv_min_mapping_quality = 5
        if options.dv_min_base_quality is None:
            options.dv_min_base_quality = 1
        if options.dv_sort_by_haplotypes is None:
            options.dv_sort_by_haplotypes = "true"
        if options.dv_parse_sam_aux_fields is None:
            options.dv_parse_sam_aux_fields = "true"
        if options.dv_add_hp_channel is None:
            options.dv_add_hp_channel = "true"
        if options.dv_use_multiallelic_mode is None:
            options.dv_use_multiallelic_mode = "true"
        options.variant_caller = "vcf_candidate_importer"
    elif options.ont_r10_q20:
        # dv chunk options
        if not options.dv_partition_size:
            options.dv_partition_size = 10000

        if not options.pepper_model:
            options.pepper_model = ModelPaths.PEPPER_MODEL_VARIANT_CALL_ONT_R10_Q20
        if not options.margin_haplotag_model:
            options.margin_haplotag_model = ModelPaths.MARGIN_HAPLOTAG_ONT
        # setup deepvariant models
        if not options.dv_model:
            if not options.dv_model_snp:
                options.dv_model_snp = ModelPaths.DEEPVARIANT_MODEL_VARIANT_CALL_ONT_R10_q20_SNP
            if not options.dv_model_indel:
                options.dv_model_indel = ModelPaths.DEEPVARIANT_MODEL_VARIANT_CALL_ONT_R10_q20_INDEL

        if not options.margin_phase_model:
            options.margin_phase_model = ModelPaths.MARGIN_PHASE_MODEL_ONT

        # DeepVariant parameters
        if not options.dv_alt_aligned_pileup:
            if not options.dv_alt_aligned_pileup_snp:
                options.dv_alt_aligned_pileup_snp = "none"
            if not options.dv_alt_aligned_pileup_indel:
                options.dv_alt_aligned_pileup_indel = "rows"

        if not options.dv_realign_reads:
            options.dv_realign_reads = "false"
        if options.dv_min_mapping_quality is None:
            options.dv_min_mapping_quality = 1
        if options.dv_min_base_quality is None:
            options.dv_min_base_quality = 1
        if options.dv_sort_by_haplotypes is None:
            options.dv_sort_by_haplotypes = "true"
        if options.dv_parse_sam_aux_fields is None:
            options.dv_parse_sam_aux_fields = "true"
        if options.dv_add_hp_channel is None:
            options.dv_add_hp_channel = "true"
        if options.dv_use_multiallelic_mode is None:
            options.dv_use_multiallelic_mode = "true"
        options.variant_caller = "vcf_candidate_importer"
    elif options.hifi:
        # dv chunk options
        if not options.dv_partition_size:
            options.dv_partition_size = 1000

        # pipeline options
        if not options.pepper_model:
            options.pepper_model = ModelPaths.PEPPER_MODEL_VARIANT_CALL_HIFI
        if not options.margin_haplotag_model:
            options.margin_haplotag_model = ModelPaths.MARGIN_HAPLOTAG_HIFI
        if not options.dv_model:
            options.dv_model = ModelPaths.DEEPVARIANT_MODEL_HIFI
        if not options.margin_phase_model:
            options.margin_phase_model = ModelPaths.MARGIN_PHASE_MODEL_HIFI

        if options.dv_model:
            if options.dv_min_mapping_quality is None:
                options.dv_min_mapping_quality = 5
            if options.dv_min_base_quality is None:
                options.dv_min_base_quality = 10
            if options.dv_vsc_min_fraction_snps is None:
                options.dv_vsc_min_fraction_snps = 0.10
            if options.dv_vsc_min_fraction_indels is None:
                options.dv_vsc_min_fraction_indels = 0.10
            if options.dv_pileup_image_width is None:
                options.dv_pileup_image_width = 221
            if options.dv_use_multiallelic_mode is None:
                options.dv_use_multiallelic_mode = "true"
            options.variant_caller = "vcf_candidate_importer"
            options.dv_use_hp_information = True
        else:
            options.variant_caller = "vcf_candidate_importer"
            options.dv_use_hp_information = True
    else:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: INVALID OPTION\n")
        exit()

    pepper_preprocess, pepper_command, pepper_postprocess, pepper_out_variant_calling, pepper_out_variant_calling_snp, pepper_out_variant_calling_indel, pepper_output_full, pepper_output_pepper = get_pepper_command(options)

    options.pepper_variant_calling = pepper_out_variant_calling
    options.pepper_variant_calling_snp = pepper_out_variant_calling_snp
    options.pepper_variant_calling_indel = pepper_out_variant_calling_indel
    options.pepper_full_vcf = pepper_output_full
    options.pepper_vcf_only_pepper = pepper_output_pepper

    commands_list = list()
    # STEP 1: RUN PEPPER SNP
    commands_list.append(pepper_preprocess)
    commands_list.append(pepper_command)
    commands_list.append(pepper_postprocess)

    options.deepvariant_output = None
    options.deepvariant_output_snp = None
    options.deepvariant_output_indel = None

    if not options.only_pepper:
        # STEP 2: RUN MARGIN HAPLOTAG
        margin_command, margin_output = get_margin_haplotag_command(options)
        commands_list.append(margin_command)
        options.margin_output = margin_output

    if not options.only_pepper and not options.gvcf:
        options.proposed_variants = options.pepper_variant_calling
    else:
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: GVCF OPTION IS ON, FULL RE-GENOTYPING WILL BE PERFORMED.\n")
        options.proposed_variants = options.pepper_full_vcf

    if not options.only_pepper and not options.only_haplotag:
        # STEP 5: RUN DEEPVARIANT
        if options.dv_model:
            deepvariant_command, deepvariant_output = get_deepvariant_command(options, "DEEPVARIANT_OUTPUT", "5_DeepVariant")
            commands_list.append(deepvariant_command)
            options.deepvariant_output_filename = "DEEPVARIANT_OUTPUT"
            options.deepvariant_output = deepvariant_output
            options.final_vcf_name = options.output_prefix

            merge_command, post_process_command, merged_vcf = get_merge_vcf_command(options, snp_indel_separate=False)
            commands_list.append(merge_command)
            commands_list.append(post_process_command)
            options.final_output_vcf = merged_vcf
        elif not options.gvcf:
            # setup SNP calling
            options.proposed_variants = options.pepper_variant_calling_snp
            options.dv_model = options.dv_model_snp
            options.dv_alt_aligned_pileup = options.dv_alt_aligned_pileup_snp

            deepvariant_command_snp, deepvariant_output_snp = get_deepvariant_command(options, "DEEPVARIANT_OUTPUT_SNP", "5.1_DeepVariant_SNP")
            options.deepvariant_output_filename_snp = "DEEPVARIANT_OUTPUT_SNP"
            commands_list.append(deepvariant_command_snp)
            options.deepvariant_output_snp = deepvariant_output_snp

            # setup indel calling
            options.proposed_variants = options.pepper_variant_calling_indel
            options.dv_model = options.dv_model_indel
            options.dv_alt_aligned_pileup = options.dv_alt_aligned_pileup_indel

            deepvariant_command_indel, deepvariant_output_indel = get_deepvariant_command(options, "DEEPVARIANT_OUTPUT_INDEL", "5.2_DeepVariant_INDEL")
            options.deepvariant_output_filename_indel = "DEEPVARIANT_OUTPUT_INDEL"
            commands_list.append(deepvariant_command_indel)
            options.deepvariant_output_indel = deepvariant_output_indel
            options.final_vcf_name = options.output_prefix

            merge_command, post_process_command, merged_vcf = get_merge_vcf_command(options, snp_indel_separate=True)
            commands_list.append(merge_command)
            commands_list.append(post_process_command)
            options.final_output_vcf = merged_vcf
        else:
            # use indel model
            options.proposed_variants = options.pepper_full_vcf
            options.dv_model = options.dv_model_snp
            options.dv_alt_aligned_pileup = options.dv_alt_aligned_pileup_indel

            deepvariant_command, deepvariant_output = get_deepvariant_command(options, "DEEPVARIANT_OUTPUT", "5_DeepVariant")
            options.deepvariant_output_filename = "DEEPVARIANT_OUTPUT"
            commands_list.append(deepvariant_command)
            options.deepvariant_output = deepvariant_output
            options.final_vcf_name = options.output_prefix

            merge_command, post_process_command, merged_vcf = get_merge_vcf_command(options, snp_indel_separate=False)
            commands_list.append(merge_command)
            commands_list.append(post_process_command)
            options.final_output_vcf = merged_vcf

        # STEP 6: RUN MARGIN PHASE
        if options.phased_output:
            options.final_output_vcf_phased = options.final_vcf_name
            margin_phase_command = get_margin_phase_command(options)
            commands_list.append(margin_phase_command)

    post_processing_command = get_post_process_variant_call(options)
    commands_list.append(post_processing_command)

    # execute command
    for i, command in enumerate(commands_list):
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: " + "[" + str(i + 1) + "/" + str(len(commands_list)) + "] " + "RUNNING THE FOLLOWING COMMAND\n-------\n" + command + "\n-------\n")
        sys.stderr.flush()

        if not options.dry:
            try:
                ret_code = subprocess.check_call(command, shell=True, executable='/bin/bash')
            except subprocess.CalledProcessError as e:
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] ERROR: " + str(e.output) + "]\n")
                sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] THE FOLLOWING COMMAND FAILED: " + str(command) + "]\n")
                exit(1)

    end_time = time.time()
    mins = int((end_time - start_time) / 60)
    secs = int((end_time - start_time)) % 60
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: TOTAL ELAPSED TIME FOR VARIANT CALLING: " + str(mins) + " Min " + str(secs) + " Sec\n")


def add_variant_call_arguments(parser):
    """
    Add arguments to a parser for sub-command "variant call"
    :param parser: argparse object
    :return:
    """
    required = parser.add_argument_group('Required Arguments')
    required.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        help="Alignment containing mapping between reads and a reference."
    )
    required.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=True,
        help="A reference file in FASTA format."
    )
    required.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory."
    )
    required.add_argument(
        "-t",
        "--threads",
        required=True,
        type=int,
        help="Number of threads to use."
    )

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument(
        "-r",
        "--region",
        type=str,
        default=None,
        help="Region in [contig_name:start-end] format. Default is None."
    )
    optional.add_argument(
        "-g",
        "--gpu",
        default=False,
        action='store_true',
        required=False,
        help="If set then will use GPUs for inference. CUDA required."
    )
    optional.add_argument(
        "--gvcf",
        default=False,
        dest='gvcf',
        action='store_true',
        required=False,
        help="If set then a gVCF output will be generated."
    )
    optional.add_argument(
        "-s",
        "--sample_name",
        type=str,
        required=False,
        default="Sample",
        help="Name of the sample to put in the output VCF."
    )
    optional.add_argument(
        "--only_pepper",
        default=False,
        action='store_true',
        required=False,
        help="If set then pipeline will stop after PEPPER."
    )
    optional.add_argument(
        "--only_haplotag",
        default=False,
        action='store_true',
        required=False,
        help="If set then pipeline will stop after Margin haplotag stage."
    )
    optional.add_argument(
        "--phased_output",
        default=False,
        action='store_true',
        required=False,
        help="If set then Margin phase will generate a phased VCF and BAM at the end when --phased_output is set."
    )
    optional.add_argument(
        "--skip_final_phased_bam",
        default=False,
        action='store_true',
        required=False,
        help="If true with phased output then the final output will not have a bam file when --phased_output is set."
    )
    optional.add_argument(
        "-p",
        "--output_prefix",
        type=str,
        required=False,
        default="PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT",
        help="Prefix for output filename. Do not include extension in prefix."
    )
    optional.add_argument(
        "-k",
        "--keep_intermediate_bam_files",
        default=False,
        action='store_true',
        required=False,
        help="If true then intermediate haplotagged bam files will be kept in the intermediate files directory. Default: False"
    )

    optional.add_argument(
        "--dry",
        default=False,
        action='store_true',
        required=False,
        help="If true then only the commands will be printed. [Debugging]"
    )

    pepper = parser.add_argument_group('Arguments for PEPPER')
    pepper.add_argument(
        "--pepper_model",
        type=str,
        help="Use a custom PEPPER model."
    )
    pepper.add_argument(
        "--pepper_quantized",
        default=False,
        action='store_true',
        help="PEPPER: Use quantization for inference while on CPU inference mode. Speeds up inference. Default is True."
    )
    pepper.add_argument(
        "--no_pepper_quantized",
        dest='pepper_quantized',
        default=False,
        action='store_false',
        help="Do not use quantization for inference while on CPU inference mode. Speeds up inference."
    )
    pepper.add_argument(
        "--pepper_downsample_rate",
        type=float,
        default=None,
        help="Downsample rate of reads while generating images. Default is 1.0"
    )
    pepper.add_argument(
        "--pepper_region_size",
        type=int,
        required=False,
        default=None,
        help="Region size in bp used to chunk the genome. Default is 100000."
    )
    pepper.add_argument(
        "--pepper_include_supplementary",
        default=False,
        action='store_true',
        help="If true then supplementary reads will be used. Default is False."
    )
    pepper.add_argument(
        "--pepper_min_mapq",
        type=int,
        required=False,
        default=None,
        help="Minimum mapping quality for read to be considered valid. Default is 5"
    )
    pepper.add_argument(
        "--pepper_min_snp_baseq",
        type=int,
        required=False,
        default=None,
        help="Minimum base quality for base to be considered valid for SNP. Default is 1"
    )
    pepper.add_argument(
        "--pepper_min_indel_baseq",
        type=int,
        required=False,
        default=None,
        help="Minimum base quality for base to be considered valid for INDELs."
    )
    pepper.add_argument(
        "--pepper_snp_frequency",
        type=float,
        required=False,
        default=None,
        help="Minimum SNP frequency for a site to be considered to have a variant."
    )
    pepper.add_argument(
        "--pepper_insert_frequency",
        type=float,
        required=False,
        default=None,
        help="Minimum insert frequency for a site to be considered to have a variant."
    )
    pepper.add_argument(
        "--pepper_delete_frequency",
        type=float,
        required=False,
        default=None,
        help="Minimum delete frequency for a site to be considered to have a variant."
    )
    pepper.add_argument(
        "--pepper_min_coverage_threshold",
        type=int,
        required=False,
        default=None,
        help="Minimum delete frequency for a site to be considered to have a variant."
    )
    pepper.add_argument(
        "--pepper_candidate_support_threshold",
        type=int,
        required=False,
        default=None,
        help="Minimum number of reads supporting a variant to be considered as a candidate."
    )
    pepper.add_argument(
        "--pepper_snp_candidate_frequency_threshold",
        type=float,
        required=False,
        default=None,
        help="Minimum frequency for a SNP candidate to be considered to be a variant."
    )
    pepper.add_argument(
        "--pepper_indel_candidate_frequency_threshold",
        type=float,
        required=False,
        default=None,
        help="Minimum frequency for an INDEL candidate to be considered to be a variant."
    )
    pepper.add_argument(
        "--pepper_skip_indels",
        default=False,
        action='store_true',
        help="If set then INDEL calling will be skipped."
    )
    pepper.add_argument(
        "--pepper_allowed_multiallelics",
        required=False,
        type=int,
        default=None,
        help="Max allowed multiallelic variants per site."
    )
    pepper.add_argument(
        "--pepper_snp_p_value",
        required=False,
        type=float,
        default=None,
        help="Predicted value used for a SNP to be considered a candidate."
    )
    pepper.add_argument(
        "--pepper_insert_p_value",
        required=False,
        type=float,
        default=None,
        help="Predicted value used for a insert to be considered a candidate."
    )
    pepper.add_argument(
        "--pepper_delete_p_value",
        required=False,
        type=float,
        default=None,
        help="Predicted value used for a delete to be considered a candidate."
    )
    pepper.add_argument(
        "--pepper_snp_p_value_in_lc",
        required=False,
        type=float,
        default=None,
        help="Predicted value used for a SNP to be considered a candidate in low complexity regions."
    )
    pepper.add_argument(
        "--pepper_insert_p_value_in_lc",
        required=False,
        type=float,
        default=None,
        help="Predicted value used for an insert to be considered a candidate in low complexity regions."
    )
    pepper.add_argument(
        "--pepper_delete_p_value_in_lc",
        required=False,
        type=float,
        default=None,
        help="Predicted value used for a delete to be considered a candidate in low complexity regions."
    )
    pepper.add_argument(
        "--pepper_snp_q_cutoff",
        required=False,
        type=float,
        default=None,
        help="GQ cutoff for a SNP variant to be re-genotyped with DeepVariant. Variants with GQ below this will be re-genotyped."
    )
    pepper.add_argument(
        "--pepper_indel_q_cutoff",
        required=False,
        type=float,
        default=None,
        help="GQ cutoff for an INDEL variant to be re-genotyped with DeepVariant. Variants with GQ below this will be re-genotyped."
    )
    pepper.add_argument(
        "--pepper_snp_q_cutoff_in_lc",
        required=False,
        type=float,
        default=None,
        help="GQ cutoff for a SNP variant in low complexity region to be re-genotyped with DeepVariant. Variants with GQ below this will be re-genotyped."
    )
    pepper.add_argument(
        "--pepper_indel_q_cutoff_in_lc",
        required=False,
        type=float,
        default=None,
        help="GQ cutoff for an INDEL variant in low complexity region to be re-genotyped with DeepVariant. Variants with GQ below this will be re-genotyped."
    )
    pepper.add_argument(
        "--pepper_report_snp_above_freq",
        required=False,
        type=float,
        default=None,
        help="Report all SNPs above frequency for re-genotyping even if the predicted value is low. Set 0 to disable this."
    )
    pepper.add_argument(
        "--pepper_report_indel_above_freq",
        required=False,
        type=float,
        default=None,
        help="Report all INDELs above frequency for re-genotyping even if the predicted value is low. Set 0 to disable this."
    )

    margin = parser.add_argument_group('Arguments for Margin')
    margin.add_argument(
        "--margin_haplotag_model",
        type=str,
        help="Custom margin model."
    )
    margin.add_argument(
        "--margin_phase_model",
        type=str,
        help="Custom margin model."
    )

    deepvariant = parser.add_argument_group('Arguments for DeepVariant')
    deepvariant.add_argument(
        "--dv_model",
        default=None,
        type=str,
        help="Custom DeepVariant model. If this parameter is set, then one model will be used for both SNP and INDEL calling."
    )
    deepvariant.add_argument(
        "--dv_model_snp",
        default=None,
        type=str,
        help="Custom DeepVariant model for calling SNP calling."
    )
    deepvariant.add_argument(
        "--dv_model_indel",
        default=None,
        type=str,
        help="Custom DeepVariant model for calling INDEL calling."
    )
    deepvariant.add_argument(
        "--dv_alt_aligned_pileup",
        type=str,
        default=None,
        required=False,
        help="DeepVariant alt_align_pileup used for make_examples associated model: dv_model."
    )
    deepvariant.add_argument(
        "--dv_alt_aligned_pileup_snp",
        type=str,
        default=None,
        required=False,
        help="DeepVariant alt_align_pileup used for make_examples for snp calling associated model: dv_model_snp."
    )
    deepvariant.add_argument(
        "--dv_alt_aligned_pileup_indel",
        type=str,
        default=None,
        required=False,
        help="DeepVariant alt_align_pileup used for make_examples for indel calling associated model: dv_model_indel."
    )
    deepvariant.add_argument(
        "--dv_pileup_image_width",
        type=int,
        default=None,
        required=False,
        help="DeepVariant image width."
    )
    deepvariant.add_argument(
        "--dv_realign_reads",
        type=str,
        default=None,
        required=False,
        help="If true then local read alingment will be performed. [set: true/false]"
    )
    deepvariant.add_argument(
        "--dv_partition_size",
        type=int,
        default=None,
        required=False,
        help="DeepVariant partition_size used for make_examples."
    )
    deepvariant.add_argument(
        "--dv_min_mapping_quality",
        type=int,
        default=None,
        required=False,
        help="DeepVariant minimum mapping quality."
    )
    deepvariant.add_argument(
        "--dv_min_base_quality",
        type=int,
        default=None,
        required=False,
        help="DeepVariant minimum base quality."
    )
    deepvariant.add_argument(
        "--dv_vsc_min_fraction_indels",
        type=float,
        default=None,
        required=False,
        help="DeepVariant minimum vsc fraction for indels."
    )
    deepvariant.add_argument(
        "--dv_vsc_min_fraction_snps",
        type=float,
        default=None,
        required=False,
        help="DeepVariant minimum vsc fraction for snps."
    )
    deepvariant.add_argument(
        "--dv_sort_by_haplotypes",
        type=str,
        default=None,
        required=False,
        help="If true then haplotype sorting will be used. [set: true/false]"
    )
    deepvariant.add_argument(
        "--dv_parse_sam_aux_fields",
        type=str,
        default=None,
        required=False,
        help="If true then auxiliary field parsing is enabled. [set: true/false]"
    )
    deepvariant.add_argument(
        "--dv_add_hp_channel",
        type=str,
        default=None,
        required=False,
        help="If true then hp channel will be added. [set: true/false]"
    )
    deepvariant.add_argument(
        "--dv_use_hp_information",
        type=str,
        default=None,
        required=False,
        help="If true then hp information will be properly used. [set: true/false]"
    )
    deepvariant.add_argument(
        "--dv_use_multiallelic_mode",
        type=str,
        default=None,
        required=False,
        help="If true multiallelic model will be used during post-processing. [set: true/false]"
    )

    profile_group = required.add_mutually_exclusive_group(required=True)
    profile_group.add_argument("--ont_r9_guppy5_sup",
                               default=False,
                               action='store_true',
                               help="Set to call variants on R9.4.1 Guppy 5+/6+ sup/hac Oxford Nanopore reads.")
    profile_group.add_argument("--ont_r10_q20",
                               default=False,
                               action='store_true',
                               help="Set to call variants on R10.4 Q20 Oxford Nanopore reads.")
    profile_group.add_argument("--hifi",
                               default=False,
                               action='store_true',
                               help="Set to call variants on PacBio HiFi reads.")

    optional.add_argument(
        '-h',
        '--help',
        action='help',
        default=argparse.SUPPRESS,
        help='Show this text and exit.'
    )

    return parser


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s: string holding boolean value
    :return:
    """
    if s.lower() not in {'false', 'true', '1', 't', '0', 'f'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true' or s.lower() == 't' or s.lower() == '1'


def main():
    """
    Main interface for PEPPER-MARGIN-DeepVariant. The submodules supported as of now are these:
    1) Variant call
    2) Phase
    3) Polish
    """
    parser = argparse.ArgumentParser(description='Run PEPPER-Margin-DeepVariant for variant calling.\n'
                                                 'Example run: run_pepper_margin_deepvariant -b <BAM> -f <FASTA> -o <OUTPUT_DIR> -t <THREADS> <--ont_r9_guppy5_sup/--ont_r10_q20/--hifi>'
                                                 'OUTPUT FILE DESCRIPTION:\n'
                                                 '1) OUTPUT_DIR/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz: Final output file. Filename can be changed with -p option.\n'
                                                 '2) OUTPUT_DIR/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.visual_report.html: Basic variant statistics of the output file.\n'
                                                 'Intermediate output files:\n'
                                                 '1) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz: Full set of variants from PEPPER.\n'
                                                 '2) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz: Variants from PEPPER that will NOT be re-genotyped with DeepVariant.\n'
                                                 '3) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz: All variants from PEPPER that will be re-genotyped with DeepVariant.\n'
                                                 '4) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs.vcf.gz: SNP variants from PEPPER that will be re-genotyped with DeepVariant.\n'
                                                 '5) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_INDEL.vcf.gz: INDEL variants from PEPPER that will be re-genotyped with DeepVariant.\n'
                                                 '6) OUTPUT_DIR/intermediate_files/DEEPVARIANT_OUTPUT_SNP.vcf.gz: SNP variants reported by PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs re-genotyped with DeepVariant\n'
                                                 '7) OUTPUT_DIR/intermediate_files/DEEPVARIANT_OUTPUT_INDEL.vcf.gz: INDEL variants reported by PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs re-genotyped with DeepVariant\n'
                                                 '8) OUTPUT_DIR/intermediate_files/DEEPVARIANT_OUTPUT.vcf.gz: Variants reported by PEPPER_VARIANT_OUTPUT_VARIANT_CALLING re-genotyped with DeepVariant\n'
                                                 'Note 1: DEEPVARIANT_OUTPUT_SNP.vcf.gz and DEEPVARIANT_OUTPUT_INDEL.vcf.gz is subset of PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz.\n'
                                                 'Note 2: DEEPVARIANT_OUTPUT.vcf.gz only exists if we use PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz as candidates and do not do SNP and INDEL calling separately.\n'
                                                 'Note 3: The final PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz file is generated by merging PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz and DEEPVARIANT_OUTPUT* files.\n',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )

    subparsers = parser.add_subparsers(dest='sub_command')
    subparsers.required = False

    parser_call_variant = subparsers.add_parser('call_variant', description='Run PEPPER-Margin-DeepVariant for variant calling.\n'
                                                                            'Example run: run_pepper_margin_deepvariant -b <BAM> -f <FASTA> -o <OUTPUT_DIR> -t <THREADS> <--ont_r9_guppy5_sup/--ont_r10_q20/--hifi>\n'
                                                                            'Output file description:\n'
                                                                            '1) OUTPUT_DIR/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz: Final output file. Filename can be changed with -p option.\n'
                                                                            '2) OUTPUT_DIR/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.visual_report.html: Basic variant statistics of the output file.\n'
                                                                            'Intermediate output files:\n'
                                                                            '1) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_FULL.vcf.gz: Full set of variants from PEPPER.\n'
                                                                            '2) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz: Variants from PEPPER that will NOT be re-genotyped with DeepVariant.\n'
                                                                            '3) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz: All variants from PEPPER that will be re-genotyped with DeepVariant.\n'
                                                                            '4) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs.vcf.gz: SNP variants from PEPPER that will be re-genotyped with DeepVariant.\n'
                                                                            '5) OUTPUT_DIR/intermediate_files/PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_INDEL.vcf.gz: INDEL variants from PEPPER that will be re-genotyped with DeepVariant.\n'
                                                                            '6) OUTPUT_DIR/intermediate_files/DEEPVARIANT_OUTPUT_SNP.vcf.gz: SNP variants reported by PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs re-genotyped with DeepVariant\n'
                                                                            '7) OUTPUT_DIR/intermediate_files/DEEPVARIANT_OUTPUT_INDEL.vcf.gz: INDEL variants reported by PEPPER_VARIANT_OUTPUT_VARIANT_CALLING_SNPs re-genotyped with DeepVariant\n'
                                                                            '8) OUTPUT_DIR/intermediate_files/DEEPVARIANT_OUTPUT.vcf.gz: Variants reported by PEPPER_VARIANT_OUTPUT_VARIANT_CALLING re-genotyped with DeepVariant\n'
                                                                            'Note 1: DEEPVARIANT_OUTPUT_SNP.vcf.gz and DEEPVARIANT_OUTPUT_INDEL.vcf.gz is subset of PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz.\n'
                                                                            'Note 2: DEEPVARIANT_OUTPUT.vcf.gz only exists if we use PEPPER_VARIANT_OUTPUT_VARIANT_CALLING.vcf.gz as candidates and do not do SNP and INDEL calling separately.\n'
                                                                            'Note 3: The final PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz file is generated by merging PEPPER_VARIANT_OUTPUT_PEPPER.vcf.gz and DEEPVARIANT_OUTPUT* files.\n'

                                                , add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    add_variant_call_arguments(parser_call_variant)

    options, unparsed = parser.parse_known_args()

    if options.sub_command == 'call_variant':
        sys.stderr.write("[" + str(datetime.now().strftime('%m-%d-%Y %H:%M:%S')) + "] INFO: VARIANT CALLING MODULE SELECTED\n")
        sys.stderr.flush()
        run_variant_calling(options)
    elif options.version is True:
        print("VERSION: ", __version__)
    else:
        sys.stderr.write("ERROR: NO SUBCOMMAND SELECTED. PLEASE SELECT ONE OF THE AVAILABLE SUB-COMMANDS.\n")
        parser.print_help()
        print("Example run: run_pepper_margin_deepvariant -b <BAM> -f <FASTA> -o <OUTPUT_DIR> -t <THREADS> <--ont_r9_guppy5_sup/--ont_r10_q20/--hifi>\n")


if __name__ == '__main__':
    main()

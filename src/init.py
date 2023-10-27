#!/usr/bin/env python 
"""
Project: TRgenotype
Script: paras.py
Author: Peng Jia
E-mail: pengjia@stu.xjtu.edu.cn
Time : 2023/5/31
Description: TODO
"""

import argparse
import os
import sys


curpath = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.append(os.path.dirname(curpath))
from src.units import *
from src.region import Region
from src.repeat import Repeat


def args_process():
    """
    argument procress
    """
    # logger.info(" ".join(sys.argv))
    defaultPara = get_value("default")
    commands = []
    commandsParser = {}
    tool_name=get_value("tools_name")
    parser = argparse.ArgumentParser(description=f'{tool_name}: Tandem repeat analysis toolbox.'
                                     # + ".show help of subcommand with '"
                                     # + get_value("tools_name") + " <subcommand> -h'"
                                     )
    parser.usage =f"{tool_name} <command> [options]"
    parser.add_argument('-V', '--version', action='version',
                        version=get_value("tools_name") + get_value("tools_version"))
    subparsers = parser.add_subparsers(title="command", metavar="", dest='command')

    ###################################################################################################################
    # add arguments for genotype module
    parser_gt = subparsers.add_parser('genotype', help='Tandem repeat genotyping')
    parser_gt.description = 'Tandem repeat genotyping.'
    commands.append("genotype")
    defaultPara_gt = defaultPara["genotype"]
    ##################################################################################
    # group input and output
    input_and_output = parser_gt.add_argument_group(title="Input and output")
    input_and_output.add_argument('-i', '--input', required=True, type=str, nargs=1,
                                  help="The path of input bam/cram file, 'HP' tag is recommended for long reads[required]")
    input_and_output.add_argument('-r', '--repeat', required=True, type=str, nargs=1,
                                  help="The path of the repeat regions, e.g. repeat.bed [required]")
    input_and_output.add_argument('-o', '--output', required=True, type=str, nargs=1,
                                  help="The path of output file prefix [required]")
    # input_and_output.add_argument('-ref', '--reference', required=False, type=str, nargs=1,
    #                               help="The path of reference file [required]")
    input_and_output.add_argument('-tech', '--technology', type=str, nargs=1, choices=["hifi", "clr", "ont", "ilm"],
                                  required=False,default="hifi",
                                  help='Sequencing technology [required]')
    # input_and_output.add_argument('-tech', '--technology', type=str, nargs=1, choices=["ccs", "clr", "ont", "ilm"],
    #                               default=[defaultPara_gt["tech"]],
    #                               help='Sequencing technology [default:'
    #                                    + str(defaultPara_gt["tech"]) + ']')
    input_and_output.add_argument('-s', '--sample', type=str, nargs=1,
                                  # default=[defaultPara_gt["hap"]],
                                  default=["default"],
                                  help=" sample name in output vcf file [default: extract from bam file]")
    # input_and_output.add_argument("-mf", '--microsatellite_region_format', type=str, nargs=1,
    #                               choices=["bed", "json", "msisensor_scan"],
    #                               default=[defaultPara_gt["microsatellite_region_format"]],
    #                               help='Input format of microsatellites region file [default:'
    #                                    + str(defaultPara_gt["microsatellite_region_format"]) + ']')
    ##################################################################################
    # group analysis regions
    # general_realign = parser_gt.add_argument_group(title="Analysis regions")

    # general_realign.add_argument('-ks', '--kmer_size', type=int, nargs=1,
    #                              default=[defaultPara_gt["kmer_size"]],
    #                              help="Debug mode for developers [default:" +
    #                                   str(defaultPara_gt["kmer_size"]) + "]")

    ##################################################################################
    # group general option

    general_option = parser_gt.add_argument_group(title="General option")
    # general_option.add_argument('-pl', '--prefix_len', type=int, nargs=1,
    #                             default=[defaultPara_gt["prefix_len"]],
    #                             help="[prefix_len] bp upstream for mutation analysis [default:" +
    #                                  str(defaultPara_gt["prefix_len"]) + "]")
    # general_option.add_argument('-sl', '--suffix_len', type=int, nargs=1,
    #                             default=[defaultPara_gt["suffix_len"]],
    #                             help="[suffix_len] bp downstream for mutation analysis  [default:" +
    #                                  str(defaultPara_gt["suffix_len"]) + "]")
    general_option.add_argument('-d', '--debug', type=str, nargs=1, choices=["True", "False"],
                                default=[defaultPara_gt["debug"]],
                                help="Debug mode for developers [default:" +
                                     str(defaultPara_gt["debug"]) + "]")
    # general_option.add_argument('-e', '--only_exact_repeat', type=str, nargs=1, choices=["True", "False"],
    #                             default=[defaultPara_gt["only_exact_repeat"]],
    #                             help="Only analyze exact repeat regions [default:"
    #                                  + str(defaultPara_gt["only_exact_repeat"]) + "]")
    # general_option.add_argument('-ih', '--ignore_homopolymer', type=str, nargs=1, choices=["True", "False"],
    #                             default=[defaultPara_gt["ignore_homopolymer"]],
    #                             help="Ignore homopolymer regions [default:"
    #                                  + str(defaultPara_gt["ignore_homopolymer"]) + "]")
    # general_option.add_argument('-up', '--using_phasing_info', type=int, nargs=1, choices=[True, False],
    #                             default=[defaultPara_gt["using_phasing_info"]],
    #                             help="Use phased information from ccs/HiFi reads "
    #                                  "to genotype and phase microsatellite allele [default:"
    #                                  + str(defaultPara_gt["using_phasing_info"]) + "]")

    # general_option.add_argument("-minr", '--minimum_repeat_times',
    #                             default=[defaultPara_gt["minimum_repeat_times"]],
    #                             type=str, nargs=1,
    #                             help="Minimum repeat times of microsatellites [default:"
    #                                  + defaultPara_gt["minimum_repeat_times"] + "]")
    # general_option.add_argument('-maxr', '--maximum_repeat_times',
    #                             default=[defaultPara_gt["maximum_repeat_times"]], type=str, nargs=1,
    #                             help="Maximum repeat times of microsatellites [default:"
    #                                  + defaultPara_gt["maximum_repeat_times"] + "]")
    # general_option.add_argument('-minh', '--minimum_phasing_reads',
    #                             default=[defaultPara_gt["minimum_phasing_reads"]], type=str, nargs=1,
    #                             help="Minimum reads for each haplotype reporting [default:"
    #                                  + str(defaultPara_gt["minimum_phasing_reads"]) + "]")
    general_option.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1,
                                default=[defaultPara_gt["minimum_mapping_quality"]],
                                help="minimum mapping quality of read for analysis [default:" +
                                     str(defaultPara_gt["minimum_mapping_quality"]) + "]")
    general_option.add_argument('-ms', '--minimum_support_reads', type=int, nargs=1,
                                default=[defaultPara_gt["minimum_support_reads"]],
                                help="minimum support reads of an available call[default:" +
                                     str(defaultPara_gt["minimum_support_reads"]) + "]")
    # general_option.add_argument('-md', '--maximum_distance_of_two_complex_events', type=int, nargs=1,
    #                             default=[defaultPara_gt["maximum_distance_of_two_complex_events"]],
    #                             help="minimum support reads of an available microsatellite call[default:" +
    #                                  str(defaultPara_gt["maximum_distance_of_two_complex_events"]) + "]")
    general_option.add_argument('-a', '--min_allele_fraction', type=int, nargs=1,
                                default=[defaultPara_gt["min_allele_fraction"]],
                                help="minimum allele fraction report [default:" +
                                     str(defaultPara_gt["min_allele_fraction"]) + "]")
    multiple_thread = parser_gt.add_argument_group(title="Multiple thread")
    multiple_thread.add_argument('-t', '--threads', type=int, nargs=1,
                                 default=[defaultPara_gt["threads"]],
                                 help="The number of  threads to use [default:" +
                                      str(defaultPara_gt["threads"]) + "]")
    multiple_thread.add_argument('-b', '--batch', type=int, nargs=1,
                                 default=[defaultPara_gt["batch"]],
                                 help="The number of repeat one thread process [default:" +
                                      str(defaultPara_gt["batch"]) + "]")
    commandsParser["genotype"] = parser_gt
    # general_option.add_argument('--sequencing_error', type=int, nargs=1,
    #                             default=[defaultPara_gt["sequencing_error"]],
    #                             help="Sequencing error, allele frequency less than SEQUENCING ERROR with be remove "
    #                                  "[default:" + str(defaultPara_gt["sequencing_error"]) + "]")

    ##################################################################################

    # add arguments for genotype module
    parser_train = subparsers.add_parser('train', help='Heterozygosity Discrimination Model Training')
    parser_train.description = 'Heterozygosity Discrimination Model Training.'
    commands.append("train")
    defaultPara_gt = defaultPara["train"]
    ##################################################################################
    # group input and output
    input_and_output = parser_train.add_argument_group(title="Input and output")
    input_and_output.add_argument('-i', '--input', required=True, type=str, nargs=1,
                                  help="The path of input bam/cram file, 'HP' tag is required for long reads[required]")
    input_and_output.add_argument('-r', '--repeat', required=True, type=str, nargs=1,
                                  help="The path of the repeat regions, e.g. repeat.bed [required]")
    input_and_output.add_argument('-o', '--output', required=True, type=str, nargs=1,
                                  help="The path of output model e.g. output.model.pkl [required]")
    # input_and_output.add_argument('-ref', '--reference', required=False, type=str, nargs=1,
    #                               help="The path of reference file [required]")
    input_and_output.add_argument('-tech', '--technology', type=str, nargs=1, choices=["hifi", "clr", "ont", "ilm"],
                                  required=True, default="hifi",
                                  help='Sequencing technology [required]')
    # input_and_output.add_argument('-tech', '--technology', type=str, nargs=1, choices=["ccs", "clr", "ont", "ilm"],
    #                               default=[defaultPara_gt["tech"]],
    #                               help='Sequencing technology [default:'
    #                                    + str(defaultPara_gt["tech"]) + ']')
    input_and_output.add_argument('-s', '--sample', type=str, nargs=1,
                                  # default=[defaultPara_gt["hap"]],
                                  default=["default"],
                                  help=" sample name in output vcf file [default: extract from bam file]")

    general_option = parser_train.add_argument_group(title="General option")

    general_option.add_argument('-d', '--debug', type=str, nargs=1, choices=["True", "False"],
                                default=[defaultPara_gt["debug"]],
                                help="Debug mode for developers [default:" +
                                     str(defaultPara_gt["debug"]) + "]")

    general_option.add_argument('-q', '--minimum_mapping_quality', type=int, nargs=1,
                                default=[defaultPara_gt["minimum_mapping_quality"]],
                                help="minimum mapping quality of read for analysis [default:" +
                                     str(defaultPara_gt["minimum_mapping_quality"]) + "]")
    general_option.add_argument('-ms', '--minimum_support_reads', type=int, nargs=1,
                                default=[defaultPara_gt["minimum_support_reads"]],
                                help="minimum support reads of an available call[default:" +
                                     str(defaultPara_gt["minimum_support_reads"]) + "]")
    # general_option.add_argument('-md', '--maximum_distance_of_two_complex_events', type=int, nargs=1,
    #                             default=[defaultPara_gt["maximum_distance_of_two_complex_events"]],
    #                             help="minimum support reads of an available microsatellite call[default:" +
    #                                  str(defaultPara_gt["maximum_distance_of_two_complex_events"]) + "]")
    general_option.add_argument('-a', '--min_allele_fraction', type=int, nargs=1,
                                default=[defaultPara_gt["min_allele_fraction"]],
                                help="minimum allele fraction report [default:" +
                                     str(defaultPara_gt["min_allele_fraction"]) + "]")


    # group for bam2dis
    # bam2dis_option = parser_gt.add_argument_group(title="Option for bam2dis")

    # bam2dis_option.add_argument('-am', '--allow_mismatch', type=bool, nargs=1, choices=[True, False],
    #                             default=[defaultPara_gt["allow_mismatch"]],
    #                             help="allow mismatch when capture microsatellite [default:"
    #                                  + str(defaultPara_gt["allow_mismatch"]) + "]")

    ##################################################################################
    # group for multiple_thread

    multiple_thread = parser_train.add_argument_group(title="Multiple thread")
    multiple_thread.add_argument('-t', '--threads', type=int, nargs=1,
                                 default=[defaultPara_gt["threads"]],
                                 help="The number of  threads to use [default:" +
                                      str(defaultPara_gt["threads"]) + "]")
    multiple_thread.add_argument('-b', '--batch', type=int, nargs=1,
                                 default=[defaultPara_gt["batch"]],
                                 help="The number of repeat one thread process [default:" +
                                      str(defaultPara_gt["batch"]) + "]")
    commandsParser["train"] = parser_train
    ###################################################################################################################

    ###################################################################################################################
    if len(os.sys.argv) < 2:
        parser.print_help()
        return False

    if os.sys.argv[1] in ["-h", "--help", "-?"]:
        parser.print_help()
        return False
    if os.sys.argv[1] in ["-V", "-v", "--version"]:
        # parser.print_help()
        parser.parse_args()
        return False
    if os.sys.argv[1] not in commands:
        logger.error("Command Error! " + os.sys.argv[1] +
                     "is not the available command.\n"
                     "[Tips] Please input correct command such as " + ", ".join(commands) + "!")
        parser.print_help()
        # parser.parse_args()
        return False
    if len(os.sys.argv) == 2 and (os.sys.argv[1] in commands):
        commandsParser[os.sys.argv[1]].print_help()
        return False
    return parser

def args_init(args):
    """
    argument procress
    """

    paras = {}
    paras["input"] = args.input[0]
    paras["output"] = args.output[0]
    paras["repeat"] = args.repeat[0]
    paras["tech"] = args.technology[0]
    paras["debug"] = True if args.debug[0].lower() == "true" else False
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    paras["threads"] = args.threads[0]
    paras["batch"] = args.batch[0]
    paras["sample"] = args.sample[0]
    error_stat = False
    if os.path.exists(paras["input"]):
        logger.info("The input file is : " + paras["input"] + ".")
    else:
        logger.error('The input file ' + paras["input"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["repeat"]):
        logger.info("The microsatellites file  is : " + paras["repeat"])
    else:
        logger.error('The microsatellites file ' + paras["repeat"] + ' is not exist, please check again')
        error_stat = True
    if paras["input"][-4:] == "cram":
        paras["input_format"] = "cram"
        cramfile = pysam.AlignmentFile(paras["input"], mode="rb", reference_filename=paras["reference"])
        if not cramfile.has_index():
            logger.info("Build index for the input cram ...")
            pysam.index(paras["input"])
        cramfile.close()
    else:
        paras["input_format"] = "bam"
        bamfile = pysam.AlignmentFile(paras["input"], mode="rb")
        if not bamfile.has_index():
            logger.info("Build index for the input bam ...")
            pysam.index(paras["input"])
        bamfile.close()
    if not os.path.exists(paras["output"]):
        pass
    else:
        if paras["debug"]:
            pass
        else:
            logger.error('The output ' + paras["output"] +
                         ' is still exist! in case of overwrite files in this workspace, '
                         'please check your script!')
            error_stat = True
    if error_stat:
        return False
    logger.info("The output is : " + paras["output"] + ".")
    input_path = paras["input"]
    input_path = input_path[:-1] if input_path[-1] == "/" else input_path
    if paras["sample"] == "default":
        case = input_path.split("/")[-1].strip(".bam")
        case = case.strip(".cram")
    else:
        case = paras["sample"]
    paras["output_details"] = paras["output"] + ".details.info"
    paras["output_info"] = paras["output"]
    paras["depth_dict"] = get_reads_depth(paras["input"])
    set_value("case", case)
    set_value("paras", paras)
    return True
def extract_repeat_info(paras):
    logger.info(f"Exacting repeat from {paras['repeat']} ...")
    my_threads = int(paras["threads"])
    my_batch = int(paras["batch"])

    repeat_infos = {}
    repeat_infos_sorted = {}
    repeat_info_num = {}
    for line in open(paras["repeat"]):
        chrom, start, end = line[:-1].split("\t")[:3]
        start = int(start)

        if chrom not in repeat_infos:
            repeat_infos[chrom] = {}
        repeat_infos[chrom][start] = line

    for chrom, info in repeat_infos.items():

        repeat_infos_sorted[chrom] = []
        start_sorted = sorted([i for i in info])
        chunk = []
        for idx, start in enumerate(start_sorted, 1):
            if idx % my_batch == 0:
                repeat_infos_sorted[chrom].append(Region(chunk, threads=my_threads))
                chunk = []
            else:
                chunk.append(info[start])
        else:
            if len(chunk) > 0:
                repeat_infos_sorted[chrom].append(Region(chunk, threads=my_threads))
        repeat_num = idx
        repeat_info_num[chrom] = repeat_num
        logger.info(f"{chrom}: {repeat_num} repeats.")

    # print([i for i in repeat_info_num])
    total_repeat = sum([i for j, i in repeat_info_num.items()])
    logger.info(f'Total: {total_repeat} repeats.')
    # print(repeat_infos_sorted)
    set_value("total_repeat_num", total_repeat)
    set_value("repeat_info", repeat_infos_sorted)
    set_value("chrom_repeat_num", repeat_info_num)
    return 1
def get_reads_depth(path_bam, min_len=10000000, sample_point_each_contig=2):
    depths = []
    bam = pysam.Samfile(path_bam, threads=4)
    contigs = {ctg: ctg_length for ctg_length, ctg in zip(bam.header.lengths, bam.header.references) if
               ctg_length > min_len}
    np.random.seed(1)
    logger.info("Calculate sequencing depth of this sample...")
    for ctg, ctg_length in contigs.items():
        points = np.random.randint(50000, ctg_length - 50000, sample_point_each_contig)
        dps = [len([i for i in bam.fetch(ctg, i - 1, stop=i)]) for i in points]
        sorted_arr = np.sort(dps)
        remove_count = int(sample_point_each_contig * 0.05)
        if remove_count > 0:
            trimmed_dps = sorted_arr[remove_count:-remove_count]
        else:
            trimmed_dps = sorted_arr
        depths.extend(trimmed_dps)
    mean, median, std = np.mean(depths), np.median(depths), np.std(depths)
    q1 = np.percentile(depths, 25)
    q3 = np.percentile(depths, 75)
    iqr = q3 - q1
    depths_dict = {"mean": mean, "median": median, "std": std,
                   "sigma_min": mean - 3 * std, "sigma_max": mean + 3 * std,
                   "q1": q1, "q3": q3, "iqr_min": q1 - 1.5 * iqr, "iqr_max": q3 + 1.5 * iqr
                   }
    return depths_dict
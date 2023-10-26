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
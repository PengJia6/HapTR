#!/usr/bin/env python 
"""
Project: HapTR
Script: param.py
Author: Peng Jia
E-mail: pengjia@stu.xjtu.edu.cn
Time : 2023/10/27
Description: TODO
"""
from src.units import *
import argparse
import os
import sys
import pysam


class Param:
    def __int__(self, command):
        self.command = command
        pass

    def args_cmd(self):
        defaultPara = get_value("default")
        commands = []
        commandsParser = {}
        tool_name = get_value("tools_name")
        parser = argparse.ArgumentParser(description=f'{tool_name}: Tandem repeat analysis toolbox.'
                                         # + ".show help of subcommand with '"
                                         # + get_value("tools_name") + " <subcommand> -h'"
                                         )
        parser.usage = f"{tool_name} <command> [options]"
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
        input_and_output.add_argument('-ref', '--reference', required=False, type=str, nargs=1, default=[""],
                                      help="The path of reference file [required]")
        input_and_output.add_argument('-tech', '--technology', type=str, nargs=1, choices=["hifi", "clr", "ont", "ilm"],
                                      required=False, default="hifi",
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
        general_option.add_argument('-f', '--flank_size', type=int, nargs=1,
                                    default=[defaultPara_gt["flank_size"]],
                                    help="minimum mapping quality of read for analysis [default:" +
                                         str(defaultPara_gt["flank_size"]) + "]")
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
        input_and_output.add_argument('--output_info', required=True, type=str, nargs=1,
                                      help="The path of model details. e.g. output.model.txt [required]")
        input_and_output.add_argument('-ref', '--reference', required=True, type=str, nargs=1, default=[""],
                                      help="The path of reference file")
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
        general_option.add_argument('--min_phased_ratio', type=float, nargs=1,
                                    default=[defaultPara_gt["min_phased_ratio"]],
                                    help=f"minimum ratio of phased reads for {tool_name} train [default:"
                                         f"{defaultPara_gt['min_phased_ratio']}]")
        general_option.add_argument('--min_phased_reads', type=int, nargs=1,
                                    default=[defaultPara_gt["min_phased_reads"]],
                                    help=f"minimum phased reads for {tool_name} train [default:"
                                         f"{defaultPara_gt['min_phased_reads']}]")

        general_option.add_argument('-f', '--flank_size', type=int, nargs=1,
                                    default=[defaultPara_gt["flank_size"]],
                                    help="minimum mapping quality of read for analysis [default:" +
                                         str(defaultPara_gt["flank_size"]) + "]")

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
        multiple_thread.add_argument('-rs', '--region_size', type=int, nargs=1,
                                     default=[defaultPara_gt["region_size"]],
                                     help="The number of repeat one region process [default:" +
                                          str(defaultPara_gt["region_size"]) + "]")
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
        self.args = parser.parse_args()
        self.command = self.args.command
        return True

    def args_init(self):
        """
        argument procress
        """
        self.input_bam_path = self.args.input[0]
        self.input_bam_type = "bam" if self.input_bam_path.endswith(".bam") else (
            "cram" if self.input_bam_path.endswith(".cram") else "error")
        self.reference_path = self.args.reference[0]
        if self.command == "train":
            self.output_model = self.args.output[0]
            self.output_info = self.args.output_info[0]
            self.output = {"model": self.output_model, "info": self.output_info}
            self.min_phased_ratio=self.args.min_phased_ratio[0]
            self.min_phased_reads=self.args.min_phased_reads[0]


        else:
            self.output = {self.output_model}

        self.repeat_bed = self.args.repeat[0]
        self.tech = self.args.technology[0]
        self.debug = True if self.args.debug[0].lower() == "true" else False
        self.minimum_support_reads = self.args.minimum_support_reads[0]
        self.minimum_mapping_quality = self.args.minimum_mapping_quality[0]
        self.threads = self.args.threads[0]
        self.batch = self.args.batch[0]
        self.region_size = self.args.region_size[0]
        self.flank_size = self.args.flank_size[0]

        self.sample = self.args.sample[0]
        error_stat = False
        if os.path.exists(self.input_bam_path):
            logger.info(f"The input file is : {self.input_bam_path} .")
        else:
            logger.error(f'The input file {self.input_bam_path} is not exist, please check again')
            error_stat = True
        if self.input_bam_path == "cram":
            if not os.path.exists(self.reference_path):
                logger.error(f'Please input available reference file (*.fa or *.fasra)')
                error_stat = True
            else:
                if not os.path.exists(f"{self.input_bam_path}.crai"):
                    logger.info("Build index for the input cram ...")
                    pysam.index(self.input_bam_path)
        elif self.input_bam_path == "bam":
            if not os.path.exists(f"{self.input_bam_path}.bai"):
                logger.info("Build index for the input bam ...")
                pysam.index(f"{self.input_bam_path}")

        if os.path.isfile(self.repeat_bed):
            logger.info(f"The repeat file  is : {self.repeat_bed}")
        else:
            logger.error(f'The repeat file {self.repeat_bed} is not exist, please check again')
            error_stat = True
        if os.path.isfile(self.reference_path):
            logger.info(f"The reference file  is : {self.reference_path}")
        else:
            logger.error(f'The reference file {self.reference_path} is not exist, please check again')
            error_stat = True
        output_path = []
        for j, i in self.output.items():
            if i in output_path:
                logger.error(
                    (f'The output {i} is  used in your command ! in case of overwrite files in this workspace, '
                     'please check your command!'))
                output_path.append(i)
                error_stat = True
            if os.path.exists(i):
                if not self.debug:
                    logger.error(f'The output {i} is still exist! in case of overwrite files in this workspace, '
                                 'please check your command!')
                    error_stat = True

            logger.info(f"The output {j} is : {i}.")
        if error_stat:
            return False
        self.depths_dict = self.get_reads_depth()
        return True

    def get_reads_depth(self, min_len=10000000, sample_point_each_contig=2):
        depths = []
        bam = pysam.Samfile(self.input_bam_path, threads=4)
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

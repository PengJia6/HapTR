#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: MSHunter
# Script : genotype.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
import numpy as np

from src.units import *
from src.repeat import *
import os
import pysam
import multiprocessing

os.environ['OPENBLAS_NUM_THREADS'] = '1'


def get_reads_depth(path_bam, min_len=10000000, sample_point_each_contig=10):
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


def genotype_init(args):
    """
    argument procress
    """
    paras = {}
    paras["input"] = args.input[0]
    paras["output"] = args.output[0]
    paras["microsatellite"] = args.microsatellite[0]

    # paras["reference"] = args.reference[0]
    # paras["microsatellite_region_format"] = args.microsatellite_region_format[0]
    paras["tech"] = args.technology[0]
    # paras["hap"] = args.haplotype_bam[0]
    # paras["prefix_len"] = args.prefix_len[0]
    # paras["suffix_len"] = args.suffix_len[0]
    paras["debug"] = True if args.debug[0].lower() == "true" else False
    paras["only_homopolymer"] = True if str(args.only_homopolymers[0]).lower() == "true" else False
    paras["only_simple"] = True if args.only_simple[0].lower() == "true" else False
    # paras["using_phasing_info"] = args.using_phasing_info[0]
    paras["minimum_support_reads"] = args.minimum_support_reads[0]
    paras["minimum_mapping_quality"] = args.minimum_mapping_quality[0]
    # paras["allow_mismatch"] = args.allow_mismatch[0]
    paras["threads"] = args.threads[0]
    paras["batch"] = args.batch[0]
    # paras["minimum_phasing_reads"] = args.minimum_phasing_reads[0]
    # paras["min_allele_fraction"] = args.min_allele_fraction[0]
    # paras["sequencing_error"] = args.sequencing_error[0]
    # paras["maximum_distance_of_two_complex_events"] = args.maximum_distance_of_two_complex_events[0]
    paras["sample"] = args.sample[0]
    # paras["ranges_of_repeat_times"] = {}

    # for i in args.minimum_repeat_times[0].split(";"):
    #     unitRange, repeatRange = i.split(":")
    #     if "-" in unitRange:
    #         unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
    #     else:
    #         unitStart = int(unitRange)
    #         unitEnd = unitStart
    #     repeatStart = int(repeatRange)
    #     # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
    #     for ur in range(unitStart, unitEnd + 1):
    #         if ur not in paras["ranges_of_repeat_times"]:
    #             paras["ranges_of_repeat_times"][ur] = {}
    #         paras["ranges_of_repeat_times"][ur]["min"] = repeatStart
    #     for i in args.maximum_repeat_times[0].split(";"):
    #         # print(i)
    #         unitRange, repeatRange = i.split(":")
    #         if "-" in unitRange:
    #             unitStart, unitEnd = tuple(map(int, unitRange.split("-")))
    #         else:
    #             unitStart = int(unitRange)
    #             unitEnd = unitStart
    #         repeatStart = int(repeatRange)
    #         # print(unitStart,unitEnd,"  ",repeatStart, repeatEnd)
    #         for ur in range(unitStart, unitEnd + 1):
    #             if ur not in paras["ranges_of_repeat_times"]:
    #                 paras["ranges_of_repeat_times"][ur] = {}
    #             paras["ranges_of_repeat_times"][ur]["max"] = repeatStart
    error_stat = False
    if os.path.exists(paras["input"]):
        logger.info("The input file is : " + paras["input"] + ".")
    else:
        logger.error('The input file ' + paras["input"] + ' is not exist, please check again')
        error_stat = True
    if os.path.isfile(paras["microsatellite"]):
        logger.info("The microsatellites file  is : " + paras["microsatellite"])
    else:
        logger.error('The microsatellites file ' + paras["microsatellite"] + ' is not exist, please check again')
        error_stat = True
    # if os.path.isfile(paras["reference"]):
    #     logger.info("The reference file is : " + paras["reference"] + ".")
    # else:
    #     paras["reference"] = "" if paras["reference"] == "." else paras["reference"]
    #     logger.error('The reference file ' + paras["reference"] + ' is not exist, please check again')
    #     error_stat = True
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
    # output_path = paras["output"]
    # output_path = output_path if output_path[-1] == "/" else output_path + "/"

    # if not os.path.exists(output_path):
    #     os.makedirs(output_path)
    # paras["output"] = output_path
    input_path = paras["input"]
    input_path = input_path[:-1] if input_path[-1] == "/" else input_path
    if paras["sample"] == "default":
        case = input_path.split("/")[-1].strip(".bam")
        case = case.strip(".cram")
    else:
        case = paras["sample"]
    # paras["output_pre"] = paras["output"] + case + ".pre.vcf.gz"
    # paras["output_tmp"] = paras["output"] + case + "_tmp"
    # if not os.path.exists(paras["output_tmp"]):
    #     os.makedirs(paras["output_tmp"])
    #
    paras["output_details"] = paras["output"] + ".details.info"
    paras["output_info"] = paras["output"]
    # paras["output_micro"] = paras["output"] + case + "_micro.vcf.gz"
    # paras["output_indel"] = paras["output"] + case + "_indel.vcf.gz"
    # paras["output_snv"] = paras["output"] + case + "_snv.vcf.gz"
    # paras["output_complex"] = paras["output"] + case + "_complex.vcf.gz"
    paras["depth_dict"] = get_reads_depth(paras["input"])
    set_value("case", case)
    set_value("paras", paras)
    # print(paras)
    #
    # print(paras)
    # exit()
    return True


def run_chunk(trs, thread, fun):
    # print(thread)
    thread = thread if thread < 3 else (thread * 0.9)
    pool = multiprocessing.Pool(processes=int(thread))
    # print("fun",fun)
    # print("trs",trs)
    res = pool.map(fun, trs)
    pool.close()
    pool.join()
    return res


def process_one_site(line):
    line, paras = line
    path_bam = paras["input"]
    line_info = line[:-1].split("\t")
    chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content, repeat_type, \
        repeat_subtype, source, site_id, complex_repeat, annotation, up, down = line_info

    repeat = Repeat(chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content,
                    repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down,
                    depth=paras["depth_dict"])
    repeat.process_reads(path_bam)

    # reads = [i for i in pysam.Samfile(f"{path_bam}").fetch(chrom, int(start), int(end))]

    # return [repeat.get_output_len_info(), repeat.get_output_details(), repeat.get_info_debug()]
    return [repeat.get_output_len_info(), "", ""]


def process_one_site_only_dis(line):
    line, paras = line
    path_bam = paras["input"]
    line_info = line[:-1].split("\t")
    chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content, repeat_type, \
        repeat_subtype, source, site_id, complex_repeat, annotation, up, down = line_info

    repeat = Repeat(chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content,
                    repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down,
                    depth=paras["depth_dict"])
    repeat.process_reads(path_bam)

    # reads = [i for i in pysam.Samfile(f"{path_bam}").fetch(chrom, int(start), int(end))]

    return [repeat.get_output_len_info(), repeat.get_output_details(), repeat.get_info_debug()]


def write_output():
    pass


def genotype(parase):
    if not genotype_init(parase):
        logger.error("Genotype init ERROR!")
        return -1
    paras = get_value("paras")
    my_threads = int(paras["threads"])
    my_batch = int(paras["batch"])
    total_number = len(open(paras["microsatellite"]).readlines())
    lines = []
    line_num = 0
    file_info = open(f"{paras['output_info']}", "w")
    file_info_detail = open(f"{paras['output_details']}", "w")
    file_info_debug = open(f"{paras['output_info']}.debug", "w")
    debug_num = 0
    for line in open(paras["microsatellite"]):
        line_num += 1
        # if line_num < 1000: continue
        # if line_num > 1100: break
        # if "chr1" not in line or "736186" not in line:
        #     continue
        #     # chr1_736186_736195
        # print(line)
        # if line_num <=249867: continue

        lines.append([line, paras])

        if line_num % (my_threads * my_batch) == 0:
            # else:
            # print(my_threads)
            # print(lines)
            # print("========")
            new = run_chunk(lines, my_threads, process_one_site)
            for lines in new:
                file_info.write(lines[0])
                file_info_detail.write(lines[1])
                if len(lines[2]) > 5:
                    # print(lines)
                    debug_num += 1
                    file_info_debug.write(lines[2])

            finished_ratio = line_num / total_number * 100
            logger.info(f"Finish {line_num} ,  {finished_ratio}%")
            lines = []
    else:
        new = run_chunk(lines, my_threads, process_one_site)
        for lines in new:
            file_info.write(lines[0])
            file_info_detail.write(lines[1])
            if len(lines[2]) > 5:
                # print(lines)
                debug_num += 1
                file_info_debug.write(lines[2])
        finished_ratio = line_num / total_number * 100
        logger.info(f"Finish {line_num} ,  {finished_ratio}%")

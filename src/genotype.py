#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: HapTR
# Script : genotype.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
import numpy as np

from src.units import *
from src.repeat import *
from src.region import *
from src.paras import args_init
import os
import pysam
import multiprocessing

os.environ['OPENBLAS_NUM_THREADS'] = '1'







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


def read_repeat_info(paras):
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
                repeat_infos_sorted[chrom].append(Region(chunk,threads=my_threads))
                chunk = []
            else:
                chunk.append(info[start])
        else:
            if len(chunk) > 0:
                repeat_infos_sorted[chrom].append(chunk)
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


def genotype(parase):
    if not args_init(parase):
        logger.error("Genotype init ERROR!")
        return -1
    paras = get_value("paras")
    repeat_infos = read_repeat_info(paras)

#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: HapTR
# Script : train.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: TODO
=============================================================================="""
import numpy as np
from src.init import args_init, get_reads_depth, extract_repeat_info
from src.units import *
from src.repeat import *
from src.run import Run
from src.region import *
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


def extract_train_sites():
    return


class Train(Run):
    def build_benchmarks(self):

        pass

    def run_one_region(self, region):
        return region

    # def run_one_batch(self, regions, pool):
    #
    #     res =
    #
    #
    #     return res

    def extract_reads_for_regions(self):
        batch_num = self.param.batch
        pool = multiprocessing.Pool(processes=self.param.threads)
        regions = []

        for chrom in self.repeats.keys():
            logger.info(f"Processing {chrom}")
            region_idx = 0
            for region_id in self.repeats[chrom]:
                region_idx += 1
                print(region_idx)
                if region_idx % batch_num == 0:
                    regions.append(self.repeats[chrom][region_id])
                    regions_res = pool.map(self.run_one_region, regions)
                    for region in regions_res:
                        self.repeats[chrom][region.region_id] = region
                    regions = []
                    del regions_res
                else:
                    regions.append(self.repeats[chrom][region_id])
            else:
                if len(regions) > 0:
                    regions_res = pool.map(self.run_one_region, regions)
                    for region in regions_res:
                        self.repeats[chrom][region.region_id] = region
                    del regions_res, regions

        pool.join()
        pool.close()

    def run(self):
        # if not args_init(self.paras):
        #     logger.error("Genotype init ERROR!")
        #     return -1
        # paras = get_value("paras")
        self.extract_repeat_info()
        # print(self.repeats)
        self.extract_reads_for_regions()

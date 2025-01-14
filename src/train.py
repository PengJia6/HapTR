#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""==============================================================================
# Project: HapTR
# Script : train.py
# Author : Peng Jia
# Date   : 2020.07.13
# Email  : pengjia@stu.xjtu.edu.cn
# Description: training model
=============================================================================="""
import pickle

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
    thread = thread if thread < 3 else (thread * 0.9)
    pool = multiprocessing.Pool(processes=int(thread))
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
    return [repeat.get_output_len_info(), "", ""]


def process_one_site_only_dis(line):
    line, paras = line
    path_bam = paras["input"]
    line_info = line[:-1].split("\t")
    chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content, repeat_type, \
        repeat_subtype, source, site_id, complex_repeat, annotation, up, down = line_info
    repeat = Repeat(chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content,
                    repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down, depth=paras["depth_dict"])
    repeat.process_reads(path_bam)
    return [repeat.get_output_len_info(), repeat.get_output_details(), repeat.get_info_debug()]


def write_output():
    pass


def extract_train_sites():
    return


def run_extract_feature_from_reads(region):
    if region.param.vcf4hap:
        region.extract_variant_locis()
    region.extract_feature_for_deep_train()
    return region


def write_features(regions, file):
    # for region_id,region in regions.items():
    #     pickle.dump({region_id,region},file)
    pickle.dump(regions, file)

    return None


class Train(Run):
    def build_benchmarks(self):
        pass

    def extract_feature_from_reads(self):
        batch_num = self.param.batch
        regions = []
        print(batch_num, "batch_num")
        print(self.param.threads, "param.threads")
        feature_file = open(f'{self.param.output_info}', "wb")
        for chrom in self.repeats.keys():
            logger.info(f"Processing {chrom}")
            region_idx = 0
            for region_id in self.repeats[chrom]:
                region_idx += 1
                # print(region_idx,region_id,)
                if region_idx % batch_num == 0:
                    regions.append(self.repeats[chrom][region_id])
                    pool = multiprocessing.Pool(processes=self.param.threads)
                    regions_res = pool.map(run_extract_feature_from_reads, regions)
                    pool.close()
                    pool.join()
                    for region in regions_res:
                        write_features(region, feature_file)

                    regions = []
                    del regions_res
                else:
                    regions.append(self.repeats[chrom][region_id])
            else:
                if len(regions) > 0:
                    pool = multiprocessing.Pool(processes=self.param.threads)
                    regions_res = pool.map(run_extract_feature_from_reads, regions)
                    pool.close()
                    pool.join()
                    for region in regions_res:
                        write_features(region, feature_file)
                    regions = []
                    del regions_res
        feature_file.close()
        # feature_file = open(f'{self.param.output_info}', "wb")
        with open(f'{self.param.output_info}', "rb") as file:
            num = 1
            while True:
                try:
                    num += 1
                    data = pickle.load(file)  # 逐个读取
                    print("读取的数据:", num)
                except EOFError:
                    break
        print("read scccess")

    def run(self):
        self.extract_repeat_info()
        self.extract_feature_from_reads()

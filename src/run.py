#!/usr/bin/env python 
"""
Project: HapTR
Script: Run.py
Author: Peng Jia
E-mail: pengjia@stu.xjtu.edu.cn
Time : 2023/10/27
Description: TODO
"""
from src.units import *
from src.repeat import Repeat
from src.region import Region
import time
import multiprocessing


class Run():

    def __init__(self, param):
        self.param = param
        pass

    def _extract_repeat(self, repeat_line):
        chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content, repeat_type, \
            repeat_subtype, source, site_id, complex_repeat, annotation, up, down = repeat_line[1][:-1].split(
            "\t")[:17]
        repeat = Repeat(chrom, int(start), int(end), strand, int(repeat_len), motif, motif_len,
                        average_repeat_times,
                        content,
                        repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down)
        return repeat

    def extract_repeat_info(self):
        logger.info(f"Exacting repeat from {self.param.repeat_bed} ...")
        repeat_infos = {}
        repeat_infos_sorted = {}
        repeat_info_num = {}
        for line in open(self.param.repeat_bed):
            chrom, start, end = line[:-1].split("\t")[:3]
            start = int(start)
            if chrom not in repeat_infos:
                repeat_infos[chrom] = {}
            repeat_infos[chrom][start] = line

        # pool = multiprocessing.Pool(processes=int(self.param.threads))
        my_batch = self.param.batch
        for chrom, info in repeat_infos.items():

            start_sorted_info = sorted(info.items(), key=lambda item: item[0])

            # this_chrom_repeat  = pool.map(self._extract_repeat, start_sorted_info)

            this_chrom_repeat = {}
            chunk = []
            for idx, repeat_line in enumerate(start_sorted_info, 1):
                repeat = self._extract_repeat(repeat_line)
                chunk.append(repeat)
                if idx % my_batch == 0:
                    region = Region(chunk)
                    this_chrom_repeat[region.region_id] = region
                    chunk = []
            else:
                if len(chunk) > 0:
                    region = Region(chunk)
                    this_chrom_repeat[region.region_id] = region
            repeat_info_num[chrom] = idx
            repeat_infos_sorted[chrom] = this_chrom_repeat
            logger.info(f"{chrom}: {idx} repeats.")
        self.repeats = repeat_infos_sorted
        total_repeat = sum([i for j, i in repeat_info_num.items()])
        self.total_repeat = total_repeat
        self.chrom_repeat = repeat_info_num
        self.repeats = repeat_infos_sorted
        logger.info(f'Total: {total_repeat} repeats.')

    def extract_repeat_reads(self):
        self.repeats_hap = []
        self.repeats_unhap = []

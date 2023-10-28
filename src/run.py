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

        my_threads = int(self.param.threads)
        my_batch = int(self.param.batch)

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


        for chrom, info in repeat_infos.items():

            start_sorted_info = sorted(info.items(), key=lambda item: item[0])

            # this_chrom_repeat  = pool.map(self._extract_repeat, start_sorted_info)

            this_chrom_repeat=[]
            for repeat_line in start_sorted_info:
                repeat=self._extract_repeat(repeat_line)
                this_chrom_repeat.append(repeat)


            repeat_num=len(this_chrom_repeat)
            repeat_info_num[chrom]=repeat_num
            repeat_infos_sorted[chrom] = this_chrom_repeat
            logger.info(f"{chrom}: {repeat_num} repeats.")
        # pool.close()
        # pool.join()
        self.repeats=repeat_infos_sorted

        total_repeat = sum([i for j, i in repeat_info_num.items()])
        self.total_repeat = total_repeat
        self.chrom_repeat = repeat_info_num
        self.repeats = repeat_infos_sorted
        logger.info(f'Total: {total_repeat} repeats.')

    def extract_repeat_reads(self):
        self.repeats_hap = []
        self.repeats_unhap = []

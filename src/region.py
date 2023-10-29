#!/usr/bin/env python 
"""
Project: HapTR
Script: regions.py
Author: Peng Jia
E-mail: pengjia@stu.xjtu.edu.cn
Time : 2023/10/26
Description: TODO
"""
import multiprocessing
from src.repeat import *


class Region:
    """
    Description: class Region
    Member variables:

        # paras: parameters of this program, containing command input and some default values
        # ms_list: information of microsatellite in this window
        # bam_path: bam file path
        region_start: start position of this window
        region_end: end position of this window
        reads: Read Object of this window; dict: read_id:Read
        reads_num: reads numbers
        # microsatellites: Microsatellite Object of this window; dict: ms_id:Microsatellite
        # microsatellites_id: list of ms_id
        # vcf_recs: list of vcf records infomation in this window
    Member methods:

    """

    def __init__(self, repeat_list, threads=1):
        self.threads = threads
        self.chrom = repeat_list[0].chrom
        self.win_start = repeat_list[0].start
        self.win_end = repeat_list[-1].end
        self.repeats = {f"{it.chrom}_{it.start}_{it.end}": it for it in repeat_list}
        self.repeats_id = [f"{it.chrom}_{it.start}_{it.end}" for it in repeat_list]
        self.region_id = f"{self.chrom}_{self.win_start}_{self.win_end}"

    def decode_repeats(self, repeat_list):
        pool = multiprocessing.Pool(processes=int(self.threads))
        # print("fun",fun)
        # print("trs",trs)
        res = pool.map(self.decode_one_repeat, repeat_list)
        pool.close()
        pool.join()
        return res

    def decode_one_repeat(self, line):
        chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content, repeat_type, \
            repeat_subtype, source, site_id, complex_repeat, annotation, up, down = line[:-1].split("\t")[:17]
        repeat = Repeat(chrom, int(start), int(end), strand, int(repeat_len), motif, motif_len, average_repeat_times,
                        content,
                        repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down)
        return repeat

    def extract_reads(self):
        for repeat_id in self.repeats_id:
            print(repeat_id)
        self.phased_num = 222
        return


if __name__ == '__main__':
    pass

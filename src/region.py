#!/usr/bin/env python 
"""
Project: HapTR
Script: regions.py
Author: Peng Jia
E-mail: pengjia@stu.xjtu.edu.cn
Time : 2023/10/26
Description: TODO
"""
import logging
import multiprocessing

import pysam

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

    def __init__(self, repeat_list, param, threads=1):
        self.threads = threads
        self.chrom = repeat_list[0].chrom
        self.win_start = repeat_list[0].start
        self.win_end = repeat_list[-1].end
        self.repeats = {f"{it.chrom}_{it.start}_{it.end}": it for it in repeat_list}
        self.repeats_id = [f"{it.chrom}_{it.start}_{it.end}" for it in repeat_list]
        self.region_id = f"{self.chrom}_{self.win_start}_{self.win_end}"
        self.param = param

    def extract_feature_for_train(self):
        print(self.region_id)
        return self
        pass

    def extract_feature_for_predict(self):
        pass
    def extract_feature_for_train_predict(self):
        pass



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

    # def _process_one_repeat(self, repeat):
    #     repeat.get_repeat_info()
    #
    #     return repeat

    #
    # res = pool.map(self._process_one_repeat, )
    # pool.close()
    # pool.join()

    def decoded_phasing_reads_info(self):
        logging.info(f"Processing region {self.region_id}")

        reads_dict = {}
        repeats_reads_dict = {}
        pool = multiprocessing.Pool(processes=int(self.param.threads))

        # # for repeat_id, repeat in self.repeats.items():
        #     reads=_extract_reads([repeat,bam])
        reads_dict = pool.map(_extract_reads, [[repeat, self.param] for repeat_id, repeat in self.repeats.items()])

        # print(reads_dict)
        # print(len(reads_dict))
        self.feature=reads_dict
        pool.close()
        pool.join()
        # bam.close()
        self.phased = True
        return self


def _extract_reads(info):
    repeat, param = info
    bam = pysam.Samfile(f"{param.input_bam_path}") if param.input_bam_type in "bam" else \
        pysam.Samfile(f"{param.input_bam_path}", reference_filename=f"{param.reference_path}")
    reads = [read.cigarstring for read in bam.fetch(repeat.chrom, repeat.start, repeat.end)]
    return reads


if __name__ == '__main__':
    pass

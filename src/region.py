#!/usr/bin/env python 
"""
Project: HapTR
Script: regions.py
Author: Peng Jia
E-mail: pengjia@xjtu.edu.cn
Time : 2023/10/26
Description: Region processing
"""
import logging
import multiprocessing

import pysam

from src.repeat import *
from src.read import ReadForTrain


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
        pysam_reads = {}
        sam_file = pysam.AlignmentFile(self.param.input_bam_path, mode="rb",
                                       reference_filename=self.param.reference_path)
        flank_size = self.param.flank_size
        # pysam.AlignmentFile(self.bam_path, mode="rb", reference_filename=self.reference_path)
        # TODO optimize
        for repeat_id, repeat in self.repeats.items():
            for alignment in sam_file.fetch(repeat.chrom, repeat.start - flank_size, repeat.end + flank_size + 1):
                if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary or alignment.is_supplementary:
                    continue
                if alignment.reference_start > repeat.start - flank_size - 1 or alignment.reference_end < repeat.end + flank_size + 1:
                    continue
                if len(alignment.query_sequence) < 2:
                    continue
                read_id = alignment.query_name + "_" + str(alignment.reference_start)
                if read_id not in pysam_reads:
                    # read = ReadForTrain(read_id=read_id, alignment=alignment, )
                    # reads[read_id] = read
                    pysam_reads[read_id] = alignment

                repeat.add_read_id(read_id=read_id, hap=int(alignment.get_tag("HP")) if alignment.has_tag("HP") else 0)
        reads_kept = {}
        for repeat_id, repeat in self.repeats.items():
            if not repeat.pass4train(min_phased_reads=self.param.min_phased_reads,
                                     min_phased_ratio=self.param.min_phased_ratio,
                                     min_depth=self.param.depths_dict["iqr_min"],
                                     max_depth=self.param.depths_dict["iqr_max"]):
                continue
            # print(repeat.support_read_nubmer_hap)
            for read_id in repeat.support_reads[1]+repeat.support_reads[2]:
                if read_id not in reads_kept:
                    read = ReadForTrain(read_id=read_id, )
                    reads_kept[read_id] = read
                else:
                    read = reads_kept[read_id]
                read.add_repeat(repeat_id)
        # print(reads_kept)
        for read_id, read in reads_kept.items():
            # print(read_id)

            features = read.extract_features(alignment=pysam_reads[read_id])
            for repeat_id, feature in features.items():
                self.repeats[repeat_id].set_train_features(feature)
        # TODO finish
        # TODO 提取特征，提取真实值

        # self.reads = reads

        # self.reads2 = reads2
        # self.reads_num = len(self.reads)

    # def extract_feature_for_train(self):
    #     self.init_reads_for_train()
    #
    #     # print(self.region_id)
    #     return
    #     pass

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
        self.feature = reads_dict
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

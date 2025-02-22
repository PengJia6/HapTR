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
from operator import index

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from collections import Counter
import pysam

from src.repeat import *
from src.read import ReadForTrain
from src.units import logger


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
        self.repeats = {f"TR_{it.chrom}_{it.start}_{it.end}": it for it in repeat_list}
        self.repeats_id = [f"TR_{it.chrom}_{it.start}_{it.end}" for it in repeat_list]
        self.region_id = f"{self.chrom}_{self.win_start}_{self.win_end}"
        self.repeat_num = len(repeat_list)
        self.param = param
        self.variants_positions = None
        self.repeats_fails = []
        self.variants_info = {}

    def extract_variant_locis(self):
        vcf_file = pysam.VariantFile(self.param.vcf4hap)
        sample = vcf_file.header.samples[0]
        positions = []
        for record in vcf_file.fetch(self.chrom, max(self.win_start - self.param.size4hap, 0), self.win_end + self.param.size4hap):
            gt = record.samples[sample]["GT"]
            if None in gt or gt[0] == gt[1] or len(record.alts) > 1: continue
            if len(record.alts) > 1 or len(record.alts[0]) != len(record.ref): continue
            positions.append(record.pos)
        self.variants_positions = positions if len(positions) > 0 else None

    def extract_feature_for_deep_train(self):
        reads = {}
        sam_file = pysam.AlignmentFile(self.param.input_bam_path, mode="rb",
                                       reference_filename=self.param.reference_path)
        flank_size = self.param.flank_size
        # print(flank_size)
        logger.info(f"Extracting features for region {self.region_id}")
        # print(aaa)
        # num = 0
        for repeat_id, repeat in self.repeats.items():
            # num += 1
            # if num % 100 == 0:
            # logger.info(f"Extracting features for region {self.region_id}, {repeat_id}, {num}, {num/self.repeat_num}")
            for alignment in sam_file.fetch(repeat.chrom, repeat.start - flank_size, repeat.end + flank_size + 1, ):
                if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary or alignment.mapping_quality < self.param.minimum_mapping_quality:
                    continue
                if alignment.reference_start > repeat.start - flank_size - 1 or alignment.reference_end < repeat.end + flank_size + 1:
                    continue
                if len(alignment.query_sequence) < 2:
                    continue
                read_id = alignment.query_name + "_" + str(alignment.reference_start)
                if read_id not in reads:
                    reads[read_id] = ReadForTrain(read_id=read_id, alignment=alignment, variant_pos=self.variants_positions)
                self.repeats[repeat_id].add_read_id(read_id=read_id, )

        read_id_kept = set()
        sam_file.close()
        # print("000")
        for repeat_id, repeat in self.repeats.items():
            if not repeat.pass4process(min_depth=self.param.depths_dict["iqr_min"],
                                       max_depth=self.param.depths_dict["iqr_max"]):
                self.repeats_fails.append(repeat_id)
                continue

            for read_id in repeat.repeat_feature:
                if read_id not in read_id_kept:
                    read_id_kept.add(read_id)
                reads[read_id].add_repeat(repeat_id, repeat)
                # print(repeat_id, read_id,)
        # print()
        # print(len(self.repeats))

        for read_id, read in reads.items():
            if read_id not in read_id_kept: continue
            if read.extract_deep_features():
                for repeat_id, repeat_str in read.repeat_str_read.items():
                    self.repeats[repeat_id].repeat_feature[read_id].set_seq(repeat_str)
                for repeat_id, repeat_mut in read.repeat_mut_read.items():
                    self.repeats[repeat_id].repeat_feature[read_id].set_mut(repeat_mut)
                for repeat_id, repeat_qual in read.repeat_qual_read.items():
                    self.repeats[repeat_id].repeat_feature[read_id].set_qual(repeat_qual)
            self.variants_info[read_id] = read.variant_info

        # for repeat_id, repeat in self.repeats.items():
        #     print(repeat_id)
        return

    def phase_TR(self):
        diff_indexs = []

        for repeat_id, repeat in self.repeats.items():
            features = pd.DataFrame()
            for read_id in repeat.repeat_feature:
                if read_id in self.variants_info:
                    read_info_df = pd.DataFrame({pos: [mut.upper() if (mut is not None and len(mut) == 1) else None] for pos, mut in self.variants_info[read_id].items()}, index=[read_id])
                    features = pd.concat([features, read_info_df])

            features = features.dropna(axis=1, thresh=3)
            features = features.dropna(axis=0, thresh=3)
            phased_reads_num = len(features.index)
            if phased_reads_num < self.param.depths_dict["iqr_min"] or phased_reads_num > self.param.depths_dict["iqr_max"]:
                repeat.set_phased_info(phased_status=False, muts_info=None, read_list=None)
            else:
                features.fillna("N", inplace=True)
                # repeat.phased = True
                muts_info = ["".join([it for it in info]) for read_id, info in features.iterrows()]
                read_list = [read_id for read_id, info in features.iterrows()]
                repeat.set_phased_info(phased_status=True, muts_info=muts_info, read_list=read_list)
        total_sites = 0
        phased_sites = 0
        for repeat_id, repeat in self.repeats.items():
            total_sites += 1
            if repeat.phased_status:
                phased_sites += 1
        self.total_sites_num = total_sites
        self.phased_sites_num = phased_sites

    def extract_motif_features_for_train(self):
        for repeat_id, repeat in self.repeats.items():
            repeat.extract_motif_features()
            repeat.combine4deep()


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

# def decoded_phasing_reads_info(self):
#     logging.info(f"Processing region {self.region_id}")
#
#     reads_dict = {}
#     repeats_reads_dict = {}
#     pool = multiprocessing.Pool(processes=int(self.param.threads))
#
#     # # for repeat_id, repeat in self.repeats.items():
#     #     reads=_extract_reads([repeat,bam])
#     reads_dict = pool.map(_extract_reads, [[repeat, self.param] for repeat_id, repeat in self.repeats.items()])
#
#     # print(reads_dict)
#     # print(len(reads_dict))
#     self.feature = reads_dict
#     pool.close()
#     pool.join()
#     # bam.close()
#     self.phased = True
#     return


#
# def _extract_reads(info):
#     repeat, param = info
#     bam = pysam.Samfile(f"{param.input_bam_path}") if param.input_bam_type in "bam" else \
#         pysam.Samfile(f"{param.input_bam_path}", reference_filename=f"{param.reference_path}")
#     reads = [read.cigarstring for read in bam.fetch(repeat.chrom, repeat.start, repeat.end)]
#     return reads


if __name__ == '__main__':
    pass

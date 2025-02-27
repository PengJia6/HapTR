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
        self.total_sites_num =0
        self.phased_sites_num = 0
        self.reads_info={}

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
    def extract_reads_info(self):
        reads = {}
        sam_file = pysam.AlignmentFile(self.param.input_bam_path, mode="rb",
                                       reference_filename=self.param.reference_path)
        flank_size = self.param.flank_size
        # print(flank_size)
        # logger.info(f"Extracting features for region {self.region_id}")
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
                    read = ReadForTrain(read_id=read_id, alignment=alignment, variant_pos=self.variants_positions)
                    if read.extract_reads_str():
                        reads[read_id] = read
                        self.variants_info[read_id] = read.extract_variant_feature()
                self.repeats[repeat_id].add_read_id(read_id=read_id, )
        sam_file.close()
        self.reads_info=reads

    def extract_feature_for_deep_train(self):
        # print(self.repeat_num,self.phased_sites_num)
        read_id_kept = set()
        process_repeats=set()
        for repeat_id, repeat in self.repeats.items():
            if not repeat.phased_status: continue
            if not repeat.pass4process(min_depth=self.param.depths_dict["iqr_min"],
                                       max_depth=self.param.depths_dict["iqr_max"]):
                self.repeats_fails.append(repeat_id)
                continue
            process_repeats.add(repeat_id)
            for read_id in repeat.repeat_feature:
                if read_id not in read_id_kept:
                    read_id_kept.add(read_id)
                self.reads_info[read_id].add_repeat(repeat_id, repeat)
                # print(repeat_id, read_id,)
        # print()
        # print(len(self.repeats))

        for read_id, read in self.reads_info.items():
            if read_id not in read_id_kept: continue
            read.extract_repeat_feature()
            # print(read.repeat_qual_read)
            for repeat_id in process_repeats:
                if repeat_id not in read.repeat_str_read: continue
                repeat_str= read.repeat_str_read[repeat_id]
                self.repeats[repeat_id].repeat_feature[read_id].set_seq(repeat_str)
                repeat_mut=read.repeat_mut_read[repeat_id]
                self.repeats[repeat_id].repeat_feature[read_id].set_mut(repeat_mut)
                repeat_qual= read.repeat_qual_read[repeat_id]
                self.repeats[repeat_id].repeat_feature[read_id].set_qual(repeat_qual)
        self.reads_info={}

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
        self.phased_sites_num = 0
        for repeat_id, repeat in self.repeats.items():
            total_sites += 1
            if repeat.phased_status:
                self.phased_sites_num += 1
    def select_features_for_haplotyping(self,model):
        for repeat_id, repeat in self.repeats.items():
            pass
        pass

    def extract_motif_features_for_train(self):
        regions=[]
        # print(self.region_id)
        num=120
        for repeat_id, repeat in self.repeats.items():
            # if repeat.phased_status:
            repeat.extract_motif_features()
            # print("1")
            repeat.combine4deep()
            # print("2")
        #     if repeat.reads_features is None or len(repeat.reads_features) <= 5: continue
        #     num+=1
        #     if num>120
        #     regions.extend(repeat.reads_features)
        # print(len(regions))
            # region_masks.extend(repeat.mask)

                # repeat.padding()


# def extract_feature_for_predict(self):
#     pass
#
#
# def extract_feature_for_train_predict(self):
#     pass
#
#
# def decode_repeats(self, repeat_list):
#     pool = multiprocessing.Pool(processes=int(self.threads))
#     # print("fun",fun)
#     # print("trs",trs)
#     res = pool.map(self.decode_one_repeat, repeat_list)
#     pool.close()
#     pool.join()
#     return res
#
#
# def decode_one_repeat(self, line):
#     chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content, repeat_type, \
#         repeat_subtype, source, site_id, complex_repeat, annotation, up, down = line[:-1].split("\t")[:17]
#     repeat = Repeat(chrom, int(start), int(end), strand, int(repeat_len), motif, motif_len, average_repeat_times,
#                     content,
#                     repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down)
#     return repeat


if __name__ == '__main__':
    pass

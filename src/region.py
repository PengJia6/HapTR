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
from wsgiref.util import request_uri

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
        self.eval_num=0

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
        self.reads = {}
        sam_file = pysam.AlignmentFile(self.param.input_bam_path, mode="rb",
                                       reference_filename=self.param.reference_path)
        flank_size = self.param.flank_size
        for repeat_id, repeat in self.repeats.items():
            for alignment in sam_file.fetch(repeat.chrom, repeat.start - flank_size, repeat.end + flank_size + 1, ):
                if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary or alignment.mapping_quality < self.param.minimum_mapping_quality:
                    continue
                if alignment.reference_start > repeat.start - flank_size - 1 or alignment.reference_end < repeat.end + flank_size + 1:
                    continue
                if len(alignment.query_sequence) < 2:
                    continue
                read_id = alignment.query_name + "_" + str(alignment.reference_start)
                if read_id not in self.reads:
                    read = ReadForTrain(read_id=read_id, alignment=alignment, variant_pos=self.variants_positions)
                    read.extract_reads_str()
                    self.reads[read_id] = read
                    self.variants_info[read_id] = read.extract_variant_feature()
                self.reads[read_id].add_repeat(repeat_id, repeat)
                self.repeats[repeat_id].add_read_id(read_id=read_id,)

        sam_file.close()
        for read_id, read in self.reads.items():
            read.extract_repeat_feature()
            for repeat_id in read.support_repeat_ids:
                repeat_str= read.repeat_str_read[repeat_id]
                repeat_mut = read.repeat_mut_read[repeat_id]
                repeat_qual = read.repeat_qual_read[repeat_id]
                self.repeats[repeat_id].add_repeat_feature(read_id,repeat_str,repeat_mut,repeat_qual)

        del self.reads

    # def set_latent_features(self,latent_features):
    #     for repeat_id, repeat in self.repeats.items():
    #         repeat.set_latent_features(latent_features)


    # def phase_TR(self):
    #     diff_indexs = []
    #     for repeat_id, repeat in self.repeats.items():
    #         features = pd.DataFrame()
    #         for read_id in repeat.repeat_feature:
    #             if read_id in self.variants_info:
    #                 read_info_df = pd.DataFrame({pos: [mut.upper() if (mut is not None and len(mut) == 1) else None] for pos, mut in self.variants_info[read_id].items()}, index=[read_id])
    #                 features = pd.concat([features, read_info_df])
    #         features = features.dropna(axis=1, thresh=3)
    #         features = features.dropna(axis=0, thresh=3)
    #         phased_reads_num = len(features.index)
    #         if phased_reads_num < self.param.depths_dict["iqr_min"] or phased_reads_num > self.param.depths_dict["iqr_max"]:
    #             repeat.set_phased_info(phased_status=False, muts_info=None, read_list=None)
    #         else:
    #             features.fillna("N", inplace=True)
    #             # repeat.phased = True
    #             muts_info = ["".join([it for it in info]) for read_id, info in features.iterrows()]
    #             read_list = [read_id for read_id, info in features.iterrows()]
    #             repeat.set_phased_info(phased_status=True, muts_info=muts_info, read_list=read_list)
    #     total_sites = 0
    #     self.phased_sites_num = 0
    #     for repeat_id, repeat in self.repeats.items():
    #         total_sites += 1
    #         if repeat.phased_status:
    #             self.phased_sites_num += 1
    def select_features_for_haplotyping(self,model):
        importance_dict={}
        for repeat_id, repeat in self.repeats.items():
            importance=repeat.evaluate_vae_output(model)
            if len(importance)>0:
                self.eval_num+=1
            for i,j in Counter(importance).items():
                if i not in importance_dict.keys():
                    importance_dict[i]=j
                else:
                    importance_dict[i]+= j
        self.importance_dict=importance_dict

            # for i in importance:
            #     if i in  importance_dict:
            #         importance_dict[i] += 1


        # pass



    def extract_features(self):
        for repeat_id, repeat in self.repeats.items():
            repeat.check_depth(min_depth=self.param.depths_dict["iqr_min"],
                               max_depth=self.param.depths_dict["iqr_max"])
            repeat.phase_TR(variants_info=self.variants_info,
                            min_reads=self.param.depths_dict["iqr_min"],
                            max_reads=self.param.depths_dict["iqr_max"])
            repeat.extract_motif_features()



if __name__ == '__main__':
    pass

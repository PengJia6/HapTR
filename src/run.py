#!/usr/bin/env python 
"""
Project: HapTR
Script: Run.py
Author: Peng Jia
E-mail: pengjia@xjtu.edu.cn
Time : 2023/10/27
Description: run ...
"""

from src.units import *
from src.repeat import Repeat
from src.region import Region
import time
import multiprocessing
import torch
from tqdm import tqdm

class Run():

    def __init__(self, param):
        self.param = param

    def _extract_repeat(self, repeat_line):
        chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content, repeat_type, \
            repeat_subtype, source, site_id, complex_repeat, annotation, up, down = repeat_line[:-1].split("\t")[:17]
        # print(repeat_line)
        repeat = Repeat(chrom, int(start), int(end), strand, int(repeat_len), motif, motif_len,
                        average_repeat_times, content, repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down, flank_size=self.param.flank_size)
        return repeat
    def run_extract_feature_from_readsxxx(self,region):
        if region.param.vcf4hap:
            region.extract_variant_locis()
        region.extract_feature_for_deep_train()
        region.phase_TR()
        region.extract_motif_features_for_train()
        return region

    def extract_repeat_info(self):
        logger.info(f"Exacting repeat from {self.param.repeat_bed} ...")
        ref_fa = pysam.FastaFile(self.param.reference_path)
        repeat_infos = {}
        repeat_infos_sorted = {}
        repeat_info_num = {}
        line_num = 0
        for line in open(self.param.repeat_bed):
            line_num += 1
            if line_num == 1: continue
            chrom, start, end = line[:-1].split("\t")[:3]

            start = int(start)
            if chrom not in repeat_infos:
                repeat_infos[chrom] = {}
            repeat_infos[chrom][start] = line

            # pool = multiprocessing.Pool(processes=int(self.param.threads))
        region_size = self.param.region_size
        for chrom, info in tqdm (repeat_infos.items(),desc="Loading repeat..."):
            # tqdm(range(self.epochs), desc="Training model ...", ):

            start_sorted_info = sorted(info.items(), key=lambda item: item[0])

            # this_chrom_repeat  = pool.map(self._extract_repeat, start_sorted_info)

            this_chrom_repeat = {}
            chunk = []
            for idx, repeat_line in enumerate(start_sorted_info, 1):
                # print(repeat_line)
                repeat = self._extract_repeat(repeat_line[1])
                repeat.content = ref_fa.fetch(repeat.chrom, repeat.start-1, repeat.end)
                repeat.up = ref_fa.fetch(repeat.chrom, repeat.start - repeat.flank_size , repeat.start )
                repeat.down = ref_fa.fetch(repeat.chrom, repeat.end, repeat.end + repeat.flank_size)
                # print("----")
                chunk.append(repeat)
                if idx % region_size == 0:
                    region = Region(chunk, self.param)
                    this_chrom_repeat[region.region_id] = region
                    chunk = []
            else:
                if len(chunk) > 0:
                    region = Region(chunk, self.param)
                    this_chrom_repeat[region.region_id] = region
            repeat_info_num[chrom] = idx
            repeat_infos_sorted[chrom] = this_chrom_repeat
            logger.info(f"{chrom}: {idx} repeats.")
        total_repeat = sum([i for j, i in repeat_info_num.items()])
        self.total_repeat = total_repeat
        self.chrom_repeat = repeat_info_num
        self.repeats = repeat_infos_sorted
        logger.info(f'Total: {total_repeat} repeats.')
        # print(self.repeats)
    # @staticmethod

    def pad_sequence_features(self, region, model_length=200, feature_dim=2):
        region_features = []
        region_masks = []
        # if feature_dim is None:
        #     feature_dim = len(region.repeats[0].)
        for repeat_id, repeat in region.repeats.items():
            if len(repeat.train_features) <= 5: continue
            ref_len = len(repeat.train_features[0])
            # for rd in repeat.train_features:
            pad_stat = False if ref_len >= model_length else True
            model_length = 200
            pad_length = model_length - ref_len

            reads_features = [torch.tensor(i + [[0] * feature_dim] * pad_length) if pad_stat else torch.tensor(i[:model_length]) for i in repeat.train_features]
            # print(reads_features.shape,)
            mask = [torch.cat((torch.zeros(ref_len), torch.ones(pad_length))) if pad_stat else torch.zeros(model_length) for i in repeat.train_features]
            # print(Counter([i.shape[0] for i in reads_features]))
            region_features.extend(reads_features)
            region_masks.extend(mask)
        return region_features, region_masks


    # def extract_repeat_reads(self):
    #     self.repeats_hap = []
    #     self.repeats_unhap = []

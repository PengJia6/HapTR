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
import torch

# from src.init import args_init, get_reads_depth, extract_repeat_info
from src.units import *
from src.repeat import *
from src.run import Run
from src.region import *
from src.model.vae import *
import os
import pysam
import multiprocessing

os.environ['OPENBLAS_NUM_THREADS'] = '1'


# def run_chunk(trs, thread, fun):
#     thread = thread if thread < 3 else (thread * 0.9)
#     pool = multiprocessing.Pool(processes=int(thread))
#     res = pool.map(fun, trs)
#     pool.close()
#     pool.join()
#     return res
#
#
# def process_one_site(line):
#     line, paras = line
#     path_bam = paras["input"]
#     line_info = line[:-1].split("\t")
#     chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content, repeat_type, \
#         repeat_subtype, source, site_id, complex_repeat, annotation, up, down = line_info
#     repeat = Repeat(chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content,
#                     repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down,
#                     depth=paras["depth_dict"])
#     repeat.process_reads(path_bam)
#     return [repeat.get_output_len_info(), "", ""]


# def process_one_site_only_dis(line):
#     line, paras = line
#     path_bam = paras["input"]
#     line_info = line[:-1].split("\t")
#     chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content, repeat_type, \
#         repeat_subtype, source, site_id, complex_repeat, annotation, up, down = line_info
#     repeat = Repeat(chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content,
#                     repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down, depth=paras["depth_dict"])
#     repeat.process_reads(path_bam)
#     return [repeat.get_output_len_info(), repeat.get_output_details(), repeat.get_info_debug()]
#



# def run_extract_feature_from_reads(region):
#     if region.param.vcf4hap:
#         region.extract_variant_locis()
#     region.extract_feature_for_deep_train()
#     region.phase_TR()
#     region.extract_motif_features_for_train()
#
#     return region


# def write_features(regions, file):
#     # for region_id,region in regions.items():
#     #     pickle.dump({region_id,region},file)
#     pickle.dump(regions, file)
#
#     return None


class Train(Run):

    def __init__(self, param):
        super().__init__(param)

    def train_vae(self):
        np.random.seed(42)
        torch.manual_seed(42)
        input_dim = 2
        sequence_len=200

        # 生成模拟数据
        # sequences, masks = torch.tensor(,self.train_feature_masks)
        # self.train_features = []
        print("sequences", self.train_features.shape)
        print("mask", self.train_feature_masks.shape)
        # 创建VAE模型实例
        vae_model = VAEModel(sequence_len=sequence_len, input_dim=input_dim)

        # 创建训练器实例
        trainer = VAETrainer(model=vae_model, data=TensorDataset(self.train_features, self.train_feature_masks), epochs=50, batch_size=32)

        # 训练VAE模型
        trainer.train()
        torch.save(trainer.model,f"{self.param.output_model}")
        print(self.param.output)

        # 提取潜在特征
        # latent_features = trainer.get_latent_features(self.train_features)
        #
        # 打印潜在特征的形状
        # print("Latent Features Shape:", latent_features.shape)

        # 获取训练过程中每个epoch的损失
        # losses = trainer.get_losses()
        # print("Losses over epochs:", losses)

    def run(self):
        self.extract_repeat_info()
        self.extract_feature_from_reads()
        self.train_vae()

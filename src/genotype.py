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
from numpy.testing.print_coercion_tables import print_new_cast_table

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


class Genotype(Run):

    def __init__(self, param):
        super().__init__(param)

    def testvae(self):
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
        print(self.param.output_model)

        # 提取潜在特征
        # latent_features = trainer.get_latent_features(self.train_features)
        #
        # 打印潜在特征的形状
        # print("Latent Features Shape:", latent_features.shape)

        # 获取训练过程中每个epoch的损失
        # losses = trainer.get_losses()
        # print("Losses over epochs:", losses)
    def load_vae_model(self):
        # self.vae_model = torch.load(self.param.model)
        self.vae_model= torch.load(f'{self.param.model}', map_location=torch.device('cpu'))
        self.vae_model.eval()
    def test_vae(self):
        for item in self.train_features:
            latent_feature = self.vae_model.generate_latent_features(torch.stack([item]))
            # z = self.vae_model.reparameterize(mu, log_var)
            print(latent_feature)
        pass
    def run(self):
        self.load_vae_model()
        self.extract_repeat_info()
        self.extract_feature_from_reads()
        self.test_vae()

        # self.testvae()

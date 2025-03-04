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

from src.run import Run
from src.region import *
from src.model.vae_liner import *
import os
import multiprocessing
import yaml

os.environ['OPENBLAS_NUM_THREADS'] = '1'

def init_dist():
    os.environ['MASTER_ADDR'] = 'localhost'
    os.environ['MASTER_PORT'] = '12346'
    # os.environ['CUDA_VISIBLE_DEVICES'] = '0,1,2'
    rank = 0
    world_size = 1
    dist.init_process_group("nccl", rank=rank, world_size=world_size)
    torch.multiprocessing.set_start_method('spawn')  # 强制安全模式
    # multiprocessing.set_start_method("spawn")
    # os.environ['PYTORCH_NO_CUDA_MEMORY_CACHING'] = '1'  # 禁用CUDA内存缓存


class Train(Run):

    def __init__(self, param):
        super().__init__(param)
        init_dist()
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.train_features = None
        self.train_feature_masks = None
        self.model=None


    def train_vae(self):
        np.random.seed(42)
        torch.manual_seed(42)
        input_dim = 2
        sequence_len=200

        # 生成模拟数据
        # sequences, masks = torch.tensor(,self.train_feature_masks)
        # self.train_features = []ll
        print("sequences", self.train_features.shape)
        print("mask", self.train_feature_masks.shape)

        # 创建VAE模型实例
        vae_model = VAEModel2(sequence_len=sequence_len, input_dim=input_dim)
        trainer = VAETrainer(model=vae_model,
                             data=TensorDataset(self.train_features.to(self.device), self.train_feature_masks.to(self.device)),
                             device=self.device, epochs=50, batch_size=128)
        trainer.train()
        self.model=trainer.save_model(self.param.output)
        self.model.share_memory()
        logger.info(f"Trained model saved to {self.param.output}")
        del self.train_features, self.train_feature_masks
        torch.cuda.empty_cache()

        # print(self.param.output)

        # 提取潜在特征
        # latent_features = trainer.get_latent_features(self.train_features)
        #
        # 打印潜在特征的形状
        # print("Latent Features Shape:", latent_features.shape)

        # 获取训练过程中每个epoch的损失
        # losses = trainer.get_losses()
        # print("Losses over epochs:", losses)
    def run_extract_train_features(self,region):
        if region.param.vcf4hap:
            region.extract_variant_locis()
        else:
            logger.error(f"The variant file was not provided")
        region.extract_reads_info( )
        region.extract_features()

        return region
    def extract_train_from_reads(self):
        total_sites = 0
        phased_sites = 0
        batch_num = self.param.batch
        regions = []
        logger.info(f"batch_num: {batch_num}")
        print(f"param.threads: {self.param.threads}")
        print(f"param.threads: {self.param.region_size}")
        # feature_file = open(f'{self.param.output_info}', "wb")
        features = []
        masks = []
        for chrom in self.repeats.keys():
            # logger.info(f"Processing {chrom}")
            region_idx = 0
            for region_id in tqdm(self.repeats[chrom],desc=f"Processing {chrom}"):
                region_idx += 1
                # print(region_idx,region_id,)
                if region_idx % batch_num == 0:
                    regions.append(self.repeats[chrom][region_id])
                    self.repeats[chrom][region_id] = None
                    pool = multiprocessing.Pool(processes=self.param.threads)
                    regions_res = pool.map(self.run_extract_train_features, regions)
                    pool.close()
                    pool.join()
                    for region in regions_res:
                        # write_features(region, feature_file)
                        total_sites += region.total_sites_num
                        phased_sites += region.phased_sites_num
                        reads_features, reads_masks,_ = self.pad_sequence_features(region,mode="all")
                        features.extend(reads_features)
                        masks.extend(reads_masks)
                    for region in regions_res:
                        self.repeats[chrom][region.region_id] =region
                    regions = []
                    del regions_res
                else:
                    regions.append(self.repeats[chrom][region_id])
                    self.repeats[chrom][region_id] = None

            else:
                if len(regions) > 0:
                    pool = multiprocessing.Pool(processes=self.param.threads)
                    regions_res = pool.map(self.run_extract_train_features, regions)
                    pool.close()
                    pool.join()
                    for region in regions_res:
                        total_sites += region.total_sites_num
                        phased_sites += region.phased_sites_num
                        reads_features, reads_masks,_ = self.pad_sequence_features(region,mode="all")
                        features.extend(reads_features)
                        masks.extend(reads_masks)
                    for region in regions_res:
                        self.repeats[chrom][region.region_id] =region
                    regions = []
                    del regions_res
        if features is None or len(features) == 0:
            logger.error("No features reads detection!")
        else:
            self.train_features = torch.stack(features, dim=0)
            self.train_feature_masks = torch.stack(masks, dim=0)

    def run_select_features_for_haplotyping(self,region):
        region.select_features_for_haplotyping(self.model)
        return region
    #todo 处理返回的importance

    def set_latent_features(self,chrom,region_id,repeat_ids,latent_features):
        for repeat_id,latent_feature in zip(repeat_ids,latent_features):
            self.repeats[chrom][region_id].repeats[repeat_id].set_latent_features(latent_feature)
        pass


    def select_features_for_haplotyping2(self):
        features=[]
        # masks=[]

        for chrom in self.repeats.keys():
            # logger.info(f"Processing {chrom}")
            region_idx = 0
            for region_id,region in tqdm(self.repeats[chrom].items(), desc=f"Processing {chrom}"):
                region_idx += 1
                reads_features, reads_masks,repeat_ids = self.pad_sequence_features(region,mode="phased")
                # print(reads_features)
                if len(reads_features) > 0:
                    latent_feature = self.model.generate_latent_features(torch.stack(reads_features).to(self.device))
                    # region.select_features_for_haplotyping(latent_feature)
                    self.set_latent_features(chrom,region_id,repeat_ids,latent_feature)

    def select_features_for_haplotyping(self):

        batch_num = self.param.batch
        regions = []
        logger.info(f"batch_num: {batch_num}")
        print(f"param.threads: {self.param.threads}")
        print(f"param.threads: {self.param.region_size}")
        # feature_file = open(f'{self.param.output_info}', "wb")
        importance={}
        eval_num=0
        for chrom in self.repeats.keys():
            # logger.info(f"Processing {chrom}")
            region_idx = 0
            for region_id in tqdm(self.repeats[chrom], desc=f"Processing {chrom}"):
                region_idx += 1
                # print(region_idx,region_id,)
                if region_idx % batch_num == 0:
                    regions.append(self.repeats[chrom][region_id])
                    self.repeats[chrom][region_id] = None
                    # pool = multiprocessing.Pool(processes=self.param.threads)
                    pool = multiprocessing.Pool(processes=2)
                    regions_res = pool.map(self.run_select_features_for_haplotyping, regions)
                    pool.close()
                    pool.join()
                    for region in regions_res:
                        self.repeats[chrom][region.region_id] = region
                        eval_num += region.eval_num
                        for i,j in region.importance_dict.items():
                            if i not in importance:
                                importance[i] = j
                            else:
                                importance[i] += j
                    regions = []
                    del regions_res
                else:
                    regions.append(self.repeats[chrom][region_id])
                    self.repeats[chrom][region_id] = None

            else:
                if len(regions) > 0:
                    # pool = multiprocessing.Pool(processes=self.param.threads)
                    pool = multiprocessing.Pool(processes=2)
                    regions_res = pool.map(self.run_select_features_for_haplotyping, regions)
                    pool.close()
                    pool.join()

                    for region in regions_res:
                        self.repeats[chrom][region.region_id] = region
                        eval_num += region.eval_num
                    # for im in regions_res:
                            #     self.repeats[chrom][region.region_id] = region
                        for i, j in region.importance_dict.items():
                            if i not in importance:
                                importance[i] = j
                            else:
                                importance[i] += j
                    regions = []
                    del regions_res
        print(importance)
        importance2={}
        for i,j in importance.items():
            importance2[int(i)] = int(j)
        with open(f"{self.param.output}.importance.yaml","w",encoding="utf-8") as f:
            data={"importance":importance2,
                       "eval_num":eval_num}
            yaml.dump(data,f,allow_unicode=True)
        # pd.DataFrame(importance).to_csv()
        # print(importance)
        # print(eval_num)
        # if features is None or len(features) == 0:
        #     logger.error("No features reads detection!")
        # else:
        #     self.train_features = torch.stack(features, dim=0)
        #     self.train_feature_masks = torch.stack(masks, dim=0)

                    # print(latent_feature.shape)
                # masks.extend(reads_masks)
        # if features is None or len(features) == 0:
        #     logger.error("No features reads detection!")
        # else:
        #     self.train_features = torch.stack(features, dim=0)
        #     self.train_feature_masks = torch.stack(masks, dim=0)
        #     for item in self.train_features:
        #         latent_feature = self.vae_model.generate_latent_features(torch.stack([item]).to(self.device))
        #         # z = self.vae_model.reparameterize(mu, log_var)
        #         print(latent_feature.shape)
            # print("eval",self.train_features.shape)

    def haplotype(self):
        pass
    def decode_repeat(self):
        pass

    def run(self):
        self.extract_repeat_info()
        self.extract_train_from_reads()
        self.train_vae()
        self.select_features_for_haplotyping()
        self.haplotype()




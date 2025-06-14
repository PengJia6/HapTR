# !/usr/bin/env python
"""
Project: HapTR
Script: find_anchor.py
Author: Peng Jia
E-mail: pengjia@xjtu.edu.cn
Time : 2023/5/25
Description: Repeat preprocessing
"""
# TODO processing phased reads
# using features of content sequence
# decode the structure of the repeat
#
import pandas as pd
import numpy as np
import pysam
import torch
from Bio.Seq import Seq
from Bio import Align
from Bio import Seq
from pandas import read_feather

from sklearn import mixture
from collections import Counter
from scipy import stats

from src.model.vae import device
from src.read import ReadFeature
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics.pairwise import euclidean_distances

model = RandomForestClassifier()


class AnchorFinder():
    def __init__(self, align_mode="local", substitution_matrix="BLASTN"):
        names = Align.substitution_matrices.load()
        # print(names)
        matrix = Align.substitution_matrices.load(name=substitution_matrix)
        aligner = Align.PairwiseAligner(substitution_matrix=matrix, mode=align_mode, )
        score_open_gap = -7
        score_extend_gap = -2
        aligner.target_internal_open_gap_score = score_open_gap
        aligner.target_internal_extend_gap_score = score_extend_gap
        aligner.target_left_open_gap_score = score_open_gap
        aligner.target_left_extend_gap_score = score_extend_gap
        aligner.target_right_open_gap_score = score_open_gap
        aligner.target_right_extend_gap_score = score_extend_gap
        aligner.query_internal_open_gap_score = score_open_gap
        aligner.query_internal_extend_gap_score = score_extend_gap
        aligner.query_left_open_gap_score = score_open_gap
        aligner.query_left_extend_gap_score = score_extend_gap
        aligner.query_right_open_gap_score = score_open_gap
        aligner.query_right_extend_gap_score = score_extend_gap
        self.aligner = aligner
        self.substitution_matrix = matrix
        self.substitution_matrix_name = substitution_matrix
        self.align_mode = align_mode

    def _get_aligner(self):
        return self.aligner

    def _set_anchor(self, up, down):
        """
        :param up: upstream sequence of the repeat sites
        :param down: upstream sequence of the repeat sites
        :return:
        """
        ""

    def _find_anchor(self, read_str, up, down, aligned_ratio=0.95):

        alignments_up = self.aligner.align(read_str, up)
        if len(alignments_up) < 1:
            up_pos = [-1, -1]
        else:
            aligned = alignments_up[0].aligned
            sub_len = len(alignments_up[0].substitutions)
            if (alignments_up[0].substitutions * np.diag([1] * sub_len)).sum() < aligned_ratio * len(up):
                up_pos = [-1, -1]
            else:
                up_pos = [aligned[0][-1][-1], aligned[1][-1][-1]]
        alignments_down = self.aligner.align(read_str, down)
        if len(alignments_down) < 1:
            down_pos = [-1, -1]
        else:
            aligned = alignments_down[0].aligned
            sub_len = len(alignments_down[0].substitutions)
            if (alignments_down[0].substitutions * np.diag([1] * sub_len)).sum() < aligned_ratio * len(down):
                down_pos = [-1, -1]
            else:
                down_pos = [aligned[0][0][0], aligned[1][0][0]]
        if up_pos[0] > down_pos[0]:
            up_pos, down_pos = [[-1, -1]] * 2
        content = read_str[up_pos[0]:down_pos[0]]

        # print([up_pos, down_pos])
        # if down_pos[0]-up_pos[0]>100:
        #     # print(read_str)
        #     print("========================")
        #     print([i.aligned for i in alignments_down])
        #     print("-------")
        #     print([i.aligned for i in alignments_up])
        # for i in alignments_down:
        #     print(i
        #           )
        # # print([i.path for i in alignments_down])
        # print("----------")
        # for i in alignments_up:
        #     print(i
        #           )
        # print([ i.aligned for i in alignments_up]
        #       )
        return [up_pos, down_pos, content]


def _get_more_times(repeat_list):
    """
    Args:
        repeat_list (list): repeat length of microsatellite site, such as [10,10,9,9,9,9,8,5,5,5,5,5]
    Returns: central value of the input list

    """
    counts = Counter(repeat_list)
    f = sorted(zip(counts.values(), counts.keys()))
    # print(f)
    if len(f) == 1:
        return f[-1][-1]
    else:
        # print(abs(f[-1][0] - f[-2][0]) / len(repeat_list))
        if abs(f[-1][0] - f[-2][0]) / len(repeat_list) < 0.3:  # TODO add on command input
            # print("ldldldl")
            return int(round(np.mean(repeat_list)))
        else:
            return f[-1][-1]


class Repeat:
    def __init__(self, chrom, start, end, strand, repeat_len, motif, motif_len, average_repeat_times, content,
                 repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down, flank_size):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.motif = motif
        # self.up_reverse_complement = self.up.reverse_complement()
        # self.down_reverse_complement = self.down.reverse_complement()
        self.strand = strand
        self.repeat_len = int(repeat_len)
        self.motif_len = motif_len
        self.average_repeat_times = average_repeat_times
        self.content2 = content
        self.up = up
        self.down = down
        self.repeat_type = repeat_type
        self.repeat_subtype = repeat_subtype
        self.source = source
        self.site_id = f"{chrom}_{start}_{end}"
        self.complex_repeat = complex_repeat
        self.annotation = annotation
        self.flank_size = int(flank_size)
        self.repeat_str = {}
        self.repeat_qual = {}
        self.repeat_mut = {}
        self.support_reads_id = []
        self.repeat_feature = {}
        self.phased_status = False
        self.phased_status_hidden = False
        self.read_id2id = {}
        self.read_id2feature = {}
        self.read_str2id = {}
        self.support_read_number = 0
        self.hidden_features_low_dim = None
        self.available_reads = []
        self.hidden_features = None
        # self.train=False

        # self.reads_cover_complete = {}
        # self.reads_cover_part = {}
        # self.anchor_finder = AnchorFinder()
        # self.anchor_len_up = len(self.up)
        # self.anchor_len_down = len(self.down)
        #
        # self.num_read_cover_repeat_compete = 0
        # self.num_read_total = 0
        # self.read_content = {}
        # self.content_consistent = []
        # self.depth = 0
        # self.qual = 0
        # self.filters = []
        # self.gt = [None, None]
        # self.af = [0, 0]
        # self.read_infos = []
        # self.debug = {"hete": {}, "homo": {}}
        # self.cluster_num = 0
        # self.dis_raw = ""
        # self.dis = []  # sorted
        # self.dis_hap1 = []

        # self.dis_hap2 = []
        # self.dis_hap0 = []
        # #
        # self.dis_str = ""
        # self.dis_str_hap1 = ""
        # self.dis_str_hap2 = ""
        # self.dis_str_hap0 = ""

    @staticmethod
    def edit_distance(str1, str2):
        # 直接统计对应位置字符不同的数量
        totol_len = 0
        dist = 0
        for c1, c2 in zip(str1, str2):
            if c1 == "N" or c2 == "N":
                continue
            else:
                totol_len += 1
                if c1 != c2:
                    dist += 1
        if totol_len == 0:
            return 1
        return dist / totol_len

    @staticmethod
    def compute_distance_matrix(sequences):
        n = len(sequences)
        dist_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                dist = Repeat.edit_distance(sequences[i], sequences[j])  # 计算两个序列之间的Levenshtein距离
                dist_matrix[i][j] = dist
                dist_matrix[j][i] = dist  # 距离矩阵是对称的
        return dist_matrix

    def phase_TR(self, variants_info, min_reads=0, max_reads=1000, min_support_site=3):
        features = pd.DataFrame()
        for read_id in self.support_reads_id:
            if read_id in variants_info:
                read_info_df = pd.DataFrame({pos: [mut.upper() if (mut is not None and len(mut) == 1) else None] for pos, mut in variants_info[read_id].items()}, index=[read_id])
                features = pd.concat([features, read_info_df])
        features = features.dropna(axis=1, thresh=min_support_site)
        features = features.dropna(axis=0, thresh=min_support_site)
        phased_reads_num = len(features.index)
        if phased_reads_num < min_reads or phased_reads_num > max_reads:
            self.phased_reads_num = phased_reads_num
            self.phased_reads_info = {}
            self.phased_bias = 0
            self.phased_status = False
        else:
            features.fillna("N", inplace=True)
            muts_info = ["".join([it for it in info]) for read_id, info in features.iterrows()]
            read_list = [read_id for read_id, info in features.iterrows()]
            self.phased_reads_num = len(read_list)
            dist_matrix = self.compute_distance_matrix(muts_info)
            dist_matrix = pdist(dist_matrix, metric='euclidean')
            Z = linkage(dist_matrix, method='average')
            labels = fcluster(Z, t=2, criterion='maxclust')
            clusters_rations = {i: j / self.phased_reads_num for i, j in Counter(labels).items()}
            if 2 not in clusters_rations:
                clusters_rations[2] = 0
            diff_ratio = np.abs(clusters_rations[1] - clusters_rations[2])
            phased_reads = {}
            for read_id, label in zip(read_list, labels):
                phased_reads[read_id] = label
            self.phased_reads_info = phased_reads
            self.phased_bias = diff_ratio
            if diff_ratio > 0.2:
                self.phased_status = False
            else:
                self.phased_status = True

    def cluster_reads(self):
        if (self.hidden_features_low_dim is None) or len(self.hidden_features_low_dim) != len(self.available_reads) or (not self.train):
            self.phased_bias_hidden = None
            return

        features = {}
        for i, j in zip(self.available_reads, self.hidden_features_low_dim):
            features[i] = j
        new_features = []
        read_list = []

        for read_id in self.available_reads:
            core_id = self.read_id2id[read_id]
            new_features.append(features[core_id])
            read_list.append(read_id)
        # dist_matrix = self.compute_distance_matrix_hidden_features(new_features)
        dist_matrix = euclidean_distances(new_features)
        print(new_features,dist_matrix)
        self.hidden_features = new_features

        dist_matrix = pdist(dist_matrix, metric='euclidean')
        Z = linkage(dist_matrix, method='average')
        labels = fcluster(Z, t=2, criterion='maxclust')
        print(Counter(labels))
        self.available_reads_num = len(self.available_reads)

        clusters_rations = {i: j / self.available_reads_num for i, j in Counter(labels).items()}
        if 2 not in clusters_rations:
            clusters_rations[2] = 0
        diff_ratio = np.abs(clusters_rations[1] - clusters_rations[2])
        phased_reads = {}
        for read_id, label in zip(read_list, labels):
            phased_reads[read_id] = label
        self.phased_reads_info_hidden = phased_reads
        self.phased_bias_hidden = diff_ratio
        if diff_ratio > 0.2:
            self.phased_status_hidden = False
        else:
            self.phased_status_hidden = True

    def merge_reads(self):
        if self.phased_status:
            cluster_reads = {}
            for read_id, label in self.phased_reads_info.items():
                if read_id in self.read_id2id:
                    seq_list = self.read_id2feature[self.read_id2id[read_id]].seq_list
                    if label not in cluster_reads:
                        cluster_reads[label] = [seq_list]
                    else:
                        cluster_reads[label].append(seq_list)
            hap_info = {}
            for label, info in cluster_reads.items():
                hap_info[label] = []
                for pos in [list(row) for row in zip(*info)]:
                    hap_info[label].append(Counter(pos).most_common(1)[0][0])
        elif self.phased_status_hidden:
            cluster_reads = {}
            for read_id, label in self.phased_reads_info_hidden.items():
                if read_id in self.read_id2id:
                    seq_list = self.read_id2feature[self.read_id2id[read_id]].seq_list
                    if label not in cluster_reads:
                        cluster_reads[label] = [seq_list]
                    else:
                        cluster_reads[label].append(seq_list)
            hap_info = {}
            for label, info in cluster_reads.items():
                hap_info[label] = []
                for pos in [list(row) for row in zip(*info)]:
                    hap_info[label].append(Counter(pos).most_common(1)[0][0])
        else:
            # TODO  hidden_features_为空
            # print(self.hidden_features_low_dim,self.hidden_features, self.phased_bias, self.phased_bias_hidden, self.normal_depth, self.available_reads)
            pass

            # print(Counter(i))
        # print(hap_info)
        pass

    def decode_structure(self):
        pass
        # if self.phased_status_hidden and self.phased_status:
        #     var_info=""
        #     hidden_info=""
        #     for i,j in self.phased_reads_info_hidden.items():
        #         if i in self.phased_reads_info:
        #             var_info+=f"{self.phased_reads_info[i]}"
        #             hidden_info+=f"{j}"
        #     print(var_info)
        #     print(hidden_info)
        #     print("=================")

    def add_read_id(self, read_id):
        # self.repeat_feature[read_id] = ReadFeature()
        self.support_reads_id.append(read_id)
        self.support_read_number += 1

    def add_repeat_feature(self, read_id, repeat_str, repeat_mut, repeat_qual):
        if repeat_str is None:
            return
        else:
            read_str = "".join(repeat_str)
            if read_str not in self.read_str2id:
                self.read_str2id[read_str] = read_id
                self.read_id2id[read_id] = read_id
                read_feature = ReadFeature()
                read_feature.set_seq(repeat_str)
                read_feature.set_mut(repeat_mut)
                read_feature.set_qual(repeat_qual)
                self.read_id2feature[read_id] = read_feature
            else:
                self.read_id2id[read_id] = self.read_str2id[read_str]

    def extract_motif_features(self, k=3):
        ref_str = f"{self.up}{self.content}{self.down}"
        ref_kmers = [ref_str[i:i + k] for i in range(len(ref_str) - k + 1)]
        ref_kmers_counter = dict(Counter(ref_kmers))
        ref_kmers_list = [ref_kmers_counter[i] for i in ref_kmers]
        available_reads_features = []
        available_reads = []

        for read_id, repeat_features in self.read_id2feature.items():
            if repeat_features.seq_list is None or len(repeat_features.seq) < 1: continue
            read_kmer3 = dict(Counter([repeat_features.seq[i:i + k] for i in range(len(repeat_features.seq) - k + 1)]))
            read2ref_kmer3 = [read_kmer3[i] / j if i in read_kmer3 else 0 for i, j in zip(ref_kmers, ref_kmers_list)]
            repeat_features.set_kmers(read2ref_kmer3)
            if repeat_features.mut_list is None or repeat_features.kmers is None: continue
            features = [[i, j] for i, j in zip(repeat_features.mut_list, repeat_features.kmers + [0])]
            available_reads_features.append(features)
            available_reads.append(read_id)
        self.train_features = available_reads_features
        self.available_reads = available_reads

    def features_scoring(self, X, y):
        model = RandomForestClassifier()
        model.fit(X, y)
        importance = model.feature_importances_
        return importance

    def set_latent_features(self, latent_features):
        # print(latent_features.shape,len(self.available_reads),len(self.phased_reads_info))
        self.hidden_features = latent_features

        pass

    def evaluate_vae_output(self, model, model_length=200, feature_dim=2, top=20):
        with torch.no_grad():
            if self.phased_status and self.normal_depth:
                ref_len = len(self.train_features[0])
                pad_stat = False if ref_len >= model_length else True
                model_length = 200
                pad_length = model_length - ref_len
                reads_features = [torch.tensor(i + [[0] * feature_dim] * pad_length) if pad_stat else torch.tensor(i[:model_length]) for i in self.train_features]
                # print(reads_features.shape,)
                # mask = [torch.cat((torch.zeros(ref_len), torch.ones(pad_length))) if pad_stat else torch.zeros(model_length) for i in repeat.train_features]
                reads_features = torch.stack(reads_features, dim=0)
                # mask=torch.stack(mask, dim=0)
                latent_features = model.generate_latent_features(reads_features.to(device))
                # print(latent_features.shape,len(self.available_reads),)
                # reads=np.sum([len(self.read_id2id[i]) for i in self.available_reads])
                # print(reads,len(self.phased_reads_info))

                torch.cuda.empty_cache()

                features = {}
                for read_id, feature in zip(self.available_reads, latent_features.to("cpu")):
                    features[read_id] = feature.numpy()
                X = []
                y = []
                for read_id, label in self.phased_reads_info.items():
                    if read_id in self.read_id2id:
                        if self.read_id2id[read_id] in features:
                            X.append(features[self.read_id2id[read_id]])
                            y.append(label)
                if len(X) < 2: return []
                importance = self.features_scoring(X, y)
                if np.sum(importance) > 0:
                    top20_indices = np.argsort(importance)[-top:][::-1]
                    return top20_indices
                else:
                    return []
            return []

    @property
    def train(self):
        return False if (not self.normal_depth) or len(self.train_features) < 2 else True

    # def padding(self,model_length=200,feature_dim=2):
    #     if  self.train_features is not None or  len(self.train_features) <= 5:
    #         self.reads_features=None
    #         self.mask=None
    #     else:
    #         ref_len = len(self.train_features[0])
    #         # for rd in repeat.train_features:
    #         pad_stat = False if ref_len >= model_length else True
    #         model_length = 200
    #         pad_length = model_length - ref_len
    #         self.reads_features = [torch.tensor(i + [[0] * feature_dim] * pad_length) if pad_stat else torch.tensor(i[:model_length]) for i in self.train_features]
    #         self.mask = [torch.cat((torch.zeros(ref_len), torch.ones(pad_length))) if pad_stat else torch.zeros(model_length) for i in self.train_features]

    def check_depth(self, min_depth=10, max_depth=10000):
        if self.support_read_number < min_depth or self.support_read_number > max_depth:
            self.normal_depth = False
        else:
            self.normal_depth = True

    def _k2_cluster(self, dis, iter_max=2, min_shift=1, cent0_ab=1, cent1_ab=1):
        dis = {int(i): (j) for i, j in dis.items()}
        sorted_value = sorted(dis.items(), key=lambda x: x[1])
        sorted_key = sorted(dis.items(), key=lambda x: x[0])
        dis_array = sorted([i for i, j in dis.items() for _ in range(j)])
        if len(dis_array) == 1:
            return [sorted_key[0]] * 2
        else:
            cent0, cent1 = sorted([sorted_value[-1][0], sorted_value[-2][0]])
            for it in range(iter_max):
                print(it)
                cent0_dis, cent1_dis = [], []
                tot_distance = 0
                for k, v in sorted_key:
                    disntance_c0 = abs(k - cent0) * cent0_ab
                    disntance_c1 = abs(k - cent1) * cent1_ab
                    tot_distance = disntance_c1 + disntance_c0
                    if disntance_c0 <= disntance_c1:
                        cent0_dis.extend([k] * v)
                    else:
                        cent1_dis.extend([k] * v)
                cent0_new = cent0_dis[len(cent0_dis) // 2]
                cent1_new = cent1_dis[len(cent1_dis) // 2]
                if abs(cent0_new - cent0) + abs(cent1_new - cent1) < min_shift:
                    break
                else:
                    cent0, cent1 = cent0_new, cent1_new
            return cent0, cent1

    def _get_read_pos(self, cigar_tuple: list, query_pos, ref_start, direction="up"):
        # pos_bias=0
        ref_pos = ref_start
        ref_end = ref_pos
        read_pos = 0
        read_end = 0
        read_query = None

        for cigar in cigar_tuple:
            if cigar[0] in [0, 7, 8]:
                ref_end = cigar[1] + ref_pos
                read_end = cigar[1] + read_pos
            elif cigar[0] in [1, 4]:
                read_end = cigar[1] + read_pos
            elif cigar[0] in [2]:
                ref_end = cigar[1] + ref_pos
            else:
                return None
            if ref_pos <= query_pos < ref_end:
                if cigar[0] not in [2]:
                    read_query = read_pos + (query_pos - ref_pos)
                elif cigar[0] not in [2]:
                    read_query = read_pos
                return read_query
            read_pos = read_end
            ref_pos = ref_end
        return read_query

    def _get_repeat_info(self, param):
        for read in pysam.Samfile(path_bam).fetch(self.chrom, self.start, self.end):
            if read.is_secondary or read.mapping_quality < 1 or len(
                    read.query_sequence) < 1:  # second alignment or low mapping quality
                continue  # TODO set mapping quality as a parameter
        #     if read.reference_end <= anchor_up or read.reference_start >= anchor_down:  # no overlap with this regions
        #         continue
        #     read_up, read_down, ref_str = None, None, ""
        #     if read.reference_start <= anchor_up - self.anchor_len_up and read.reference_end >= anchor_down + self.anchor_len_down:
        #         read_up = self.get_read_pos(cigar_tuple=read.cigartuples, query_pos=anchor_up - self.anchor_len_up,
        #                                     ref_start=read.reference_start)
        #         read_down = self.get_read_pos(cigar_tuple=read.cigartuples,
        #                                       query_pos=anchor_down + self.anchor_len_down,
        #                                       ref_start=read.reference_start)
        #         read_str = read.query_sequence[read_up:read_down]
        #
        #     elif read.reference_start <= anchor_up - self.anchor_len_up:
        #         this_cigar = read.cigartuples[0]
        #         if this_cigar[0] in [4, 2] and this_cigar[1] > self.anchor_len_up:
        #             read_up = self.get_read_pos(cigar_tuple=read.cigartuples, query_pos=anchor_up - self.anchor_len_up,
        #                                         ref_start=read.reference_start)
        #             # read_down = self.get_read_pos(cigar_tuple=read.cigartuples,
        #             #                               query_pos=anchor_down + self.anchor_len_down,
        #             #                               ref_start=read.reference_start)
        #             read_str = read.query_sequence[read_up:]
        #         else:
        #             continue
        #     elif read.reference_end >= anchor_down + self.anchor_len_down:
        #         this_cigar = read.cigartuples[-1]
        #
        #         if this_cigar[-1] in [4, 2] and this_cigar[1] > self.anchor_len_down:
        #             # read_up = self.get_read_pos(cigar_tuple=read.cigartuples, query_pos=anchor_up - self.anchor_len_up,
        #             #                             ref_start=read.reference_start)
        #
        #             read_down = self.get_read_pos(cigar_tuple=read.cigartuples,
        #                                           query_pos=anchor_down + self.anchor_len_down,
        #                                           ref_start=read.reference_start)
        #             read_str = read.query_sequence[:read_down]
        #         else:
        #             continue
        #     else:
        #         continue
        #     if len(read_str) < 5: continue
        #     self.num_read_total += 1
        #     up_info, down_info, content = self.anchor_finder.find_anchor(read_str=read_str, up=self.up,
        #                                                                  down=self.down)
        #     # if up_info[0] - down_info[0] < -100:
        #     #     print(read.reference_start, read.reference_end, read.cigarstring
        #     #           )
        #     #     print("-----")
        #
        #     if -1 in up_info or -1 in down_info: continue  # TODO add more parameters for anchor finding
        #     up_bais = self.anchor_len_up - up_info[1]
        #     down_bais = down_info[1]
        #     this_len = down_info[0] - up_info[0] - down_bais - up_bais
        #     self.read_content[content] = read.query_name
        #     self.reads_cover_complete[read.query_name] = [this_len, up_bais, down_bais]
        #     self.num_read_cover_repeat_compete += 1
        #     read_names.append(read.query_name)
        #     read_features.append([this_len - self.repeat_len, up_bais, down_bais])
        #     if this_len in reads_len_dis:
        #         reads_len_dis[this_len].append(read.query_name)
        #     else:
        #         reads_len_dis[this_len] = [read.query_name]
        #     contents.append(content)
        #     if read.has_tag("HP"):
        #         hap = f'hap{read.get_tag("HP")}'
        #     else:
        #         hap = "hap0"
        #     if this_len in reads_len_dis_hap[hap]:
        #         reads_len_dis_hap[hap][this_len].append(read.query_name)
        #     else:
        #         reads_len_dis_hap[hap][this_len] = [read.query_name]
        # read_features = np.array(read_features)
        # self.read_infos = [read_names, read_features, contents]
        # if self.num_read_cover_repeat_compete < self.depth["iqr_min"] or self.num_read_cover_repeat_compete > \
        #         self.depth["iqr_max"]:
        #     self.filters.append("Abnormal_read_depth")
        # elif self.num_read_cover_repeat_compete / self.num_read_total < 0.5:
        #     self.filters.append("Low_quality_reads")
        #
        # # print(read_features)
        #
        # # read_contents = list(set([j for i, j in self.read_content.items()]))
        # def sort_dis(dis: dict):
        #     sorted_key = sorted(dis.items(), key=lambda x: x[0])
        #     return [(i, dis[i]) for i, j in sorted_key]
        #
        # self.dis_raw = sort_dis({i: len(j) for i, j in reads_len_dis.items()})
        # self.dis_hap0 = sort_dis({i: len(j) for i, j in reads_len_dis_hap["hap0"].items()})
        # self.dis_hap1 = sort_dis({i: len(j) for i, j in reads_len_dis_hap["hap1"].items()})
        # self.dis_hap2 = sort_dis({i: len(j) for i, j in reads_len_dis_hap["hap2"].items()})
        # dis_str = ";".join([f"{i}:{j}" for i, j in self.dis_raw])
        # hap0_dis_str = ";".join([f"{i}:{j}" for i, j in self.dis_hap0])
        # hap1_dis_str = ";".join([f"{i}:{j}" for i, j in self.dis_hap1])
        # hap2_dis_str = ";".join([f"{i}:{j}" for i, j in self.dis_hap2])
        # self.dis_str = dis_str
        return

    def _process_reads(self, path_bam):
        # num_total_read = 0
        # num_cover_repeat_compete = 0
        # num_cover_repeat_part = 0
        # if self.site_id != "chr1_44835_44876":
        #     return
        anchor_up = self.start - self.anchor_len_up
        anchor_down = self.end + self.anchor_len_down
        read_names, read_features = [], []
        reads_len_dis = {}
        reads_len_dis_hap = {"hap0": {}, "hap1": {}, "hap2": {}}
        contents = []

        for read in pysam.Samfile(path_bam).fetch(self.chrom, self.start, self.end):
            if read.is_secondary or read.mapping_quality < 1 or len(
                    read.query_sequence) < 1:  # second alignment or low mapping quality
                continue  # TODO set mapping quality as a parameter
            if read.reference_end <= anchor_up or read.reference_start >= anchor_down:  # no overlap with this regions
                continue
            read_up, read_down, ref_str = None, None, ""
            if read.reference_start <= anchor_up - self.anchor_len_up and read.reference_end >= anchor_down + self.anchor_len_down:
                read_up = self.get_read_pos(cigar_tuple=read.cigartuples, query_pos=anchor_up - self.anchor_len_up,
                                            ref_start=read.reference_start)
                read_down = self.get_read_pos(cigar_tuple=read.cigartuples,
                                              query_pos=anchor_down + self.anchor_len_down,
                                              ref_start=read.reference_start)
                read_str = read.query_sequence[read_up:read_down]

            elif read.reference_start <= anchor_up - self.anchor_len_up:
                this_cigar = read.cigartuples[0]
                if this_cigar[0] in [4, 2] and this_cigar[1] > self.anchor_len_up:
                    read_up = self.get_read_pos(cigar_tuple=read.cigartuples, query_pos=anchor_up - self.anchor_len_up,
                                                ref_start=read.reference_start)
                    # read_down = self.get_read_pos(cigar_tuple=read.cigartuples,
                    #                               query_pos=anchor_down + self.anchor_len_down,
                    #                               ref_start=read.reference_start)
                    read_str = read.query_sequence[read_up:]
                else:
                    continue
            elif read.reference_end >= anchor_down + self.anchor_len_down:
                this_cigar = read.cigartuples[-1]

                if this_cigar[-1] in [4, 2] and this_cigar[1] > self.anchor_len_down:
                    # read_up = self.get_read_pos(cigar_tuple=read.cigartuples, query_pos=anchor_up - self.anchor_len_up,
                    #                             ref_start=read.reference_start)

                    read_down = self.get_read_pos(cigar_tuple=read.cigartuples,
                                                  query_pos=anchor_down + self.anchor_len_down,
                                                  ref_start=read.reference_start)
                    read_str = read.query_sequence[:read_down]
                else:
                    continue
            else:
                continue
            if len(read_str) < 5: continue
            self.num_read_total += 1
            up_info, down_info, content = self.anchor_finder.find_anchor(read_str=read_str, up=self.up,
                                                                         down=self.down)
            # if up_info[0] - down_info[0] < -100:
            #     print(read.reference_start, read.reference_end, read.cigarstring
            #           )
            #     print("-----")

            if -1 in up_info or -1 in down_info: continue  # TODO add more parameters for anchor finding
            up_bais = self.anchor_len_up - up_info[1]
            down_bais = down_info[1]
            this_len = down_info[0] - up_info[0] - down_bais - up_bais
            self.read_content[content] = read.query_name
            self.reads_cover_complete[read.query_name] = [this_len, up_bais, down_bais]
            self.num_read_cover_repeat_compete += 1
            read_names.append(read.query_name)
            read_features.append([this_len - self.repeat_len, up_bais, down_bais])
            if this_len in reads_len_dis:
                reads_len_dis[this_len].append(read.query_name)
            else:
                reads_len_dis[this_len] = [read.query_name]
            contents.append(content)
            if read.has_tag("HP"):
                hap = f'hap{read.get_tag("HP")}'
            else:
                hap = "hap0"
            if this_len in reads_len_dis_hap[hap]:
                reads_len_dis_hap[hap][this_len].append(read.query_name)
            else:
                reads_len_dis_hap[hap][this_len] = [read.query_name]
        read_features = np.array(read_features)
        self.read_infos = [read_names, read_features, contents]
        if self.num_read_cover_repeat_compete < self.depth["iqr_min"] or self.num_read_cover_repeat_compete > \
                self.depth["iqr_max"]:
            self.filters.append("Abnormal_read_depth")
        elif self.num_read_cover_repeat_compete / self.num_read_total < 0.5:
            self.filters.append("Low_quality_reads")

        # print(read_features)

        # read_contents = list(set([j for i, j in self.read_content.items()]))
        def sort_dis(dis: dict):
            sorted_key = sorted(dis.items(), key=lambda x: x[0])
            return [(i, dis[i]) for i, j in sorted_key]

        self.dis_raw = sort_dis({i: len(j) for i, j in reads_len_dis.items()})
        self.dis_hap0 = sort_dis({i: len(j) for i, j in reads_len_dis_hap["hap0"].items()})
        self.dis_hap1 = sort_dis({i: len(j) for i, j in reads_len_dis_hap["hap1"].items()})
        self.dis_hap2 = sort_dis({i: len(j) for i, j in reads_len_dis_hap["hap2"].items()})
        dis_str = ";".join([f"{i}:{j}" for i, j in self.dis_raw])
        hap0_dis_str = ";".join([f"{i}:{j}" for i, j in self.dis_hap0])
        hap1_dis_str = ";".join([f"{i}:{j}" for i, j in self.dis_hap1])
        hap2_dis_str = ";".join([f"{i}:{j}" for i, j in self.dis_hap2])
        self.dis_str = dis_str
        self.hap0_dis_str = hap0_dis_str
        self.hap1_dis_str = hap1_dis_str
        self.hap2_dis_str = hap2_dis_str

        # cluster_num = len(reads_len_dis)
        # self.cluster_num = cluster_num
        # if cluster_num == 0:
        #     self.gt = [None, None]
        #     self.qual = 0
        #     return
        # elif cluster_num == 1:
        #     # print(reads_len_dis)
        #     self.gt = [list(reads_len_dis.keys())[0]] * 2
        #     self.qual = 0
        #     # self.content_consistent = read_contents[0]
        #     return
        # else:
        #     hete_model = dpgmm = mixture.GaussianMixture(n_components=2,
        #                                                  covariance_type='full',
        #                                                  # tol=0.0001,
        #                                                  # n_init=2,
        #                                                  # reg_covar=1e-5,
        #                                                  # init_params="kmeans",
        #                                                  max_iter=10000, )
        #     homo_model = mixture.GaussianMixture(n_components=1,
        #                                          covariance_type='full',
        #                                          # tol=0.0001,
        #                                          # n_init=2,
        #                                          # reg_covar=1e-5,
        #                                          # init_params="kmeans",
        #                                          max_iter=10000, )
        #     hete_model.fit(X=read_features)
        #     homo_model.fit(X=read_features)
        #     # print(hete_model.get_params())
        #     # print(homo_model.get_params())
        #     self.debug["hete"]["weight"] = f"{hete_model.weights_[0]},{hete_model.weights_[1]}"
        #     self.debug["hete"]["bic"] = hete_model.bic(read_features)
        #     self.debug["hete"]["aic"] = hete_model.aic(read_features)
        #     self.debug["hete"]["score"] = hete_model.score(read_features)
        #     self.debug["hete"]["mean"] = f"{hete_model.means_[0][0]},{hete_model.means_[1][0]}"
        #     self.debug["hete"]["var"] = f"{hete_model.covariances_[0][0][0]},{hete_model.covariances_[1][0][0]}"
        #
        #     self.debug["homo"]["weight"] = homo_model.weights_[0]
        #     self.debug["homo"]["bic"] = homo_model.bic(read_features)
        #     self.debug["homo"]["aic"] = homo_model.aic(read_features)
        #     self.debug["homo"]["score"] = homo_model.score(read_features)
        #     self.debug["homo"]["mean"] = f"{homo_model.means_[0][0]}"
        #     self.debug["homo"]["var"] = f"{homo_model.covariances_[0][0][0]}"
        #
        #     # print(self.debug)
        #
        #     # self.debug[""]=""
        #
        #     pre_dis = {}
        #     pre_num = {}
        #     # print(dpgmm.means_)
        #     read_predict = []
        #     for read_name, read_feature in zip(read_names, read_features):
        #         k_pre = dpgmm.predict(np.array([read_feature]))
        #         read_predict.append(k_pre)
        #         if k_pre[0] not in pre_dis:
        #             pre_dis[k_pre[0]] = []
        #             pre_num[k_pre[0]] = 0
        #         pre_dis[k_pre[0]].append(read_feature[0])
        #         pre_num[k_pre[0]] += 1
        #     m = sorted(pre_num.keys(), key=(lambda x: pre_num[x]))
        #     if pre_num[m[-1]] > self.num_read_cover_repeat_compete * 0.7:
        #         # genotype = [int(round(np.mean(pre_dis[m[-1]])))] * 2
        #         genotype = [get_more_times(pre_dis[m[-1]])] * 2
        #
        #         qual = self.num_read_cover_repeat_compete / (1 + np.std(pre_dis[m[-1]]))
        #         self.af = [(0.5 * pre_num[m[-1]]) / self.num_read_cover_repeat_compete] * 2
        #     else:
        #         # hap1 = int(round(np.mean(pre_dis[m[-1]])))
        #         hap1 = get_more_times(pre_dis[m[-1]])
        #         qual1 = pre_num[m[-1]] / (1 + np.std(pre_dis[m[-1]]))
        #         # hap2 = int(round(np.mean(pre_dis[m[-2]])))
        #         hap2 = get_more_times(pre_dis[m[-2]])
        #         qual2 = pre_num[m[-2]] / (1 + np.std(pre_dis[m[-2]]))
        #         # if get_more_times(1)
        #
        #         hap1_dis = pre_dis[m[-1]]
        #         hap2_dis = pre_dis[m[-2]]
        #         # print("ks:",stats.kstest(hap1_dis,hap2_dis))
        #         # print("t:",stats.kstest(hap1_dis,hap2_dis))
        #
        #         genotype = [hap1, hap2]
        #         qual = (qual1 + qual2) / 2
        #         self.af = [pre_num[m[-1]] / self.num_read_cover_repeat_compete,
        #                    pre_num[m[-2]] / self.num_read_cover_repeat_compete]
        #     self.gt = genotype
        #     self.qual = qual
        # print(self.gt, self.qual, self.filters)
        # print(pre_dis)
        # print(pre_num)
        # print("----------------------------------")

    def _get_output_len_info(self):
        # if None in self.gt:
        #     return ""
        filters = ",".join(self.filters)
        # af_str=f"AF:{self.af[0]},{self.af[1]}"
        out_str = f"{self.chrom}\t{self.start}\t{self.end}\t{self.repeat_len}\tGT:{self.gt[0]},{self.gt[1]}\t" \
                  f"{self.qual}\t{self.num_read_cover_repeat_compete}\t{self.num_read_total}\t{self.dis_str}\t" \
                  f"{self.hap0_dis_str}\t{self.hap1_dis_str}\t{self.hap2_dis_str}\t{filters}\n"
        return out_str

    def _get_output_details(self):
        if None in self.gt:
            return ""
        filters = ",".join(self.filters)
        # af_str=f"AF:{self.af[0]},{self.af[1]}"
        out_str = f">{self.chrom}\t{self.start}\t{self.end}\t{self.repeat_len}\tGT:{self.gt[0]},{self.gt[1]}\t" \
                  f"{self.qual}\t{self.num_read_cover_repeat_compete}\t{self.num_read_total}\t{self.dis_str}\t" \
                  f"{self.hap0_dis_str}\t{self.hap1_dis_str}\t{self.hap2_dis_str}\t{filters}\n"
        out_str += f"{self.content}\t{self.repeat_len}\treference\n"
        a1, a2, a3 = self.read_infos
        for read_name, read_f, read_s in zip(a1, a2, a3):
            # names = ",".join(name)
            rd_f = ",".join([str(i) for i in read_f])
            out_str += f"{read_s}\t{rd_f}\t{read_name}\n"
        return out_str

    def _get_info_debug(self):
        filters = ",".join(self.filters)
        # af_str=f"AF:{self.af[0]},{self.af[1]}"
        if len(self.filters) > 0 or self.cluster_num < 2:
            return ""
        out_str = f"{self.chrom}\t{self.start}\t{self.end}\t{self.repeat_len}\tGT:{self.gt[0]},{self.gt[1]}\t" \
                  f"{self.qual}\t{self.num_read_cover_repeat_compete}\t{self.num_read_total}\t{self.dis_str}\t" \
                  f"{self.hap0_dis_str}\t{self.hap1_dis_str}\t{self.hap2_dis_str}\t{filters}\t" \
                  f"{self.debug['hete']['weight']}\t{self.debug['hete']['bic']}\t{self.debug['hete']['aic']}\t{self.debug['hete']['score']}\t{self.debug['hete']['mean']}\t{self.debug['hete']['var']}\t" \
                  f"{self.debug['homo']['weight']}\t{self.debug['homo']['bic']}\t{self.debug['homo']['aic']}\t{self.debug['homo']['score']}\t{self.debug['homo']['mean']}\t{self.debug['homo']['var']}" \
                  f"\n"
        return out_str

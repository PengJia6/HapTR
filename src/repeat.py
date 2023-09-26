# !/usr/bin/env python
import numpy as np
import pysam

"""
Project: PHD_programs
Script: find_anchor.py
Author: Peng Jia
E-mail: pengjia@stu.xjtu.edu.cn
Time : 2023/5/25
Description: TODO
"""
# TODO processing phased reads
# using features of content sequence
# decode the structure of the repeat
#


from Bio.Seq import Seq
from Bio import Align
from Bio import Seq
from sklearn import mixture
from collections import Counter
from scipy import stats

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

    def get_aligner(self):
        return self.aligner

    def set_anchor(self, up, down):
        """
        :param up: upstream sequence of the repeat sites
        :param down: upstream sequence of the repeat sites
        :return:
        """
        ""

    def find_anchor(self, read_str, up, down, aligned_ratio=0.95):

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


def get_more_times(repeat_list):
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
                 repeat_type, repeat_subtype, source, complex_repeat, annotation, up, down, depth):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.up = up
        self.down = down
        self.motif = motif
        # self.up_reverse_complement = self.up.reverse_complement()
        # self.down_reverse_complement = self.down.reverse_complement()
        self.strand = strand
        self.repeat_len = int(repeat_len)
        self.motif_len = motif_len
        self.average_repeat_times = average_repeat_times
        self.content = content
        self.repeat_type = repeat_type
        self.repeat_subtype = repeat_subtype
        self.source = source
        self.site_id = f"{chrom}_{start}_{end}"
        self.complex_repeat = complex_repeat
        self.annotation = annotation
        self.reads_cover_complete = {}
        self.reads_cover_part = {}
        self.anchor_finder = AnchorFinder()
        self.anchor_len_up = len(self.up)
        self.anchor_len_down = len(self.down)

        self.num_read_cover_repeat_compete = 0
        self.num_read_total = 0
        self.read_content = {}
        self.content_consistent = []
        self.depth = depth
        self.qual = 0
        self.filters = []
        self.gt = [None, None]
        self.af = [0, 0]
        self.read_infos = []
        self.debug = {"hete": {}, "homo": {}}
        self.cluster_num = 0
        self.dis_raw = ""
        self.dis = []  # sorted
        self.dis_hap1 = []
        self.dis_hap2 = []
        self.dis_hap0 = []
        #
        self.dis_str = ""
        self.dis_str_hap1 = ""
        self.dis_str_hap2 = ""
        self.dis_str_hap0 = ""

    def k2_cluster(self, dis, iter_max=2, min_shift=1, cent0_ab=1, cent1_ab=1):
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

    def get_read_pos(self, cigar_tuple: list, query_pos, ref_start,direction="up"):
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
                return  None
            if ref_pos <= query_pos < ref_end:
                if cigar[0] not in [2]:
                    read_query = read_pos + (query_pos - ref_pos)
                elif cigar[0] not in [2]:
                    read_query = read_pos
                return read_query
            read_pos = read_end
            ref_pos = ref_end
        return read_query

    def process_reads(self, path_bam):
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
            if read.is_secondary or read.mapping_quality < 1 or  len(read.query_sequence)<1:  # second alignment or low mapping quality
                continue  # TODO set mapping quality as a parameter
            if read.reference_end <= anchor_up or read.reference_start >= anchor_down:  # no overlap with this regions
                continue
            read_up,read_down,ref_str=None,None,""
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
            if len(read_str) < 5:continue
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
            return [(i, dis[i]) for i,j in sorted_key]

        self.dis_raw = sort_dis({i:len(j) for i, j in reads_len_dis.items()})
        self.dis_hap0 = sort_dis({i:len(j) for i, j in reads_len_dis_hap["hap0"].items()})
        self.dis_hap1 = sort_dis({i:len(j) for i, j in reads_len_dis_hap["hap1"].items()})
        self.dis_hap2 = sort_dis({i:len(j) for i, j in reads_len_dis_hap["hap2"].items()})
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

    def get_output_len_info(self):
        # if None in self.gt:
        #     return ""
        filters = ",".join(self.filters)
        # af_str=f"AF:{self.af[0]},{self.af[1]}"
        out_str = f"{self.chrom}\t{self.start}\t{self.end}\t{self.repeat_len}\tGT:{self.gt[0]},{self.gt[1]}\t" \
                  f"{self.qual}\t{self.num_read_cover_repeat_compete}\t{self.num_read_total}\t{self.dis_str}\t" \
                  f"{self.hap0_dis_str}\t{self.hap1_dis_str}\t{self.hap2_dis_str}\t{filters}\n"
        return out_str

    def get_output_details(self):
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

    def get_info_debug(self):
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

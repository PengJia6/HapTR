#!/usr/bin/env python 
"""
Project: HapTR
Script: Read.py
Author: Peng Jia
E-mail: pengjia@xjtu.edu.cn
Time : 2023/11/9
Description: Alignment reads processing
"""
import logging

import pysam

from src.units import logger


class Read:
    """
    Description: Read
    """

    # support_repeats = []

    def __init__(self, read_id, alignment: pysam.AlignedSegment):
        self.support_repeats = {}
        self.support_repeat_num = 0
        self.read_id = read_id
        self.alignment = alignment
        self.has_qual = True if self.alignment.query_qualities is not None else False
        self.variant_info = {}
    def add_repeat(self, repeat_id, repeat):
        self.support_repeats[repeat_id] = repeat
        self.support_repeat_num += 1

    # def set_varaint(self, varaint_positions: list):
    #     self.varaints = varaint_positions

    # def extract_features(self):
    #     # print(alignment)
    #     # self.microsatellites_num = len(self.microsatellites)
    #     # print(self.read_id, self.microsatellites_num,len(self.support_microsatellites))
    #     read_muts = {}
    #     repeat_lengths = {}
    #     for ms_id, ms in self.microsatellites.items():
    #         ms_start = ms.start
    #         ms_end = ms.end
    #         ms_start_pre = ms.start_pre
    #         ms_end_suf = ms.end_suf
    #         query_repeat_length = len(
    #             "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
    #         prefix = self.this_read_list[ms_start_pre - self.align_start:ms_start - self.align_start]
    #         suffix = self.this_read_list[ms_end - self.align_start:ms_end_suf - self.align_start]
    #         ms_content = self.this_read_list[ms_start - self.align_start:ms_end - self.align_start]
    #         mismatches = []
    #         deletions = []
    #         insertions = []
    #         pos_based_info = {}
    #         pos_deletion = []
    #         ref_pos = ms_start_pre - 2
    #         for pot in range(ms_start_pre - self.align_start - 1, ms_end_suf - self.align_start):
    #             ref_pos += 1
    #             this_read_base = self.this_read_list[pot]
    #             this_ref_base = self.this_ref_list[pot]
    #             if this_read_base == this_ref_base:
    #                 continue
    #             else:
    #                 if ref_pos < ms_start_pre:
    #                     band = [1]
    #                 elif ref_pos < ms_start:
    #                     band = [2]
    #                 elif ref_pos < ms_end:
    #                     band = [3]
    #                 elif ref_pos < ms_end_suf:
    #                     band = [4]
    #                 else:
    #                     band = [5]
    #                 this_read_base_len = len(this_read_base)
    #                 this_ref_base_len = len(this_ref_base)
    #                 if this_read_base_len == this_ref_base_len:
    #                     mismatches.append([ref_pos, this_read_base, band])
    #                     pos_based_info[ref_pos] = ["SNV", this_ref_base, this_read_base, band]
    #                 else:
    #                     if this_read_base_len < this_ref_base_len:
    #                         deletions.append([ref_pos, this_read_base, band])
    #                         # pos_based_info[ref_pos] = ["DEL", this_ref_base, ""]
    #                         pos_deletion.append([ref_pos, this_ref_base, band[0]])
    #                     else:
    #                         if ref_pos == ms_start - 1 and this_read_base_len > ms.repeat_unit_len:
    #                             band = [3]
    #                         insertions.append([ref_pos, this_read_base, band])
    #                         if band[0] != 3:
    #                             pos_based_info[ref_pos] = ["INS", "", this_read_base[1:], band]
    #
    #         if len(pos_deletion) > 0:
    #             deletion_start = pos_deletion[0][0]
    #             deletion_len = 1
    #             del_str = pos_deletion[0][1]
    #             band = [pos_deletion[0][2]]
    #             for i in range(1, len(pos_deletion)):
    #                 if pos_deletion[i - 1][0] + 1 == pos_deletion[i][0]:
    #                     deletion_len += 1
    #                     del_str += pos_deletion[i][1]
    #                     band.append(pos_deletion[i][2])
    #                 else:
    #                     pos_based_info[deletion_start] = ["DEL", del_str, "", set(band)]
    #
    #                     deletion_len = 1
    #                     deletion_start = pos_deletion[i][0]
    #                     del_str = pos_deletion[i][1]
    #                     band = [pos_deletion[i][2]]
    #             # print(del_str,deletion_start,band)
    #             if len(set(band)) > 1 or list(set(band))[0] != 3:
    #                 pos_based_info[deletion_start] = ["DEL", del_str, "", set(band)]
    #         read_str = self.this_read_list[ms_start_pre - 1 - self.align_start:ms_end_suf + 2 - self.align_start]
    #
    #         read_muts[ms_id] = Read_Mutation(self.read_id, repeat_length=query_repeat_length, strand=self.strand,
    #                                          hap=self.hap,
    #                                          mismatches=mismatches, deletions=deletions, insertions=insertions,
    #                                          prefix=prefix, suffix=suffix, ms_content=ms_content,
    #                                          pos_based_info=pos_based_info,
    #                                          read_str=read_str
    #                                          )
    #         repeat_lengths[ms_id] = query_repeat_length
    #     self.repeat_lengths = repeat_lengths
    #     self.mut_info = read_muts
    #     return {}


class ReadFeature:
    #

    def __init__(self, ):
        self.mut_list = None
        self.has_qual = None
        self.seq_list = None
        self.qual_list = None
        self.kmers = None
        # pass

    def set_seq(self, seq_list):
        self.seq_list = seq_list

    def set_mut(self, mut_list):
        self.mut_list = mut_list

    def set_qual(self, qual_list):
        if qual_list:
            self.qual_list = qual_list
            self.has_qual = True
        else:
            self.has_qual = False

    def set_kmers(self, kmers):
        self.kmers = kmers


class ReadForTrain(Read):
    def __init__(self, read_id, alignment, variant_pos, ):
        super().__init__(read_id, alignment)
        if variant_pos:
            self.variant_pos = sorted([i for i in variant_pos if (self.alignment.reference_start < i < self.alignment.reference_end)])
            self.support_variant_num = len(variant_pos)
        else:
            self.support_variant_num = 0

        # print("------")
        # self.start = alignment.reference_start

        # print(self.start)

        # self.phased = False if self.hap == 0 else True
        #
        # self.chrom = chrom
        # self.read_name = alignment.query_name
        # self.align_start = alignment.reference_start
        # self.align_end = alignment.reference_end
        # self.this_read_str = alignment.query_sequence.upper()
        # self.tech = tech
        # # print(alignment)
        # if tech == "contig":
        #     self.this_read_quals = []
        # else:
        #     self.this_read_quals = "".join([chr(i + 33) for i in alignment.query_qualities])
        # self.strand = False if alignment.is_reverse else True  # True: forward False: reverse
        # self.this_read_list = []
        # self.this_quals_list = []
        # self.this_ref_str = ""
        # self.this_ref_list = []
        # self.read_id = read_id
        # self.reference = reference
        # self.support_microsatellites = []
        # # self.alignment = alignment
        # self.microsatellites = {}
        # self.direction = False if alignment.is_reverse else True
        # self.hap = int(alignment.get_tag("HP")) if alignment.has_tag("HP") else 0
        # # print(self.hap)
        # self.cigartuples = alignment.cigartuples
        # self.mut_info = {}
        # self.repeat_lengths = {}

    # def get_microsatellite_detail(self, ms_info):
    #     self. = ms_info

        # return self.variant_info
        # print(self.read_id, self.support_variant_num, self.support_repeat_num)
        # repeat_info = {}
        # variant_info = {}
        #
        # for repeat_id, repeat in self.support_repeats.items():
        #     print(repeat_id)

        # return

    def extract_reads_str(self):
        sub_read_str = []
        sub_read_quals = []
        sub_read_insertion_deletion = []
        # sub_read_deletion = []

        read_pos = 0
        this_read_str = self.alignment.query_sequence
        if self.has_qual:
            this_read_qual = "".join([chr(i + 33) for i in self.alignment.query_qualities])

        for cigartuple in self.alignment.cigartuples:
            if cigartuple[0] in [0, 7, 8]:
                match_read = list(this_read_str[read_pos:read_pos + cigartuple[1]])
                sub_read_str.extend(match_read)
                if self.has_qual:
                    match_quals = list(this_read_qual[read_pos:read_pos + cigartuple[1]])
                    sub_read_quals.extend(match_quals)
                sub_read_insertion_deletion.extend([0] * cigartuple[1])
                read_pos += cigartuple[1]
            elif cigartuple[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
                if cigartuple[0] == 1:
                    if len(sub_read_str) == 0:
                        continue
                    else:
                        sub_read_str[-1] += this_read_str[read_pos:read_pos + cigartuple[1]]
                        # sub_read_insertion[]
                        if self.has_qual:
                            sub_read_quals[-1] += this_read_qual[read_pos:read_pos + cigartuple[1]]
                        sub_read_insertion_deletion[-1] = cigartuple[1]
                elif cigartuple[0] == 5:
                    continue
                read_pos += cigartuple[1]
            elif cigartuple[0] in [2, 3]:  # 2:D; 3:N: skip region of reference
                sub_read_str.extend([""] * cigartuple[1])
                if self.has_qual:
                    sub_read_quals.extend([""] * cigartuple[1])
                sub_read_insertion_deletion.extend([-1] * cigartuple[1])
            else:
                print("MMMMM==============")
                return None
        self.read_str = sub_read_str
        self.read_quals = sub_read_quals
        self.read_muts = sub_read_insertion_deletion
        return True

        # self.readsub_read_str, sub_read_quals

    def extract_variant_feature(self):
        varaints_info = {}
        # range_pos = [self.variant_pos[0], self.variant_pos[-1]]
        pass_var_num = 0
        for pos in self.variant_pos:
            # print(pos, range_pos)
            if pos < self.variant_pos[0] or pos > self.variant_pos[-1]:
                varaints_info[pos] = "N"
                logger.warn("The variant position is out of read aligned range.")
            else:
                varaints_info[pos] = self.read_str[pos - self.alignment.reference_start - 1]
                pass_var_num += 1
        # pass_var_num = pass_var_num
        return varaints_info

    def extract_repeat_feature(self):
        repeat_str_read = {}
        repeat_quals_read = {}
        repeat_mut_read = {}
        for repeat_id, repeat in self.support_repeats.items():

            repeat_str_read[repeat_id] = self.read_str[repeat.start - self.alignment.reference_start - repeat.flank_size:repeat.end - self.alignment.reference_start + repeat.flank_size]
            repeat_mut_read[repeat_id] = self.read_muts[repeat.start - self.alignment.reference_start - repeat.flank_size:repeat.end - self.alignment.reference_start + repeat.flank_size]
            if self.has_qual:
                repeat_quals_read[repeat_id] = self.read_quals[repeat.start - self.alignment.reference_start - repeat.flank_size:repeat.end - self.alignment.reference_start + repeat.flank_size]
        self.repeat_str_read = repeat_str_read
        self.repeat_qual_read = repeat_quals_read
        self.repeat_mut_read = repeat_mut_read

    #######################################################################################
    def extract_true_hete(self):
        # TODO finish
        return {}

    def get_read_str(self):
        pass
        # TODO:  To be optimized

        self.this_ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, start=self.align_start,
                                                                  end=self.align_end).upper()
        self.this_ref_list = list(self.this_ref_str)
        sub_read_str = []
        sub_read_quals = []
        # read_mut_info = ReadInfo()
        # read_mut_info.direction = self.direction
        # read_mut_info.hap = self.hap
        read_pos = 0
        for cigartuple in self.cigartuples:
            if cigartuple[0] in [0, 7, 8]:  # 0 : M : match or mishmatch ; 7: :=:match; 8:X:mismatch
                match_read = list(self.this_read_str[read_pos:read_pos + cigartuple[1]])
                match_quals = list(self.this_read_quals[read_pos:read_pos + cigartuple[1]])
                sub_read_str.extend(match_read)
                sub_read_quals.extend(match_quals)
                read_pos += cigartuple[1]
            elif cigartuple[0] in [1, 4, 5]:  # 1:I:inserion ;4:S:soft clip 5:H:hardclip
                if cigartuple[0] == 1:
                    if len(sub_read_str) == 0:
                        continue
                    else:
                        sub_read_str[-1] += self.this_read_str[read_pos:read_pos + cigartuple[1]]
                        if self.tech == "contig": continue
                        sub_read_quals[-1] += self.this_read_quals[read_pos:read_pos + cigartuple[1]]
                elif cigartuple[0] == 5:
                    continue
                read_pos += cigartuple[1]
            elif cigartuple[0] in [2, 3]:  # 2:D; 3:N: skip region of reference
                sub_read_str.extend([""] * cigartuple[1])
                if self.tech == "contig": continue
                sub_read_quals.extend([""] * cigartuple[1])
            else:
                print("MMMM===============")
                return -1
        self.this_read_list = sub_read_str
        # self.this_read_str = ""
        self.this_quals_list = sub_read_quals
        self.this_read_quals = ""

    def get_repeat_length(self, ms_start, ms_end):  # give a start and end
        query_repeat_length = len(
            "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
        return query_repeat_length

    def get_repeat_length_all_ms(self):  # return all microsatellite covered
        self.microsatellites_num = len(self.microsatellites)
        repeat_lengths = {}
        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = ms.end
            query_repeat_length = len(
                "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
            repeat_lengths[ms_id] = query_repeat_length
        self.repeat_lengths = repeat_lengths

    def get_quals(self, q_start, q_end):
        quals_list = self.this_quals_list[q_start - self.align_start - 1:q_end - self.align_start - 1]
        return quals_list, self.strand

    def get_ms_info_one_read(self):
        self.microsatellites_num = len(self.microsatellites)
        # print(self.read_id, self.microsatellites_num,len(self.support_microsatellites))
        read_muts = {}
        repeat_lengths = {}
        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = ms.end
            ms_start_pre = ms.start_pre
            ms_end_suf = ms.end_suf
            query_repeat_length = len(
                "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
            prefix = self.this_read_list[ms_start_pre - self.align_start:ms_start - self.align_start]
            suffix = self.this_read_list[ms_end - self.align_start:ms_end_suf - self.align_start]
            ms_content = self.this_read_list[ms_start - self.align_start:ms_end - self.align_start]
            mismatches = []
            deletions = []
            insertions = []
            pos_based_info = {}
            pos_deletion = []
            ref_pos = ms_start_pre - 2
            for pot in range(ms_start_pre - self.align_start - 1, ms_end_suf - self.align_start):
                ref_pos += 1
                this_read_base = self.this_read_list[pot]
                this_ref_base = self.this_ref_list[pot]
                if this_read_base == this_ref_base:
                    continue
                else:
                    if ref_pos < ms_start_pre:
                        band = [1]
                    elif ref_pos < ms_start:
                        band = [2]
                    elif ref_pos < ms_end:
                        band = [3]
                    elif ref_pos < ms_end_suf:
                        band = [4]
                    else:
                        band = [5]
                    this_read_base_len = len(this_read_base)
                    this_ref_base_len = len(this_ref_base)
                    if this_read_base_len == this_ref_base_len:
                        mismatches.append([ref_pos, this_read_base, band])
                        pos_based_info[ref_pos] = ["SNV", this_ref_base, this_read_base, band]
                    else:
                        if this_read_base_len < this_ref_base_len:
                            deletions.append([ref_pos, this_read_base, band])
                            # pos_based_info[ref_pos] = ["DEL", this_ref_base, ""]
                            pos_deletion.append([ref_pos, this_ref_base, band[0]])
                        else:
                            if ref_pos == ms_start - 1 and this_read_base_len > ms.repeat_unit_len:
                                # TODO
                                band = [3]
                            insertions.append([ref_pos, this_read_base, band])
                            if band[0] != 3:
                                pos_based_info[ref_pos] = ["INS", "", this_read_base[1:], band]

            if len(pos_deletion) > 0:
                deletion_start = pos_deletion[0][0]
                deletion_len = 1
                del_str = pos_deletion[0][1]
                band = [pos_deletion[0][2]]
                for i in range(1, len(pos_deletion)):
                    if pos_deletion[i - 1][0] + 1 == pos_deletion[i][0]:
                        deletion_len += 1
                        del_str += pos_deletion[i][1]
                        band.append(pos_deletion[i][2])
                    else:
                        pos_based_info[deletion_start] = ["DEL", del_str, "", set(band)]

                        deletion_len = 1
                        deletion_start = pos_deletion[i][0]
                        del_str = pos_deletion[i][1]
                        band = [pos_deletion[i][2]]
                # print(del_str,deletion_start,band)
                if len(set(band)) > 1 or list(set(band))[0] != 3:
                    pos_based_info[deletion_start] = ["DEL", del_str, "", set(band)]
            read_str = self.this_read_list[ms_start_pre - 1 - self.align_start:ms_end_suf + 2 - self.align_start]
            # read_muts[ms_id] = Read_Mutation(self.read_id, repeat_length=query_repeat_length, strand=self.strand,
            #                                  hap=self.hap,
            #                                  mismatches=mismatches, deletions=deletions, insertions=insertions,
            #                                  prefix=prefix, suffix=suffix, ms_content=ms_content,
            #                                  pos_based_info=pos_based_info,
            #                                  read_str=read_str
            #                                  )
            repeat_lengths[ms_id] = query_repeat_length
        self.repeat_lengths = repeat_lengths
        self.mut_info = read_muts

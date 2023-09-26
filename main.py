#!/usr/bin/env python 
"""
Project: TRgenotype
Script: main.py
Author: Peng Jia
E-mail: pengjia@stu.xjtu.edu.cn
Time : 2023/5/31
Description: TODO
"""
import os
import sys

curpath = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.append(os.path.dirname(curpath))
#
# from mshunter.genotype import *
# from mshunter.qc import *
from src.paras import args_process
from src.units import *
from src.genotype import genotype

# logger.info(" ".join(sys.argv))


def main():
    """
    Main function.
    :return:
    """
    global_init()
    arg = args_process()
    # print("000000000000",arg,)
    if arg:
        parase = arg.parse_args()
        if parase.command == "genotype":
            genotype(parase)
            # genotype_ngs(parase)
        # if parase.command == "qc":
        #     qc(parase)

        # if parase.command == "ngs":
        #     # genotype(parase)
        #     genotype_ngs(parase)


if __name__ == "__main__":
    main()
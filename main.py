#!/usr/bin/env python 
"""
Project: HapTR
Script: main.py
Author: Peng Jia
E-mail: pengjia@stu.xjtu.edu.cn
Time : 2023/5/31
Description: TODO
"""
import os
import sys
import time

curpath = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.append(os.path.dirname(curpath))
#
# from mshunter.genotype import *
# from mshunter.qc import *
from src.init import args_process
from src.units import *
from src.genotype2 import genotype
from src.train import Train
from src.genotype import Genotype
from src.param import Param



def main():
    """
    Main function.
    :return:
    """
    global_init()
    param = Param()
    # print("000000000000",arg,)
    if param.args_cmd() and param.args_init():

        if param.command == "genotype":
            genotype = Genotype(param)
            genotype.run()
        elif param.command == "train":
            train = Train(param)
            train.run()


if __name__ == "__main__":
    a = time.time()

    main()
    b = time.time()
    logger.info(f"Total cost: {round(b - a)} s")

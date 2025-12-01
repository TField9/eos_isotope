import numpy as np
import torch
from torch import nn
from torch.nn import TransformerEncoder, TransformerEncoderLayer
import os
import torch.optim as optim
import matplotlib.pyplot as plt
import ROOT
from torch.utils.data import TensorDataset, DataLoader, random_split
from tqdm import tqdm  # 引入tqdm用于显示进度条
import sys  # 添加sys模块用于重定向输出
import datetime  # 添加datetime用于时间戳
import traceback  # 添加traceback用于错误处理
import warnings
import argparse
import time
import numpy as np
from array import array

def norm_check():
    # 打开 ROOT 文件
    f = ROOT.TFile.Open("/eos/ams/user/s/selu/mdst/tianye/root/Be_MC_ML_residual_norm/sdst_be7_61_70_residual_norm.root", "READ")
    if not f or f.IsZombie():
        print("Error: Failed to open ROOT file")
        return
    
    # 获取树
    t = f.Get("tree")
    if not t:
        print("Error: Tree 'tree' not found")
        f.Close()
        return
    

    ek = ROOT.EK()


    # 设置分支地址

    t.SetBranchAddress("ek", ek)
    y_norm = np.zeros((3, 7), dtype=np.float32)  # 3x7 的二维数组
    
    # 设置分支地址
    t.SetBranchAddress("y_residual", y_norm)
    
    # 创建直方图
    h1 = ROOT.TH1D("h1", "y_residual_norm", 100, -5, 5)
    
    # 遍历所有条目
    nentries = t.GetEntries()
    for i in range(nentries):
        t.GetEntry(i)
        
        # 应用能量截断条件
        if ek.Ek[0] < 0.25 or ek.Ek[0] > 1.5:
            continue
        
        # 填充直方图 (detector=0, layer=6)
        h1.Fill(y_norm[0][6])
    
    # 创建画布并绘制直方图
    c1 = ROOT.TCanvas("c1", "c1", 800, 600)
    h1.Draw()
    
    # 保存为 PDF
    c1.Print("/eos/ams/user/s/selu/mdst/tianye/transformer/pdf/norm_check.pdf")
    
    # 清理资源
    f.Close()

# 执行函数
if __name__ == "__main__":
    ROOT.gInterpreter.ProcessLine('#include "/afs/cern.ch/work/s/seludev/private/tianye/CC/sdst.h"')
    norm_check()
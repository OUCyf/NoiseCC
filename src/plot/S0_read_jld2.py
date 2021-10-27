#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 20:09:08 2021

######################
##### read h5 ########
######################
# 1.read h5-file
h5_file = h5py.File(files[1],'r')

# 2.show all keys in h5-file
h5_file.keys()

# 3.循环读取所有 keys in h5-file
for key in h5_file.keys():
    onekey = key
    onekey_name = h5_file[key].name

# 4.已知某个group的 key "NN"
h5_file["NN"]
h5_file["NN"].keys()
f_dict = dict(h5_file["NN"])
f_dict.keys() # 所有的keyword


# 5.读取 group 的 datasets
data = f_dict["data"][()] # 建议 
data = f_dict["data"].value # data 是 numpy 的 ndarray 多维数组模式
trace = data[0] # 某一道数据

# 6.读取 group 的 Int Float 类型
baz = f_dict["baz"].value
baz = h5_file["NN"]["baz"].value

# 7.读取 group 的 字符串
# encode的作用是将unicode编码转换成其他编码的字符串，如str2.encode(‘utf8’)，表示将unicode编码的字符串str2转换成utf8编码
comp = h5_file["NN"]["comp"].value[0].decode('utf-8')

# 8. 关闭文件
f_dict.close()

######################
##### write h5 ########
######################

@author: yf
"""
#%%
import numpy as np
import h5py
import os
import glob



#%% 1. set parameter
file = "../../data/BJ.081_BJ.084__2020_04_11_00_00_00T2021_04_13_00_00_00__all.jld2"
chan = "NN"
dt = 0.005


#%% 2. read h5
# open file
f = h5py.File(file,'r')

# read data
data = f[chan]["data"][0] 

# read parameters
azi = f[chan]["azi"][()]
baz = f[chan]["baz"][()]
maxlag = f[chan]["maxlag"][()]
cc_len = f[chan]["cc_len"][()]
cc_step = f[chan]["cc_step"][()]
corr_type = f[chan]["corr_type"][()]
comp = f[chan]["comp"][()]
dist = f[chan]["dist"][()]   # dist = f[chan]["dist"].value
lat = f[chan]["lat"][()]
lon = f[chan]["lon"][()]
N_glob = f[chan]["N_glob"][()]
N_read = f[chan]["N_read"][()]
N_good = f[chan]["N_good"][()]
name = f[chan]["name"][()][0].decode('utf-8')

# close h5-file
f.close()



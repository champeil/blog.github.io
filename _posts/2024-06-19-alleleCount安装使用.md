---
layout:     post
title:      alleleCount安装使用
date:       2024-06-19
author:     champeil
catalog: true
tags:
    - alleleCount
    - htslib
    - 软件安装使用
    - linux
---

# 前言
- 在ASCAT检测CNV的过程中，需要用到alleleCount这么一个软件来计数基因组数据（bam文件等）中等位基因，有两种安装方法：
  - `mamba install cancerit-allelecount`安装的allelecount软件
  - 直接源码编译
```bash

```
- 在源码编译alleleCount的过程中其中一步htslib，里面有多种库需要依赖，报错最多的也是htslib配置过程
  - 在conda环境中如果安装了对应的配置文件，如果直接setup的话，则默认寻找的是自身的配置文件而不是conda环境配置好的文件，所以我们需要调整一下配置文件的位置
  - 可以直接使用conda安装了htslib以后重定向到htslib对应环境的库

# 安装过程
## alleleCount安装
```bash
# 1. download alleleCount
https://github.com/cancerit/alleleCount/releases/tag/v4.2.1
# 2. de-compress and enter the dir
tar -zxvf v4.2.1.tar.gz
# 3. modify the file of configure and setup
vim ./build/opt-build.sh
./configure --enable-plugins --enable-libcurl --with-libdeflate --prefix=$INST_PATH \
  CPPFLAGS="-I PATH/miniconda/include -I PATH/htslib/include" \ 
  LDFLAGS="-L PATH/miniconda/miniconda/lib -L PATHminiconda/htslib/lib  -Wl,-R${INST_PATH}/lib" \ 
  --disable-libcurl
# change the CPPFLAGS and LDFLAGS like this, to re-position the source file
# libcurl is for connection of the Internet for htslib, which is not used here and cannot installed by conda? so disabled
# 4. configure the finished 
vim ~/.bashrc
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/alleleCount/alleleCount-4.2.1/lib # add the export and save
source ~/.bashrc
```

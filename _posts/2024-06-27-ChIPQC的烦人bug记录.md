---
layout:     post
title:      ChIPQC的烦人bug记录
date:       2023-06-27
author:     champeil
catalog: true
tags:
    - chipseq
    - atacseq
    - chipqc
    - r
    - software
    - bug
---

# introduction
- 在R4.2.2版本中，我安装了ChIPQC，版本是1.40.0，在运行的过程中输入了csv文件，但是一直在报错
- 经过各种查资料，包括版本回退等都没有用，都会报各种错，所以以1.40.0为例，通过修改R代码暴力拆解错误原因

# 错误寻找
## 错误一
`Error in h(simpleError(msg, call)) : 
  error in evaluating the argument 'x' in selecting a method for function 't': 'names' attribute [9] must be the same length as the vector [7]`
- 默认所有的文件都是按照教程准备好并且都没有修改内容
- 在使用csv三部分一部分一部分排除以后，发现去除实验组、对照组的bam文件以及相对应的ID会报错，而去除peak以后则不报错，提示是peak文件的问题
- 在源代码的`R/ChIPQC_IF.R`中，ChIPQC函数的`res = new("ChIPQCexperiment",Samples=samples,DBA=experiment,annotation=annotation)`出错
  - 对应的代码在`R/ChIPQCexperiment-class.R`里面的`showChIPQCexperiment`函数，而出现错误的地方为`print(QCmetrics(object))`
    - 对应的代码在`R/ChIPQCexperiment-class.R`里面的`setMethod("QCmetrics", "ChIPQCsample", function(object)`函数
      - 在输出res的过程中发现有两个输出的是空值，导致长度跟名字不同，所以报错，具体是`fragmentlength(object,width=readlength(object))`与`signif(RelativeCrossCoverage(object),3)`这两个值
- 尝试解决
  - 查看函数
    - 对应到`R/ChIPQCexperiment-class.R`里面的`setMethod("RelativeCrossCoverage", signature(object="ChIPQCsample"), function(object)`
    - 对应到`R/ChIPQCexperiment-class.R`里面的`setMethod("fragmentlength", "ChIPQCsample", function(object,width)`
      - 由于readlength大小超出了crosscoverage的长度，所以返回的是空值，也就是`crosscoverage(object)[-seq(1:(2*readlength(object)))]`这里返回空值
      - 进而导致了MaxShift为空值
      - 尝试直接将`fragmentlength(object,width=readlength(object))`与`signif(RelativeCrossCoverage(object),3)`直接设置为0，结果成功
  - 关于`crosscoverage(object)[-seq(1:(2*readlength(object)))]`返回空值
  

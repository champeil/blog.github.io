---
layout:     post
title:      ChIPQC的烦人bug记录
date:       2024-06-27
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
          - 问题确定：fragmentlength无法进行检测，也就是空值
  - 关于`crosscoverage(object)[-seq(1:(2*readlength(object)))]`返回空值
      - 在转录因子或者是窄的表观遗传标记中，watson与crick reads，也就是正链与负链的reads往往在结合位点附近堆积，这种堆积现象可以衡量chipseq效率，也就是堆积越多，效率越高
      - 如果正链的峰从5'往3'方向移动，两个峰之间的重叠越来越多，而其峰的总体基因组覆盖度便越来越低
      - chipqc通过每一次正链从5'到3'移动1bp来计算总体基因组覆盖度（负链不移动），并转换成cross-coverage分数
          - 第n次移动的crosscoverage分数=（未移动前基因组覆盖度-第n次移动的基因组覆盖度）/未移动前基因组覆盖度
          - 以移动的nbp为x轴，crosscoverage分数为y轴，便会得到一条曲线，这条曲线往往在fragmentlength处会有一个最大峰，而在readlength处会有一个小峰
          - fragment处的峰代表着chipseq在结合位点处成功富集，而readlength处的峰代表幻影峰
              - 假阳性的幻影蜂推测为转录因子与转录机制的一般粘性，也就是跟任何DNA结合并获得结合物
              - 或者是一些富集程度低的背景峰（下图）
          - 灰色区域则是识别fragmentlength对应的peak的时候需要排除的区域，在代码中也就是2倍的readlength，排除了幻影蜂以后，crosscoverage最大值为fragmentlength
      - 解释 
          - 测序的过程
              - 首先将DNA通过不同的方法，例如超声等，切分成各种大小的DNA片段，也就是fragment，insert length就是fragment length，然后我们将DNA小段，也就是fragment两端添加接头进行测序
              - 单端测序只从一端出发，双端则从fragment两端出发进行测序
              - 通常在二代测序中，我们建库的时候会对fragment进行跑胶质控，通常取300-500bp之间的fragment条带，而reads通常为150bp，也就是fragment的长度大于等于2*readlength的长度
          - 根据代码来看，crosscoverage的长度由shift决定，默认为300bp，大于等于2倍的readlength，所以理论上来说是有值的，双端都可以，那就更不用说50bp左右的单端了
          - 而在单端测序中，通常crosscoverage是用来衡量峰位移的长度的，因为chip结合下来的双链，在正链与负链测序后，峰都在结合位点附近聚集，所以需要峰位移，所以理论上crosscoverage plot在单端测序也管用，fragment length可以推断位移距离
          - 所以在双端测序中，fragment length为fragment长度，而在单端测序中，则为位移距离
      - 所以问题就出在了readlength(object)中
  - 关于`readlength(object)`返回read length中，我查看了代码`R/sampleQC.r`，里面使用了`GenomicAlignments`R包读取bam文件，并且将width列的前1000个值的均值用作read长度
      - 背景：我使用STAR进行比对的
      - 使用这个R包读取发现，里面的width出现17000长度的现象，这个很明显不合理，但是我也发现了cigar中以及njun列中分别是比对情况以及有多少个可变剪接，而width是reads在基因组中的覆盖长度
      - 所以意识到，可能是STAR这个软件将可变剪接等事件包含进去导致的width长度拉长，而width指的是比对以后包括中间的可变剪接、插入、缺失等事件的比对情况总长度，所以导致readlength大于预期的reads长度
      - 通过去除可变剪接，或者是直接将readlength设置成50，是跑成功的，进一步验证想法，所以为了兼容STAR的比对情况，需要再想一个办法
      - 而留意到，qwidth列是将这些事件去除以后，单纯的readlength长度，也就是reads的原始长度，而这个似乎比较符合我的想法，所以建议是将width转换成qwidth列
  - 解决办法
      - `R/sampleQC.r`中的`readlength=round(mean(width(temp[1:tocheckforreads])))`换成`readlength=round(mean(GenomicAlignments::qwidth(temp[1:tocheckforreads])))`即可
      - 然后`R CMD build ./ChIPQC --no-build-vignettes`编译，然后当前文件夹下的包就行了

![image](https://github.com/champeil/champeil.github.io/assets/33405808/63831b79-e2bb-4daa-8245-07d112844391)

## 错误二
- 严格来说也不算是错误，只不过是一个适配的东西
- 在后续的图片绘制过程中，有时候想要将input样本纳入进去查看chip实验的富集效率，但是按照DBA的那种csv的文件读入方法的话，在后面会将input样本给去掉，也就是说如果按照DBA的那种输入的话，input样本是没有peaks信息，也就是无法进行后续的查看的
- 看代码的时候发现samplelist与controllist进行读取与构建对象以后，samplelist列表中结构是这样的：`samplelist：samplename：bam + peak`以及`samplelist：controlname：bam`，所以在后面成功读入以后发现input那边有很多关于peak富集的参数都是NA值，这就很难搞了
- 其实关于这个chipqc应该是提供了一个方法解决，就是`taximofen`的那个变量，他是没有control列的，就是将input用做ID里面，当做sample来看了，但是觉得很麻烦
- 所以尝试改了一下代码，位置在`ChIPQC/R/ChIPQC_IF.R`里面95行开始

```r
# raw
#  samplelist = NULL
#  controlist = NULL
#  for(i in 1:nrow(meta)) {
#    newrec = NULL
#    newrec$peaks = experiment$peaks[[i]]
#    if(nrow(newrec$peaks)==0) {
#      newrec$peaks = NULL
#    }
#    newrec$bam   = as.character(meta$bamRead[i])
#    samplelist   = listadd(samplelist,newrec)
#    if(!is.null(meta$bamControl[i])) {
#      if(!is.na(meta$bamControl[i])) {
#        if(meta$bamControl[i]!="") {
#          savenames = names(controlist) 
#          controlist = listadd(controlist,as.character(meta$bamControl[i]))
#          names(controlist) = c(savenames,as.character(meta$Control[i]))
#        }
#      }
#    }
#  }
#  controlist = controlist[!duplicated(controlist)]
#  controls = 0
#  for(cfile in controlist) {
#    if (!(cfile %in% as.character(meta$bamRead))) {
#      newrec = NULL
#      newrec$bam = cfile
#      newrec$peaks=NULL
#      samplelist = listadd(samplelist,newrec)
#      controls = controls+1
#    }
#  }
#  
#  names(samplelist) = unique(c(rownames(meta),names(controlist)))

# change to
  samplelist = NULL
  controlist = NULL
  for(i in 1:nrow(meta)) {
    newrec = NULL
    newrec$peaks = experiment$peaks[[i]]
    if(nrow(newrec$peaks)==0) {
      newrec$peaks = NULL
    }
    newrec$bam   = as.character(meta$bamRead[i])
    savenames_sample = names(samplelist)
    samplelist   = listadd(samplelist,newrec)
    names(samplelist) = c(savenames_sample,as.character(meta$ID[i]))
    if(!is.null(meta$bamControl[i])) {
      if(!is.na(meta$bamControl[i])) {
        if(meta$bamControl[i]!="") {
          newrec_cont = NULL
          newrec_cont$peaks = experiment$peaks[[i]]
          if(nrow(newrec_cont$peaks)==0) {
            newrec_cont$peaks = NULL
          }
          newrec_cont$bam   = as.character(meta$bamControl[i])
          savenames = names(controlist)
          controlist = listadd(controlist,newrec_cont)
          names(controlist) = c(savenames,as.character(meta$Control[i]))
        }   
      }   
    }
  }
  samplelist <- c(samplelist,controlist)

# 原始代码的input是不会读取到peak列的，就只有bam变量，所以我们需要对samplelist结构进行修饰，也就是control也设置成bam与peak变量
```

## 错误三
- 检测peak的时候，并不是所有染色体都有peak，一旦在指定的染色体中没有检测到peak，则会报错`Error in data.frame(Counts, bedRangesSummitsTemp): arguments imply differing number of rows: 1, 0`，提示bedRangesSummitsTemp也就是peak的数目为0，所以长度对不上
- 所以需要针对peak进行检查，如果没有peak的话，则返回的是空的GRange对象
- 代码位置：`ChIPQC/R/sampleQC.R`的302行
```r
# raw
if(!is.null(bedFile)){
    AvProfile <- colMeans(CoverageMatrix)
    NormAvProfile <- (AvProfile/FlagTagCounts[4])*1e6
    elementMetadata(bedRangesTemp) <- data.frame(Counts,bedRangesSummitsTemp)

# change to
if(!is.null(bedFile)){
    AvProfile <- colMeans(CoverageMatrix)
    NormAvProfile <- (AvProfile/FlagTagCounts[4])*1e6
    if(length(bedRangesSummitsTemp)==0){
        bedRangesTemp <- GRanges()
    }else{
         elementMetadata(bedRangesTemp) <- data.frame(Counts,bedRangesSummitsTemp)
    }

```













  

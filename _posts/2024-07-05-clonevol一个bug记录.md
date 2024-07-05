---
layout:     post
title:      clonevol一个bug记录
date:       2024-07-01
author:     champeil
catalog: true
tags:
    - WES
    - WGS
    - MPTevol
    - clonevol
    - r
    - software
    - bug
---

# introduction
- 在使用青建师兄的metevol构建进化树的过程中，发现引用的`clonevol`R包里面`inferClonalTrees`函数不断地报错，错误为`Error in generate.boot(variants, vaf.col.names = vaf.col.names, depth.col.names = depth.col.names,  : 
  Input error: cluster column name does not appear in variants file.`
- 意思就是`cluster.col.name`的列找不到了，寻思输入文件的cluster列的名字就是这个，为啥会报错

# process
- 上网查看了clonevol的github中，有人也提出来这个问题，他说如果将Cluster换成默认的小写就会成功
- 所以尝试将Cluster列名换成cluster，发现成功了
- 意义不明的小bug......

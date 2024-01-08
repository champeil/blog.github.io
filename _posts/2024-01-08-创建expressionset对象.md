---
layout:     post
title:      创建expressionset对象
date:       2023-01-08
author:     champeil
catalog: true
tags:
    - expressionset
    - Biobase
    - R_package
    - expression_data
    - R_object
---

# 创建expressionset对象 [来源:Biobase](https://www.bioconductor.org/packages/release/bioc/html/Biobase.html)
## 介绍
- Biobase 包含标准化数据结构来表示基因组数据。 ExpressionSet 对象旨在将多个不同的信息源组合成一个方便的结构
- expressionset主要包含以下内容：
  - assayData: 基因/探针表达矩阵
  - phenoData: 记录样本表型信息的数据框
  - featureData: 记录芯片或者是技术等信息
  - experimentData: 记录实验相关信息

## 创建expressionset对象
``` R
# the script is for constructing the expressionset objects
# necessary data
- assayData: expression data matrix with the row as gene names and col as sample names
- phenoData: dataframe to record each samples' pheno with the rownames are the same with assayData and col names are pheno name

# script
expressionset2 <- ExpressionSet(
  assayData = matrix,
  phenoData = new(
    "AnnotatedDataFrame",
    data=tibble(
      SampleID=colnames(matrix)
    ) %>%
      dplyr::left_join(clinical_data2,by=c("SampleID"="Tumor_Sample_Barcode")) %>%
      tibble::column_to_rownames(var="SampleID"),
        varMetadata=data.frame(labelDescription=c(colnames(clinical_data2)[-1]),row.names=c(colnames(clinical_data2)[-1]))
      ),
  annotation=set_2
)
```

## 参考
- Huber W, Carey VJ, Gentleman R, Anders S, Carlson M, Carvalho BS, Bravo HC, Davis S, Gatto L, Girke T, Gottardo R, Hahne F, Hansen KD, Irizarry RA, Lawrence M, Love MI, MacDonald J, Obenchain V, Ole's AK, Pag'es H, Reyes A, Shannon P, Smyth GK, Tenenbaum D, Waldron L, Morgan M (2015). “Orchestrating high-throughput genomic analysis with Bioconductor.” Nature Methods, 12(2), 115–121. http://www.nature.com/nmeth/journal/v12/n2/full/nmeth.3252.html.

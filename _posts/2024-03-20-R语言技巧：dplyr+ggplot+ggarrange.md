---
layout:     post
title:      R语言技巧：dplyr+ggplot+ggarrange
date:       2024-03-20
author:     champeil
catalog: true
tags:
    - R
    - dplyr
    - ggplot
    - R语言技巧
---

# 前言
- 在数据框的基础上根据某一列进行分组，并各自绘制ggplot图，最后使用ggarrange合并图谱

# 例子

```r
# author: laojp
# time: 2024.03.20
# position: SYSUCC bioinformatic platform

rnaseq %>%
  dplyr::select(Tumor_Sample_Barcode,external_gene_name,fpkm) %>%
  dplyr::distinct(Tumor_Sample_Barcode,external_gene_name,.keep_all = TRUE) %>%
  dplyr::filter(external_gene_name!="") %>%
  tidyr::pivot_wider(values_from = fpkm,names_from = Tumor_Sample_Barcode,values_fill = 0) %>%
  dplyr::left_join(pathway,by=c("external_gene_name"="encode_gene")) %>%
  dplyr::filter(!is.na(pathway)) %>%
  dplyr::full_join(
    TCGA_rnaseq %>%
      dplyr::select(Tumor_Sample_Barcode,gene_name,fpkm) %>%
      dplyr::mutate(Tumor_Sample_Barcode=str_sub(Tumor_Sample_Barcode,1,12)) %>%
      dplyr::distinct(Tumor_Sample_Barcode,gene_name,.keep_all = TRUE) %>%
      dplyr::filter(gene_name!="") %>%
      tidyr::pivot_wider(values_from = fpkm,names_from = Tumor_Sample_Barcode,values_fill = 0) %>%
      dplyr::left_join(pathway,by=c("gene_name"="encode_gene")) %>%
      dplyr::filter(!is.na(pathway)),
    by=c("external_gene_name"="gene_name","pathway","Molecule")
  ) %>% distinct() %>%
  dplyr::arrange(pathway,Molecule) %>%
  tidyr::gather(key = "Tumor_Sample_Barcode",value="fpkm",-external_gene_name,-Molecule,-pathway) %>%
  dplyr::left_join(common_clin,by="Tumor_Sample_Barcode") %>%
  dplyr::filter(!is.na(group)) %>%
  dplyr::group_by(Molecule) %>%
  do(
    plots=ggplot(.)+
      aes(x=external_gene_name,y=fpkm,color=group)+
      geom_boxplot()+
      coord_flip()+
      stat_compare_means(method = "t.test",aes(label=..p.signif..),label.x=1.5,hide.ns = TRUE)
  ) %>% ungroup %>%

rnaseq %>%
  dplyr::select(Tumor_Sample_Barcode,external_gene_name,fpkm) %>%
  dplyr::distinct(Tumor_Sample_Barcode,external_gene_name,.keep_all = TRUE) %>%
  dplyr::filter(external_gene_name!="") %>%
  tidyr::pivot_wider(values_from = fpkm,names_from = Tumor_Sample_Barcode,values_fill = 0) %>%
  dplyr::left_join(pathway,by=c("external_gene_name"="encode_gene")) %>%
  dplyr::filter(!is.na(pathway)) %>%
  dplyr::full_join(
    TCGA_rnaseq %>%
      dplyr::select(Tumor_Sample_Barcode,gene_name,fpkm) %>%
      dplyr::mutate(Tumor_Sample_Barcode=str_sub(Tumor_Sample_Barcode,1,12)) %>%
      dplyr::distinct(Tumor_Sample_Barcode,gene_name,.keep_all = TRUE) %>%
      dplyr::filter(gene_name!="") %>%
      tidyr::pivot_wider(values_from = fpkm,names_from = Tumor_Sample_Barcode,values_fill = 0) %>%
      dplyr::left_join(pathway,by=c("gene_name"="encode_gene")) %>%
      dplyr::filter(!is.na(pathway)),
    by=c("external_gene_name"="gene_name","pathway","Molecule")
  ) %>% distinct() %>%
  dplyr::arrange(pathway,Molecule) %>%
  tidyr::gather(key = "Tumor_Sample_Barcode",value="fpkm",-external_gene_name,-Molecule,-pathway) %>%
  dplyr::left_join(common_clin,by="Tumor_Sample_Barcode") %>%
  dplyr::filter(!is.na(group)) %>%
  dplyr::group_split(Molecule) %>%
  purrr::map(function(df){
    plots=ggplot(df)+
      aes(x=external_gene_name,y=fpkm,fill=factor(group,levels=c("PTC+HT","PTC-HT")))+
      geom_bar(stat="summary",fun="mean",position=position_dodge(),alpha=0.8)+
      labs(fill="group",title=unique(df$Molecule),x="")+
      stat_summary(fun.data='mean_sd',geom="errorbar",colour="black",position=position_dodge(.9),width=0.15)+
      stat_compare_means(method = "t.test",aes(label=..p.signif..),label.x=1.5,hide.ns = TRUE)+
      theme_bw()+
      theme(
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        panel.grid = element_blank(),
        plot.title = element_text(size=15,face = "bold",hjust = 0.5)
      )
    return(plots)
  }) %>%
  ggarrange(plotlist = .,common.legend = TRUE,align="hv",ncol=3,nrow=4)
  pull(plots) %>%
  ggarrange(plotlist = .)

```

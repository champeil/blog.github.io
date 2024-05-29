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
  pull(plots) %>%
  ggarrange(plotlist = .,common.legend = TRUE)

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

list.files("/annotate_peaks/",pattern="*.txt",recursive = FALSE,full.names = TRUE) %>%
  as.list() %>%
  purrr::set_names(basename(list.files("/annotate_peaks/",pattern="*.txt",recursive = FALSE)) %>% str_remove_all(pattern="_summit_50_anno.txt")) %>%
  purrr::imap(function(x,name){
    return(read.table(x,header=FALSE,sep="\t",skip=1) %>%
             dplyr::rename(peak=V1,seqnames=V2,start=V3,end=V4,strand=V5,height=V6,annotation=V8,distanceToTSS=V10) %>%
             dplyr::rename(V4=peak,V5=height) %>%
             dplyr::select(seqnames,start,end,strand,V4,V5,annotation,distanceToTSS) %>%
             dplyr::mutate(annotation=case_when(
               str_detect(annotation,pattern="intron") ~ "Intron",
               str_detect(annotation,pattern="exon") ~ str_extract(annotation, "exon \\d+ of \\d+"),
               str_detect(annotation,pattern="promoter-TSS") ~ "Promoter",
               str_detect(annotation,pattern="Intergenic") ~ "Intergenic",
               TRUE ~ "TSS"
             )) %>%
             tidyr::separate(annotation,sep="( | of )",into = c("type","first","of","second")) %>%
             dplyr::mutate(annotation=ifelse(is.na(second),type,ifelse(first==1,paste("First",type,sep=" "),ifelse(first==second,paste("Last",type,sep=" "),"Internal exon"))),
                           sample=name))
  }) %>% 
  purrr::reduce(bind_rows) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(all_num=length(sample)) %>% ungroup() %>%
  dplyr::group_by(sample,annotation) %>%
  dplyr::mutate(anno_num=length(sample),
                freq=anno_num/all_num) %>% ungroup() %>%
  {
    df <- .
    list(
      plot1=ggplot(data=df,aes(x=sample,y=freq,fill=annotation)) +
        geom_col(position="fill",width=1)+
        theme_bw()+
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              legend.text = element_text(size=12),
              legend.title = element_text(size=12)),
      plot2=ggplot(data=df %>% dplyr::filter(str_detect(annotation,pattern="exon")),
                   aes(x=sample,y=anno_num,fill=annotation))+
        geom_col(position = "fill",width = 1)+
        theme_bw()+
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              legend.text = element_text(size=12),
              legend.title = element_text(size=12))
    )
  } %>%
  ggarrange(plotlist = .,ncol=2,nrow=1,align="hv")
  

```

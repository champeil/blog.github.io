---
layout:     post
title:      创建GATK mm10 resource bundle
date:       2023-06-26
author:     champeil
catalog: true
tags:
    - wes
    - wgs
    - wxs
    - software
    - gatk
    - reference
    - mm10
    - somatic
---

# 前言
- GATK resource bundle针对人类测序数据建立一系列的标准文件，包括已知的snp以及indel信息以用于BaseRecalibreator，RealignerTargetCreator，IndelRealigner。该方法用于使用UCSC mm10文件创建相似的文件
- BQSR步骤文件只需要突变位点用于跳过真实的突变
- getpileupsummaries步骤需要gemline位点并带有allele frequency
- PON步骤为了消除人工与技术误差以及germline，并且根据gatk教程中，其是通过mutect2对normal构建来的

# 文件
- bqsr dbsnp：NCBI
	- known snp sites: ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/
	- known indel sites: https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz
- germline source：sanger mouse
	- https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
- pon
	- 参考[5]
	- 首先小鼠没有太多的normal样本供我们做
	- 其次网上没有一个比较可信的pon文件，甚至连pon文件都没有
	- 所以我们只能使用germline当做pon进行过滤

## bqsr snp and indel source from ncbi and sanger

```shell
# download dbsnp vcf file from ncbi, combine and modify the name
wget --recursive --no-parent --no-directories --accept vcf*vcf.gz ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/
# add chr prefix to each chromosome （convert GRCm38 to mm10 format）
for vcf in $(ls -1 vcf_chr_*.vcf.gz) ; do
  vcf_new=${vcf/.vcf.gz/.vcf}
  echo $vcf
  zcat $vcf | sed 's/^\([0-9XY]\)/chr\1/' > $vcf_new
  rm -fv $vcf
done
# combine all dbsnp vcf file
vcf_file_string=""
for vcf in $(ls -1 vcf_chr_*.vcf) ; do
  vcf_file_string="$vcf_file_string -V $vcf"
done
echo $vcf_file_string

java -Xms16G -Xmx16G -cp ${gatk_path}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
-R genome.fa $vcf_file_string -out dbsnp.vcf
```

- mouse indel文件则从Sanger MGP（Sanger Mouse Genetics Programme）下载处理

```shell
# download all MGP indels (this file is 5/2015 released)
wget ftp://ftp-mouse.sanger.ac.uk/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz 

### first state : change to mm10 format
# filter for passing variants with chr add , adjust header
zcat mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 \
| grep -v "#contig" | grep -v "#source" \
> mgp.v5.merged.indels.dbSNP142.normed.pass.vcf

# keep only passing and adjust chromosome name
zcat mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | grep -v "^#" | cut -f 1-8 \
| grep -w "PASS" | sed 's/^\([0-9MXY]\)/chr\1/' \
>> mgp.v5.merged.indels.dbSNP142.normed.pass.vcf

# sort and generate index
java -Xms16G -Xmx16G -jar ${PICARD_ROOT}/picard.jar SortVcf VERBOSITY=WARNING \
	I=mgp.v5.indels.pass.chr.vcf \
	O=mgp.v5.indels.pass.chr.sort.vcf
```

## mouse germline
- 因为ncbi中dbsnp下载的snp文件没有GT信息，也就是无法计算突变的af
- 方法
	- 从mouse sanger中下载snp与indel的germline文件，在vcf文件中选择对应的种系小鼠+pass filter
	- sort+index以后使用bcftools concat合并两者文件
	- sort+index文件
	- 使用bcftools增加AF，并且使用bgzip+tabix压缩与index文件

```shell
wget -c -t 0 https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
wget -c -t 0 mgp.v5.merged.indels.dbSNP142.normed.vcf.gz
```

```perl
#!/usr/bin/perl -w 
use strict;
use warnings;
use File::Basename;

=head1
# this script is for modifying the germline vcf file from mouse sanger source
# author: laojp
# time: 2023.04.07
# position: SYSUCC bioinformatic platform
# usage: perl mouse_sanger_germline.pl [input snp.gz file] [input indel.gz file] [mouse line] [output snp.vcf file] [output indel.vcf file]
# ps: only for mouse sanger germline variation file
# ps: please copy the pl and next shell file to the reference and run
# update: will considerate the transferation from grcm38 to mm10 [format] 
# update: will considerate the combination of different part of ncbi dbsnp file [ncbi dbsnp file]
=cut

# judge the parameter is exist or not
my $vcf_file = shift or die "Please input the vcf file \n $!";
my $indel_file = shift or die "Please input the indel file \n $!";
my $mouse_line = shift or die "Please refer the mouse line \n $!";
my $out_file = shift or die "Please refer the output file \n $!";
my $out_indel = shift or die "Please refer the output indel file \n $!";
my $script_dir = dirname(__FILE__);

# read file --- outputhead --- filter the PASS and choose the mouse line --- zip and create index
sub modify_vcf {
	my %header;
	open(VCF, "-|", "gunzip -c $_[0]") or die "$!";
	open(OUT, ">", $_[1]) or die "$!";
	while(<VCF>){
		chomp;
		if($_ =~ /^#CHROM/){
			my @F = split(/\t/, $_);
			for(0..$#F){
				$header{$F[$_]} = $_;
			}
			print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${mouse_line}\n";
			last;
		}
		if($_ =~ /^#/){
			print OUT "$_\n";
		}
	}
	while(<VCF>){
		chomp;
		my @F = split(/\t/, $_);
		my @G = split(/:/, $F[$header{"$mouse_line"}]);
		if($F[$header{"FILTER"}] eq "PASS" && $G[0] ne './.'){
			my $output = join("\t", $F[$header{"#CHROM"}], $F[$header{"POS"}], $F[$header{"ID"}], $F[$header{"REF"}], $F[$header{"ALT"}], $F[$header{"QUAL"}], $F[$header{"FILTER"}], $F[$header{"INFO"}], $F[$header{"FORMAT"}], $F[$header{"$mouse_line"}]);	
			print OUT "$output\n";
		}
	}
	close(VCF);
	close(OUT);	
}

modify_vcf($vcf_file, $out_file); 
modify_vcf($indel_file, $out_indel);

`bash $script_dir/mouse_sanger_germline_next.sh $out_file`;
`bash $script_dir/mouse_sanger_germline_next.sh $out_indel`;
`/home/laojp/software/bcftools/bin/bcftools concat $out_file $out_indel -O z -a -o germline_merge.vcf.gz`
`bash $script_dir/mouse_sanger_germline_next.sh germline_merge.vcf.gz`
`/home/laojp/software/bcftools/bin/bcftools +fill-tags germline_merge_sort.vcf.gz -- -t AF | bgzip > germline_merge_sort_af.vcf.gz`
`tabix germline_merge_sort_af.vcf.gz`

```

```shell
# mouse_sanger_germline_next.sh
#!/bin/bash
# this script is for handle the sanger_mouse vcf file which outputed from mouse_sanger_germline.pl
# author: laojp
# time: 2023.04.07
# position: SYSUCC bioinformatic platform
# usage: bash mouse_sanger_germline_next.sh [input vcf file] 

# update: will update together with mouse_sanger_germline.pl

#sort --- generate index
filedir=$(dirname $(readlink -f ${1}))
if [[ ${1} =~ "gz" ]]; then
	name=$(basename ${1} .vcf.gz)
	/home/laojp/software/gatk_4.2.6.1/gatk-4.2.6.1/gatk SortVcf \
		-I ${1} \
		-O ${filedir}/${name}_sort.vcf.gz
else
	name=$(basename ${1} .vcf.gz)
	gzip ${1}
	/home/laojp/software/gatk_4.2.6.1/gatk-4.2.6.1/gatk SortVcf \
		-I ${1}.gz \
		-O ${filedir}/${name}_sort.vcf.gz
fi

```
# 支持一哈
If you like my tutorial, please cite:
- [Zhang M, Wen H, Liang M, Qin Y, Zeng Q, Luo D, Zhong X, Li S. Diagnostic Value of Sylvian Fissure Hyperechogenicity in Fetal SAH. AJNR Am J Neuroradiol. 2022 Apr;43(4):627-632. doi: 10.3174/ajnr.A7449. Epub 2022 Mar 10. PMID: 35272984; PMCID: PMC8993207.](https://pubmed.ncbi.nlm.nih.gov/35272984/)

# 参考信息
\[1]. https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md
\[2]. https://www.biostars.org/p/84501/#84582
\[3]. https://www.biostars.org/p/300996/
\[4]. https://gatk.broadinstitute.org/hc/en-us/community/posts/360078092732-Common-germline-variants-sites-VCF-for-GetpileupSummaries-for-mouse-data
\[5].  Lange S, Engleitner T, Mueller S, Maresch R, Zwiebel M, González-Silva L, Schneider G, Banerjee R, Yang F, Vassiliou GS, Friedrich MJ, Saur D, Varela I, Rad R. Analysis pipelines for cancer genome sequencing in mice. Nat Protoc. 2020 Feb;15(2):266-315. doi: 10.1038/s41596-019-0234-7IF: 17.021 Q1 . Epub 2020 Jan 6. PMID: 31907453IF: 17.021 Q1 .

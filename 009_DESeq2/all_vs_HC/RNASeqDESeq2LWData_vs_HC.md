---
title: "DESeq2_RNASeq_LWData_vs_HC"
output: 
  html_document:
    #css: styles.css
    toc: true
    toc_depth: 4
    number_sections: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    fig_caption: true
    df_print: paged
    #code_folding: hide
    keep_md: true
editor_options: 
  markdown: 
    wrap: 72
---



# Overview

In this analysis we will perform differential analysis using DESeq2. We
will look in some details of DESeq2.

# Read in CountTables

assume you are in directory 009_DESeq2 (project folder structure
needed), you need the following path to read in count tables:


``` r
ct.dir <- "../../008_FeatureCountsout/ENSEMBLE/" # set input dir = counttables
```

now list files to get vector of path to count files. I have also
unstranded und stranded counted count files in my directory, so I need
to define that I only want to list the reverse stranded count files,
which end with \_\_S2.


``` r
v.f <- list.files(path = ct.dir, 
                  pattern = "__S2$", 
                  full.names = T) # full path to Count tables

# Identify the indices of files 4, 5, and 6 for GMSM (would be 7,8,9 for HC)
indices <- c(7, 8, 9) # i want GMSM in the first columns! as my reference for analysis 1                       # later for the second analysis i want HC as my reference 

# Order the files in v.f (see comments above)
v.f <- c(v.f[indices], v.f[-indices])

rm(indices)
v.f
```

```
##  [1] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#HC1_S109919_1_SE__S2"   
##  [2] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#HC2_S109922_1_SE__S2"   
##  [3] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#HC3_S109921_1_SE__S2"   
##  [4] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#DM1_S109918_1_SE__S2"   
##  [5] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#DM2_S109917_1_SE__S2"   
##  [6] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#DM3_S109920_1_SE__S2"   
##  [7] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#GMSM1_S109914_1_SE__S2" 
##  [8] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#GMSM2_S109916_1_SE__S2" 
##  [9] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#GMSM3_S109915_1_SE__S2" 
## [10] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#OA1789_S109912_1_SE__S2"
## [11] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#OA1887_S109913_1_SE__S2"
## [12] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#OA1966_S109910_1_SE__S2"
## [13] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#OA2022_S109909_1_SE__S2"
## [14] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#OA2066_S109908_1_SE__S2"
## [15] "../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#OA2272_S109911_1_SE__S2"
```

extract short name from file -\> will be used as new column name in
count table


``` r
# string modification to extract the phenotype adapted to Lukas Data (inkl. nr. of biological replicate)
v.n <- gsub("_S1099.._1_SE__S2", "", 
            gsub("../../008_FeatureCountsout/ENSEMBLE//BSF_1140_HNJ5KBBXY_1#", "", v.f))
```

create emtpy list for for-loop


``` r
l.f <- vector(mode = "list", 
              length = 0)
```

read in count tables, subset to EnsembleIds and Counts, make named
vector


``` r
for (i in seq_along(v.f)){
  df <- read.table(v.f[i], sep = "\t", header = T)[, c(1, 7)]
  v <- df[, 2]
  names(v) <- df[, 1]
  l.f[[v.n[i]]] <- v
}

# after this the variable df remains from the last file (not needed, can be ignored, the loop processed each file (sample) and cached the sample data in df until the next file was processed, so basically just needed for processing)
# l.f is needed!! here the values for each file (samples are stored)
```

l.f is list of integer vectors

convert list to data.frame, and then to matrix


``` r
# here the list is converted to different data type, information is unchanged

df.c <- data.frame(l.f)

# convert df to matrix
m.c <- as.matrix(df.c)
# head = show first lines for quick consistency check (subsequently used)
head(m.c)
```

```
##                     HC1  HC2  HC3  DM1  DM2  DM3 GMSM1 GMSM2 GMSM3 OA1789
## ENSRNOG00000066169    0    0    0    0    0    0     0     0     0      0
## ENSRNOG00000070168    0    0    0    0    0    0     0     0     0      0
## ENSRNOG00000070901 1950 2140 1587 4004 3907 2207   984  2367  1632   1629
## ENSRNOG00000018029   26   19   12   34   14   54     6    30    22     29
## ENSRNOG00000031391    0    1    0    0    0    1     0     0     0      0
## ENSRNOG00000055342    0    0    0    0    0    0     0     0     0      0
##                    OA1887 OA1966 OA2022 OA2066 OA2272
## ENSRNOG00000066169      0      0      0      0      0
## ENSRNOG00000070168      0      0      0      0      0
## ENSRNOG00000070901   1614   1703   1498   2178   2325
## ENSRNOG00000018029     10      7     18     27     22
## ENSRNOG00000031391      0      0      1      0      1
## ENSRNOG00000055342      0      0      0      0      0
```

# create different tables for summarized experiment object

Here is how a summarized experiment object looks like
![link](https://www.researchgate.net/profile/Pascal-Martin-6/publication/274256021/figure/fig3/AS:294713816043527@1447276657416/Structure-of-summarizedExperiment-objects.png)

## Create colData Table for Summarized experiment object

colData Table contains metadata of samples, e.g. genotype, batch, cell
type, ...\
rownames of colData must be equal to colnames of assay (=count-) table


``` r
# grouping of biological replicates of the same phenotypes, needs to be adapted later to group the OA apropriately

df.m <- data.frame(geno=substr(v.n, start = 1, stop = 2), row.names = v.n)
```

Define reference level of genotype -\> important, because base level of
GLM will be reference level of genotype, and comparisions will be
relative to reference level.


``` r
# define ref. level the order of the level should be adapted for second analysis with HC
# as the reference (just to make sure to not mix something up) plus additional levels
# for the apropriate OA phenotypes e.g. OA1/OA2/OA3
df.m$geno <- factor(df.m$geno, levels = c("HC","DM","GM","OA"))
```

## Create rowData Table for Summarized experiment Object

### import and subset GTF file to columns you want to have in rowData

This is the GTF file we used also for mapping. So identifiers in count
table and GTF file fit together.\
e.g. gene symbol, ensemble gene id, ....


``` r
# GTF file used for STAR aligner, change if analysis is performed with the data from Star and feature counts made with RefSeq genome
gm <- readGFFAsGRanges("/home/fhwn.ac.at/201138/RNASeq_LW/000_referenceseq/ENSEMBL/Rattus_norvegicus.mRatBN7.2.109.gtf")
gm
```

```
## GRanges object with 1284429 ranges and 21 metadata columns:
##                      seqnames            ranges strand |   source       type
##                         <Rle>         <IRanges>  <Rle> | <factor>   <factor>
##         [1]                 1 36112690-36122387      - |  ensembl gene      
##         [2]                 1 36112690-36122387      - |  ensembl transcript
##         [3]                 1 36122324-36122387      - |  ensembl exon      
##         [4]                 1 36122324-36122387      - |  ensembl CDS       
##         [5]                 1 36121478-36121512      - |  ensembl exon      
##         ...               ...               ...    ... .      ...        ...
##   [1284425] JACYVU010000665.1       13375-13602      - |  ensembl exon      
##   [1284426] JACYVU010000665.1       13535-14898      - |  ensembl transcript
##   [1284427] JACYVU010000665.1       14554-14898      - |  ensembl exon      
##   [1284428] JACYVU010000665.1       13726-13831      - |  ensembl exon      
##   [1284429] JACYVU010000665.1       13535-13602      - |  ensembl exon      
##                 score     phase            gene_id gene_version gene_source
##             <numeric> <integer>        <character>  <character> <character>
##         [1]        NA      <NA> ENSRNOG00000066169            1     ensembl
##         [2]        NA      <NA> ENSRNOG00000066169            1     ensembl
##         [3]        NA      <NA> ENSRNOG00000066169            1     ensembl
##         [4]        NA         0 ENSRNOG00000066169            1     ensembl
##         [5]        NA      <NA> ENSRNOG00000066169            1     ensembl
##         ...       ...       ...                ...          ...         ...
##   [1284425]        NA      <NA> ENSRNOG00000071103            1     ensembl
##   [1284426]        NA      <NA> ENSRNOG00000071103            1     ensembl
##   [1284427]        NA      <NA> ENSRNOG00000071103            1     ensembl
##   [1284428]        NA      <NA> ENSRNOG00000071103            1     ensembl
##   [1284429]        NA      <NA> ENSRNOG00000071103            1     ensembl
##               gene_biotype      transcript_id transcript_version
##                <character>        <character>        <character>
##         [1] protein_coding               <NA>               <NA>
##         [2] protein_coding ENSRNOT00000101581                  1
##         [3] protein_coding ENSRNOT00000101581                  1
##         [4] protein_coding ENSRNOT00000101581                  1
##         [5] protein_coding ENSRNOT00000101581                  1
##         ...            ...                ...                ...
##   [1284425]         lncRNA ENSRNOT00000106124                  1
##   [1284426]         lncRNA ENSRNOT00000096769                  1
##   [1284427]         lncRNA ENSRNOT00000096769                  1
##   [1284428]         lncRNA ENSRNOT00000096769                  1
##   [1284429]         lncRNA ENSRNOT00000096769                  1
##             transcript_source transcript_biotype               tag exon_number
##                   <character>        <character>       <character> <character>
##         [1]              <NA>               <NA>              <NA>        <NA>
##         [2]           ensembl     protein_coding Ensembl_canonical        <NA>
##         [3]           ensembl     protein_coding Ensembl_canonical           1
##         [4]           ensembl     protein_coding Ensembl_canonical           1
##         [5]           ensembl     protein_coding Ensembl_canonical           2
##         ...               ...                ...               ...         ...
##   [1284425]           ensembl             lncRNA Ensembl_canonical           2
##   [1284426]           ensembl             lncRNA              <NA>        <NA>
##   [1284427]           ensembl             lncRNA              <NA>           1
##   [1284428]           ensembl             lncRNA              <NA>           2
##   [1284429]           ensembl             lncRNA              <NA>           3
##                        exon_id exon_version         protein_id protein_version
##                    <character>  <character>        <character>     <character>
##         [1]               <NA>         <NA>               <NA>            <NA>
##         [2]               <NA>         <NA>               <NA>            <NA>
##         [3] ENSRNOE00000618632            1               <NA>            <NA>
##         [4]               <NA>         <NA> ENSRNOP00000083062               1
##         [5] ENSRNOE00000610554            1               <NA>            <NA>
##         ...                ...          ...                ...             ...
##   [1284425] ENSRNOE00000593698            1               <NA>            <NA>
##   [1284426]               <NA>         <NA>               <NA>            <NA>
##   [1284427] ENSRNOE00000616959            1               <NA>            <NA>
##   [1284428] ENSRNOE00000606447            1               <NA>            <NA>
##   [1284429] ENSRNOE00000543065            2               <NA>            <NA>
##               gene_name transcript_name projection_parent_transcript
##             <character>     <character>                  <character>
##         [1]        <NA>            <NA>                         <NA>
##         [2]        <NA>            <NA>                         <NA>
##         [3]        <NA>            <NA>                         <NA>
##         [4]        <NA>            <NA>                         <NA>
##         [5]        <NA>            <NA>                         <NA>
##         ...         ...             ...                          ...
##   [1284425]        <NA>            <NA>                         <NA>
##   [1284426]        <NA>            <NA>                         <NA>
##   [1284427]        <NA>            <NA>                         <NA>
##   [1284428]        <NA>            <NA>                         <NA>
##   [1284429]        <NA>            <NA>                         <NA>
##   -------
##   seqinfo: 48 sequences from an unspecified genome; no seqlengths
```

gm is a GRanges Object. We want to extract only rows with type ==
"gene". (introns are not interesting)


``` r
gm <- gm[gm$type == "gene"] 
```

keep only colums which are not NA. With mcols you access the data.frame
that is next to the ranges columns.


``` r
# this part somehow performs data munging of the genome annotation
# hard to understand since the genome annotation has a complicated format
# mainly: this code keeps the important stuff and deletes the rest from the genome 
# annotation for rattus norvegicus
meta.keep <- apply(mcols(gm), 2, function(x) sum(!is.na(x)) > 0) 

mcols(gm) <- mcols(gm)[, meta.keep] # subset colnames of mcols of gr
```

add rownames to GRanges Object.


``` r
names(gm) <- gm$gene_id
```

subset to genes present in both, count table and annotation


``` r
gm <- gm[gm$gene_id %in% rownames(df.c)] 
```

match gm and df_ct by EnsembleId


``` r
df.c <- df.c[match(gm$gene_id, rownames(df.c)),]
identical(gm$gene_id, rownames(df.c)) # check if order in EnsembleId is 
```

```
## [1] TRUE
```

``` r
# this should output TRUE
```

# make SummarizedExperiment Object and save as RDS file


``` r
e <- SummarizedExperiment(rowRanges = gm, assays = list("counts" = m.c), 
                          colData = df.m )
e
```

```
## class: RangedSummarizedExperiment 
## dim: 30560 15 
## metadata(0):
## assays(1): counts
## rownames(30560): ENSRNOG00000066169 ENSRNOG00000070168 ...
##   ENSRNOG00000069230 ENSRNOG00000071103
## rowData names(7): source type ... gene_biotype gene_name
## colnames(15): HC1 HC2 ... OA2066 OA2272
## colData names(1): geno
```

``` r
colnames(e)
```

```
##  [1] "HC1"    "HC2"    "HC3"    "DM1"    "DM2"    "DM3"    "GMSM1"  "GMSM2" 
##  [9] "GMSM3"  "OA1789" "OA1887" "OA1966" "OA2022" "OA2066" "OA2272"
```

``` r
colData(e)
```

```
## DataFrame with 15 rows and 1 column
##            geno
##        <factor>
## HC1          HC
## HC2          HC
## HC3          HC
## DM1          DM
## DM2          DM
## ...         ...
## OA1887       OA
## OA1966       OA
## OA2022       OA
## OA2066       OA
## OA2272       OA
```

``` r
saveRDS(e, "e_summExp.rds")
```

# make Deseq2 Object


``` r
dds <- DESeqDataSet(e, design = ~ geno)
dds
```

```
## class: DESeqDataSet 
## dim: 30560 15 
## metadata(1): version
## assays(1): counts
## rownames(30560): ENSRNOG00000066169 ENSRNOG00000070168 ...
##   ENSRNOG00000069230 ENSRNOG00000071103
## rowData names(7): source type ... gene_biotype gene_name
## colnames(15): HC1 HC2 ... OA2066 OA2272
## colData names(1): geno
```

# remove genes which have zero counts for all sample


``` r
dds <- dds[ rowSums( counts(dds) ) > 0 , ]
class(dds)
```

```
## [1] "DESeqDataSet"
## attr(,"package")
## [1] "DESeq2"
```

# run DESeq : Modeling counts against genotype


``` r
dds <- DESeq(dds)
class(dds)
```

```
## [1] "DESeqDataSet"
## attr(,"package")
## [1] "DESeq2"
```

this function performed the following steps:

-   estimate size factors

-   estimate gene-wise dispersion

-   fit curve to gene-wise disperation estimates

-   shrink gene-wise dispersion estimates

-   GLM fit for each gene

take a detailed look at each of this steps
[linkToHbc](https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html)

## some additional comments on estimating size factors

Check the size factors


``` r
sizeFactors(dds)
```

```
##       HC1       HC2       HC3       DM1       DM2       DM3     GMSM1     GMSM2 
## 0.9064680 1.2867470 0.7876660 1.3637949 1.0412807 1.1190991 0.5547837 1.2662609 
##     GMSM3    OA1789    OA1887    OA1966    OA2022    OA2066    OA2272 
## 0.9696542 0.8996994 0.9251896 0.9823004 0.8094461 1.1683839 1.3132017
```

Take a look at the total number of reads for each sample


``` r
colSums(counts(dds))
```

```
##      HC1      HC2      HC3      DM1      DM2      DM3    GMSM1    GMSM2 
## 14605237 21336704 13266649 24548865 18966108 17492986  9465863 20863915 
##    GMSM3   OA1789   OA1887   OA1966   OA2022   OA2066   OA2272 
## 16340637 14970134 15824119 16784689 13369590 19705148 21829214
```

How do the numbers correlate with the size factor?

-   the smaller the total number of reads, the smaller the library sizes

Now take a look at the total depth after normalization using:


``` r
colSums(counts(dds, normalized=T))
```

```
##      HC1      HC2      HC3      DM1      DM2      DM3    GMSM1    GMSM2 
## 16112248 16581895 16842988 18000408 18214213 15631311 17062259 16476790 
##    GMSM3   OA1789   OA1887   OA1966   OA2022   OA2066   OA2272 
## 16852025 16639040 17103650 17087125 16516961 16865302 16622895
```

Normalized reads are obtained by dividing the total reads by the size
factors


``` r
colSums(counts(dds))/sizeFactors(dds)
```

```
##      HC1      HC2      HC3      DM1      DM2      DM3    GMSM1    GMSM2 
## 16112248 16581895 16842988 18000408 18214213 15631311 17062259 16476790 
##    GMSM3   OA1789   OA1887   OA1966   OA2022   OA2066   OA2272 
## 16852025 16639040 17103650 17087125 16516961 16865302 16622895
```

## some additional comments on estimating gene-wise dispersion

try to plot mean count vs. variance counts
[link](https://avikarn.com/2020-07-02-RNAseq_DeSeq2/)


``` r
# Calculating mean for each gene
mean_readCounts <- apply(m.c[,1:3], 1, mean) #change [,1:3] to columns you want to pick

# Calculating variance for each gene
var_readCounts <- apply(m.c[,1:3], 1, var) #here too, in this exp. only for GMSM1-GMSM3

df <- data.frame(mean_readCounts, var_readCounts)

ggplot(df) +
  geom_point(aes(x=mean_readCounts, y= var_readCounts)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene") +
  labs(title = "DESeq2 model - Dispersion")
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/modeldispersion-1.pdf)<!-- -->

## some additional comments on fitting a curve to gene-wise disperation estimates


``` r
plotDispEsts(dds)
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/genewisedispersion-1.pdf)<!-- -->

# Extract Results

The model coefficients returned can be obtained like this:


``` r
resultsNames(dds)
```

```
## [1] "Intercept"     "geno_DM_vs_HC" "geno_GM_vs_HC" "geno_OA_vs_HC"
```

The result function is used to extract the estimated coefficients. We
need to pass the name of the comparison from which we want to get the
results, and the significance level (alpha defines how the "independent
filtering" of lowly expressed genes will work).


``` r
res.DM.vs.HC <- results(dds, name = "geno_DM_vs_HC", alpha = 0.05)
```


``` r
res.GM.vs.HC <- results(dds, name = "geno_GM_vs_HC", alpha = 0.05)
```


``` r
res.OA.vs.HC <- results(dds, name = "geno_OA_vs_HC", alpha = 0.05)
```

remove values where padjust or pvalue == NA


``` r
res.DM.vs.HC.noNa <- res.DM.vs.HC[!is.na(res.DM.vs.HC$pvalue) | !is.na(res.DM.vs.HC$padj), ]
```


``` r
res.GM.vs.HC.noNa <- res.GM.vs.HC[!is.na(res.GM.vs.HC$pvalue) | !is.na(res.GM.vs.HC$padj), ]
```


``` r
res.OA.vs.HC.noNa <- res.OA.vs.HC[!is.na(res.OA.vs.HC$pvalue) | !is.na(res.OA.vs.HC$padj), ]
```

order values by p values


``` r
res.DM.vs.HC.noNaOrd <- res.DM.vs.HC.noNa[order(res.DM.vs.HC.noNa$pvalue), ]
head(res.DM.vs.HC.noNaOrd, n=4)
```

```
## log2 fold change (MLE): geno DM vs HC 
## Wald test p-value: geno DM vs HC 
## DataFrame with 4 rows and 6 columns
##                     baseMean log2FoldChange     lfcSE      stat      pvalue
##                    <numeric>      <numeric> <numeric> <numeric>   <numeric>
## ENSRNOG00000036859   2626.48       -1.55923  0.129904  -12.0029 3.43084e-33
## ENSRNOG00000010645   1011.26        4.60585  0.387109   11.8981 1.21123e-32
## ENSRNOG00000017962   2698.28        1.48782  0.126075   11.8010 3.85534e-32
## ENSRNOG00000006019   7215.61       -1.83594  0.164275  -11.1760 5.34220e-29
##                           padj
##                      <numeric>
## ENSRNOG00000036859 4.74210e-29
## ENSRNOG00000010645 8.37081e-29
## ENSRNOG00000017962 1.77628e-28
## ENSRNOG00000006019 1.84600e-25
```


``` r
res.GM.vs.HC.noNaOrd <- res.GM.vs.HC.noNa[order(res.GM.vs.HC.noNa$pvalue), ]
head(res.GM.vs.HC.noNaOrd, n=4)
```

```
## log2 fold change (MLE): geno GM vs HC 
## Wald test p-value: geno GM vs HC 
## DataFrame with 4 rows and 6 columns
##                     baseMean log2FoldChange     lfcSE      stat      pvalue
##                    <numeric>      <numeric> <numeric> <numeric>   <numeric>
## ENSRNOG00000032832 6612.7746        4.62445   1.06006   4.36244 1.28617e-05
## ENSRNOG00000037931  203.4574        5.98096   1.43510   4.16763 3.07778e-05
## ENSRNOG00000032609   13.2718        3.85027   1.09740   3.50854 4.50580e-04
## ENSRNOG00000011702   14.3770        3.55510   1.04045   3.41688 6.33438e-04
##                         padj
##                    <numeric>
## ENSRNOG00000032832  0.250687
## ENSRNOG00000037931  0.299945
## ENSRNOG00000032609  0.999519
## ENSRNOG00000011702  0.999519
```


``` r
res.OA.vs.HC.noNaOrd <- res.OA.vs.HC.noNa[order(res.OA.vs.HC.noNa$pvalue), ]
head(res.OA.vs.HC.noNaOrd, n=4)
```

```
## log2 fold change (MLE): geno OA vs HC 
## Wald test p-value: geno OA vs HC 
## DataFrame with 4 rows and 6 columns
##                     baseMean log2FoldChange     lfcSE      stat      pvalue
##                    <numeric>      <numeric> <numeric> <numeric>   <numeric>
## ENSRNOG00000007720   28.8254        3.34008  1.009654   3.30814 0.000939172
## ENSRNOG00000032609   13.2718        3.28207  1.001882   3.27590 0.001053247
## ENSRNOG00000003585    3.2015        4.68686  1.468079   3.19251 0.001410419
## ENSRNOG00000020295   24.9110        1.09631  0.348481   3.14597 0.001655397
##                         padj
##                    <numeric>
## ENSRNOG00000007720  0.999759
## ENSRNOG00000032609  0.999759
## ENSRNOG00000003585  0.999759
## ENSRNOG00000020295  0.999759
```

show nr. sig. genes


``` r
summary(res.DM.vs.HC)
```

```
## 
## out of 19494 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 1780, 9.1%
## LFC < 0 (down)     : 1403, 7.2%
## outliers [1]       : 3, 0.015%
## low counts [2]     : 5669, 29%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```


``` r
summary(res.GM.vs.HC)
```

```
## 
## out of 19494 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 0, 0%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 3, 0.015%
## low counts [2]     : 0, 0%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```


``` r
summary(res.OA.vs.HC)
```

```
## 
## out of 19494 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 0, 0%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 3, 0.015%
## low counts [2]     : 0, 0%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

# add shrinked log2FC to result tables

for DM vs HC


``` r
res.DM.vs.HC.lfc <- lfcShrink(dds, coef = "geno_DM_vs_HC",  type="apeglm")
dge <- cbind(res.DM.vs.HC.lfc, 
                          "stat"=res.DM.vs.HC$stat, 
             "lfcSE_not_shrunk" = res.DM.vs.HC$lfcSE, 
             "log2FoldChange_not_shrunk" = res.DM.vs.HC$log2FoldChange)
df_res.DM.vs.HC <- as.data.frame(dge)
df_res.DM.vs.HC$symbol <- rowData(dds)$gene_name
df_res.DM.vs.HC <- df_res.DM.vs.HC[, c("symbol", head(names(df_res.DM.vs.HC), n=-1))]
# this will keep the ENSEMBLE IDs in the excel results table file
df_res.DM.vs.HC$ENSEMBLEID <- rownames(df_res.DM.vs.HC)
# order by pvalue
df_res.DM.vs.HC <- df_res.DM.vs.HC[order(df_res.DM.vs.HC$pvalue), ]
# remove NA from pvalue and padj
df_res.DM.vs.HC.noNa <- df_res.DM.vs.HC[!is.na(df_res.DM.vs.HC$pvalue) & !is.na(df_res.DM.vs.HC$padj), ]
dim(df_res.DM.vs.HC)
```

```
## [1] 19494    10
```

``` r
dim(df_res.DM.vs.HC.noNa)
```

```
## [1] 14200    10
```

``` r
sum(is.na(df_res.DM.vs.HC.noNa$pvalue))
```

```
## [1] 0
```

``` r
sum(duplicated(df_res.DM.vs.HC.noNa$pvalue))
```

```
## [1] 0
```


``` r
res.GM.vs.HC.lfc <- lfcShrink(dds, coef = "geno_GM_vs_HC",  type="apeglm")
dge <- cbind(res.GM.vs.HC.lfc, 
                          "stat"=res.GM.vs.HC$stat, 
             "lfcSE_not_shrunk" = res.GM.vs.HC$lfcSE, 
             "log2FoldChange_not_shrunk" = res.GM.vs.HC$log2FoldChange)
df_res.GM.vs.HC <- as.data.frame(dge)
df_res.GM.vs.HC$symbol <- rowData(dds)$gene_name
df_res.GM.vs.HC <- df_res.GM.vs.HC[, c("symbol", head(names(df_res.GM.vs.HC), n=-1))]
# this will keep the ENSEMBLE IDs in the excel results table file
df_res.GM.vs.HC$ENSEMBLEID <- rownames(df_res.GM.vs.HC)
# order by pvalue
df_res.GM.vs.HC <- df_res.GM.vs.HC[order(df_res.GM.vs.HC$pvalue), ]
# remove NA from pvalue and padj
df_res.GM.vs.HC.noNa <- df_res.GM.vs.HC[!is.na(df_res.GM.vs.HC$pvalue) & !is.na(df_res.GM.vs.HC$padj), ]
dim(df_res.GM.vs.HC)
```

```
## [1] 19494    10
```

``` r
dim(df_res.GM.vs.HC.noNa)
```

```
## [1] 19491    10
```

``` r
sum(is.na(df_res.GM.vs.HC.noNa$pvalue))
```

```
## [1] 0
```

``` r
sum(duplicated(df_res.GM.vs.HC.noNa$pvalue))
```

```
## [1] 2434
```

Unclear how the duplicated stuff happens, literature research shows that
it is common in RNA Seq Data that when there is no significant
difference lots of adjusted p-values can be identical (therefore the
duplicated values in this case 2443)


``` r
res.OA.vs.HC.lfc <- lfcShrink(dds, coef = "geno_OA_vs_HC",  type="apeglm")
dge <- cbind(res.OA.vs.HC.lfc, 
                          "stat"=res.OA.vs.HC$stat, 
             "lfcSE_not_shrunk" = res.OA.vs.HC$lfcSE, 
             "log2FoldChange_not_shrunk" = res.OA.vs.HC$log2FoldChange)
df_res.OA.vs.HC <- as.data.frame(dge)
df_res.OA.vs.HC$symbol <- rowData(dds)$gene_name
df_res.OA.vs.HC <- df_res.OA.vs.HC[, c("symbol", head(names(df_res.OA.vs.HC), n=-1))]
# this will keep the ENSEMBLE IDs in the excel results table file
df_res.OA.vs.HC$ENSEMBLEID <- rownames(df_res.OA.vs.HC)
# order by pvalue
df_res.OA.vs.HC <- df_res.OA.vs.HC[order(df_res.OA.vs.HC$pvalue), ]
# remove NA from pvalue and padj
df_res.OA.vs.HC.noNa <- df_res.OA.vs.HC[!is.na(df_res.OA.vs.HC$pvalue) & !is.na(df_res.OA.vs.HC$padj), ]
dim(df_res.OA.vs.HC)
```

```
## [1] 19494    10
```

``` r
dim(df_res.OA.vs.HC.noNa)
```

```
## [1] 19491    10
```

``` r
sum(is.na(df_res.OA.vs.HC.noNa$pvalue))
```

```
## [1] 0
```

``` r
sum(duplicated(df_res.OA.vs.HC.noNa$pvalue))
```

```
## [1] 2435
```

Plot histograms to visualize the p value distribution


``` r
# added this line to visualize the p value distribution
# literature research shows that it is quite normal for the FDR
# to be the same for comparisons that are not significant this explains the many 
# duplicates
hist(df_res.DM.vs.HC$pvalue)
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-39-1.png)<!-- -->

``` r
hist(df_res.GM.vs.HC$pvalue)
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-39-2.png)<!-- -->

``` r
hist(df_res.OA.vs.HC$pvalue)
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-39-3.png)<!-- -->

ignore next code block, it was commented and kept for bugfixing



# write out all and subsetted results into excel file

list with all (unfiltered and filtered) results


``` r
l_df.res <- list(DM.vs.HC=df_res.DM.vs.HC,
                 GM.vs.HC=df_res.GM.vs.HC,
                 OA.vs.HC=df_res.OA.vs.HC)
l_df.res.noNa <- list(DM.vs.HC.noNa=df_res.DM.vs.HC.noNa,
                      GM.vs.HC.noNa=df_res.GM.vs.HC.noNa,
                      OA.vs.HC.noNa=df_res.OA.vs.HC.noNa)
```

list with genes with sig. padj \< 0.05, separated by up and down


``` r
# noNa.sigPadjUp
l_df.res.noNa.sigPadjUp <- list()
v_nam <- c("DM.vs.HC.noNa.sigPadjUp",
           "GM.vs.HC.noNa.sigPadjUp",
           "OA.vs.HC.noNa.sigPadjUp")
for(i in seq_along(l_df.res.noNa)) {
  l_df.res.noNa.sigPadjUp[[v_nam[i]]] <- l_df.res.noNa[[i]][l_df.res.noNa[[i]]$padj < 0.05 & 
                                                              l_df.res.noNa[[i]]$log2FoldChange > 0, ]
}

# noNa.sigPadjDown
l_df.res.noNa.sigPadjDo <- list()
v_nam <- c("DM.vs.HC.noNa.sigPadjDo",
           "GM.vs.HC.noNa.sigPadjDo",
           "OA.vs.HC.noNa.sigPadjDo")
for(i in seq_along(l_df.res.noNa)) {
  l_df.res.noNa.sigPadjDo[[v_nam[i]]] <- l_df.res.noNa[[i]][l_df.res.noNa[[i]]$padj < 0.05 & 
                                                              l_df.res.noNa[[i]]$log2FoldChange < 0, ]
}
```

list with genes with sig. padj \< 0.05 & foldChange \> 2, separated by
up and down


``` r
# noNa.sigPadjUp
l_df.res.noNa.sigPadjFcUp <- list()
v_nam <- c("DM.vs.HC.noNa.sigPadjFcUp",
           "GM.vs.HC.noNa.sigPadjFcUp",
           "OA.vs.HC.noNa.sigPadjFcUp")
for(i in seq_along(l_df.res.noNa)) {
  l_df.res.noNa.sigPadjFcUp[[v_nam[i]]] <- l_df.res.noNa[[i]][l_df.res.noNa[[i]]$padj < 0.05 & 
                                                                l_df.res.noNa[[i]]$log2FoldChange > 1, ]
}

l_df.res.noNa.sigPadjFcDo <- list()
v_nam <- c("DM.vs.HC.noNa.sigPadjFcDo",
           "GM.vs.HC.noNa.sigPadjFcDo",
           "OA.vs.HC.noNa.sigPadjFcDo")
for(i in seq_along(l_df.res.noNa)) {
  l_df.res.noNa.sigPadjFcDo[[v_nam[i]]] <- l_df.res.noNa[[i]][l_df.res.noNa[[i]]$padj < 0.05 & 
                                                                l_df.res.noNa[[i]]$log2FoldChange < -1, ]
}
```

concatenate all lists, and write out


``` r
l_df.res.all <- c(l_df.res, l_df.res.noNa, 
                  l_df.res.noNa.sigPadjUp, 
                  l_df.res.noNa.sigPadjDo, 
                  l_df.res.noNa.sigPadjFcUp, 
                  l_df.res.noNa.sigPadjFcDo)
write.xlsx(l_df.res.all, "l_df.res.all_vs_HC.xlsx")
```

save RDS to read in easily for later pathway analysis


``` r
saveRDS(l_df.res.all, "l_df.res.all_vs_HC.rds")
```

# visualize Results

## plot MA Plot

for DM vs GM


``` r
# par(mfrow=c(2,2))

plotMA(res.DM.vs.HC, ylim=c(-2,2), alpha=0.05)
plotMA(res.DM.vs.HC.lfc, ylim=c(-2,2), alpha=0.05)
# dev.off()
```

<img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-44-1.png" width="50%" /><img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-44-2.png" width="50%" />


``` r
# par(mfrow=c(2,2))

plotMA(res.GM.vs.HC, ylim=c(-2,2), alpha=0.05)
plotMA(res.GM.vs.HC.lfc, ylim=c(-2,2), alpha=0.05)
# dev.off()
```

<img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-45-1.png" width="50%" /><img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-45-2.png" width="50%" />


``` r
# par(mfrow=c(2,2))

plotMA(res.OA.vs.HC, ylim=c(-2,2), alpha=0.05)
plotMA(res.OA.vs.HC.lfc, ylim=c(-2,2), alpha=0.05)
# dev.off()
```

<img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-46-1.png" width="50%" /><img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-46-2.png" width="50%" />



## plot normalized counts


``` r
# par(mfrow=c(1,2))
plotCounts(dds, gene=which.min(res.DM.vs.HC$padj), intgroup="geno") #this gene has the smallest p value
plotCounts(dds, gene=which.min(res.GM.vs.HC$padj), intgroup="geno")

plotCounts(dds, gene=which.min(res.OA.vs.HC$padj), intgroup="geno")

# dev.off()
```

<img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-47-1.png" width="50%" /><img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-47-2.png" width="50%" /><img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-47-3.png" width="50%" />

## plot Volcano plots

for DM vs GM


``` r
#reset par
#par(mfrow=c(1,2))
# Make a basic volcano plot
with(res.DM.vs.HC, plot(log2FoldChange, -log10(pvalue), pch=20, main="DM.vs.HC", xlim=c(-10,10)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res.DM.vs.HC, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.DM.vs.HC, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#dev.off()
```

<embed src="RNASeqDESeq2LWData_vs_HC_files/figure-html/DMvsGM-1.pdf" width="50%" type="application/pdf" />


``` r
#reset par
#par(mfrow=c(1,2))
# Make a basic volcano plot
with(res.GM.vs.HC, plot(log2FoldChange, -log10(pvalue), pch=20, main="GM.vs.HC", xlim=c(-10,10)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res.GM.vs.HC, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.GM.vs.HC, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#dev.off()
```

<embed src="RNASeqDESeq2LWData_vs_HC_files/figure-html/HCvsGM-1.pdf" width="50%" type="application/pdf" />


``` r
#reset par
# par(mfrow=c(1,2))
# Make a basic volcano plot
with(res.OA.vs.HC, plot(log2FoldChange, -log10(pvalue), pch=20, main="OA.vs.HC", xlim=c(-10,10)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res.OA.vs.HC, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.OA.vs.HC, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#dev.off()
```

<embed src="RNASeqDESeq2LWData_vs_HC_files/figure-html/OAvsGM-1.pdf" width="50%" type="application/pdf" />


## perform variance stabilizing transformation


``` r
vst.b.f <- vst(dds, blind=FALSE)
vst.b.t <- vst(dds, blind=TRUE)
vst.b.f
```

```
## class: DESeqTransform 
## dim: 19494 15 
## metadata(1): version
## assays(1): ''
## rownames(19494): ENSRNOG00000070901 ENSRNOG00000018029 ...
##   ENSRNOG00000062442 ENSRNOG00000071103
## rowData names(37): source type ... maxCooks dispFit
## colnames(15): HC1 HC2 ... OA2066 OA2272
## colData names(2): geno sizeFactor
```

``` r
head(assay(vst.b.t), 3)
```

```
##                          HC1       HC2       HC3       DM1       DM2       DM3
## ENSRNOG00000070901 11.121336 10.764381 11.030155 11.556787 11.902706 11.000388
## ENSRNOG00000018029  6.478058  6.148604  6.162030  6.400313  6.109494  6.809121
## ENSRNOG00000031391  5.269703  5.474298  5.269703  5.269703  5.269703  5.489061
##                        GMSM1     GMSM2     GMSM3    OA1789    OA1887    OA1966
## ENSRNOG00000070901 10.853322 10.926054 10.780855 10.881857 10.830393 10.821789
## ENSRNOG00000018029  6.024920  6.373194  6.350632  6.546367  6.024700  5.885085
## ENSRNOG00000031391  5.269703  5.269703  5.269703  5.269703  5.269703  5.269703
##                       OA2022    OA2066    OA2272
## ENSRNOG00000070901 10.912165 10.922214 10.850837
## ENSRNOG00000018029  6.340304  6.360151  6.204006
## ENSRNOG00000031391  5.527533  5.269703  5.472230
```

Variance stabilization shrinks standard deviation of lowly expressed
genes.\
It is similar, but not identical, to log transformation.\
this gives log2(n + 1)


``` r
ntd <- normTransform(dds)
```

now plot distribution of log2 transformed, vst transformed
(transformatin across all samples), vst transformed (transformation
within each group/genotype)


``` r
# par(mfrow=c(1,4))
meanSdPlot(assay(ntd))
meanSdPlot(assay(vst.b.t))
meanSdPlot(assay(vst.b.f))
# dev.off()
```

<img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-50-1.png" width="33%" /><img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-50-2.png" width="33%" /><img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-50-3.png" width="33%" />

## plot heatmap of sig. genes


``` r
m_vst.b.f <- assay(vst.b.f)
```

### heatmap DM vs HC

plot sig. genes based on padj \< 0.05 cutoff


``` r
pheatmap(m_vst.b.f[rownames(m_vst.b.f) %in% 
                     c(rownames(l_df.res.all$DM.vs.HC.noNa.sigPadjUp),
                       rownames(l_df.res.all$DM.vs.HC.noNa.sigPadjDo)), 1:6], 
         scale = "row", 
         cluster_cols = TRUE, 
         cluster_rows = TRUE, 
         show_rownames = FALSE)
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-52-1.png)<!-- -->

there are similar numbers of up and down, when cutoff is only padj \<
0.05

plot sig. genes based on padj \< 0.05 and foldChange \> 2 cutoff


``` r
pheatmap(m_vst.b.f[rownames(m_vst.b.f) %in% 
                     c(rownames(l_df.res.all$DM.vs.HC.noNa.sigPadjFcUp), 
                     rownames(l_df.res.all$DM.vs.HC.noNa.sigPadjFcDo)), 1:6], 
         scale = "row", 
         cluster_cols = TRUE, 
         cluster_rows = TRUE, 
         show_rownames = FALSE)
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-53-1.png)<!-- -->

with additional log2fc cutoff we see that there is a stronger effect on
down-regulated genes, than on up-regulated genes



``` r
'pheatmap(m_vst.b.f[rownames(m_vst.b.f) %in% 
                     c(rownames(l_df.res.all$HC.vs.GM.noNa.sigPadjUp),
                       rownames(l_df.res.all$HC.vs.GM.noNa.sigPadjDo)), 1:6], 
         scale = "row", 
         cluster_cols = TRUE, 
         cluster_rows = TRUE, 
         show_rownames = FALSE)'
```

```
## [1] "pheatmap(m_vst.b.f[rownames(m_vst.b.f) %in% \n                     c(rownames(l_df.res.all$HC.vs.GM.noNa.sigPadjUp),\n                       rownames(l_df.res.all$HC.vs.GM.noNa.sigPadjDo)), 1:6], \n         scale = \"row\", \n         cluster_cols = TRUE, \n         cluster_rows = TRUE, \n         show_rownames = FALSE)"
```

### heatmap GM vs HC

plot sig. genes based on padj \< 0.05 cutoff


``` r
'pheatmap(m_vst.b.f[rownames(m_vst.b.f) %in% 
                     c(rownames(l_df.res.all$GM.vs.HC.noNa),
                       rownames(l_df.res.all$GM.vs.HC.noNa)), 1:6], 
         scale = "row", 
         cluster_cols = TRUE, 
         cluster_rows = TRUE, 
         show_rownames = FALSE)'
```

```
## [1] "pheatmap(m_vst.b.f[rownames(m_vst.b.f) %in% \n                     c(rownames(l_df.res.all$GM.vs.HC.noNa),\n                       rownames(l_df.res.all$GM.vs.HC.noNa)), 1:6], \n         scale = \"row\", \n         cluster_cols = TRUE, \n         cluster_rows = TRUE, \n         show_rownames = FALSE)"
```

there are similar numbers of up and down, when cutoff is only padj \<
0.05


plot sig. genes based on padj \< 0.05 and foldChange \> 2 cutoff




# 

# Do some Quality Control

## QC: Heatmap of the sample-to-sample distances

calculate distances of samples to each other (distance matrix)


``` r
sampleDists <- dist(t(assay(vst.b.t)))
sampleDists
```

```
##             HC1      HC2      HC3      DM1      DM2      DM3    GMSM1    GMSM2
## HC2    29.45164                                                               
## HC3    49.89047 43.84366                                                      
## DM1    78.84607 78.15695 83.72138                                             
## DM2    79.30831 80.22016 79.79253 32.86280                                    
## DM3    59.37176 55.64602 70.59555 78.20150 80.02698                           
## GMSM1  45.91974 37.45318 29.87584 80.76991 79.96818 67.71866                  
## GMSM2  32.99149 30.11617 50.99798 77.40696 78.52518 54.99919 47.62008         
## GMSM3  52.95641 41.76522 31.79093 82.83726 80.45531 62.90778 34.22337 44.48469
## OA1789 33.48947 29.15466 41.66672 75.36939 77.81974 60.19484 37.10119 40.45320
## OA1887 53.99363 43.22294 31.31188 83.40256 80.12778 66.04946 36.97334 45.34075
## OA1966 46.89929 36.42437 25.71718 81.01046 78.54540 65.75861 29.04665 45.42264
## OA2022 47.24226 41.83682 54.41366 91.70697 90.78323 62.28314 53.26701 38.32364
## OA2066 48.04007 38.61223 26.12228 78.44320 75.89547 65.13916 31.70554 42.76390
## OA2272 39.28097 29.56503 30.14697 78.12347 77.18184 61.24884 32.39462 33.95288
##           GMSM3   OA1789   OA1887   OA1966   OA2022   OA2066
## HC2                                                         
## HC3                                                         
## DM1                                                         
## DM2                                                         
## DM3                                                         
## GMSM1                                                       
## GMSM2                                                       
## GMSM3                                                       
## OA1789 45.90251                                             
## OA1887 23.19761 47.11779                                    
## OA1966 27.30167 38.79195 25.68551                           
## OA2022 45.84229 55.23234 46.26066 47.59006                  
## OA2066 25.08032 39.60303 22.59722 22.21587 48.26986         
## OA2272 29.35806 32.90880 28.88959 24.98173 43.10810 23.10234
```

plot heatmap of sample distances


``` r
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vst.b.t$geno
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-57-1.png)<!-- -->

## QC: Principal component plot of the samples


``` r
plotPCA(vst.b.t, intgroup=c("geno"))
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-58-1.png)<!-- -->

## QC: p-value distribution


``` r
# par(mfrow=c(1, 2))
hist(l_df.res.noNa$DM.vs.HC.noNa$pvalue, main = "p-value distribution: DM.vs.HC", xlab = NULL)

# dev.off()
```

<img src="RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-59-1.png" width="50%" />

## QC: distribution of vst transformed values


``` r
l_dens <- list()
for(i in 1:ncol(m_vst.b.f)){
  l_dens[[i]] <- density(m_vst.b.f[, i])
}
v_col <- sample(colors(), size = 8, replace = F)

plot(l_dens[[1]], lwd = 2, col = "red",
     main = "plotDensityOfDifferentSamples", xlab = "")
for(i in 2:ncol(m_vst.b.f)){
  lines(l_dens[[i]], col = i-1, lwd = 2)
}
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-60-1.png)<!-- -->

``` r
# dev.off()
```

# Independent filtering of results

for DM.vs.HC


``` r
plot(res.DM.vs.HC$baseMean+1, -log10(res.DM.vs.HC$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
abline(h=-10*log10(0.05), col="red")
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-61-1.png)<!-- -->

The p-value cutoff set when calling the "results" - function is stored
here:


``` r
metadata(res.DM.vs.HC)$alpha
```

```
## [1] 0.05
```

The thresholds set by independent filtering are stored here:


``` r
metadata(res.DM.vs.HC)$filterThreshold
```

```
## 29.08163% 
##  2.286906
```

The 29% of lowest value will be filtered out from the data.\
All genes with a mean read count \< 2.28 will be filtered out.

Plot number of rejections of the Null hypothesis vs. percentiles of data
filtered out


``` r
plot(metadata(res.DM.vs.HC)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="fractionFilteredOut")
lines(metadata(res.DM.vs.HC)$lo.fit, col="red")
abline(v=metadata(res.DM.vs.HC)$filterTheta)
```

![](RNASeqDESeq2LWData_vs_HC_files/figure-html/unnamed-chunk-64-1.png)<!-- -->

rank genes by their mean of normalized counts, filter out increasing
fraction, so starting from the most lowly expressed genes. removing the
very lowly expressed genes at the beginning increases the number of sig.
genes. At higher expressed genes, filtering out genes results in a
decrease of sig. genes


``` r
save.image("RNASeq_LW_all_vs_GM")
```

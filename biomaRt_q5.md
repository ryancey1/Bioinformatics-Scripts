---
title: "Sample code for Unit 3-4 HW Q4 using biomaRt"
author: "Ryan Yancey"
date: "2/17/2021"
output:
  html_document: 
    toc: yes
    df_print: default
    keep_md: yes
---



The purpose of this document is to demonstrate how to retrieve a TSV file with data mined from BioMart using Bioconductor's "biomaRt" package. This is not a replacement for the method asked for in the homework, but rather a supplement!

First, let's load necessary libraries.


```r
library("biomaRt")
```

## Setting up biomaRt host

Now, let's structure our query to BioMart. First, we need to set up the host we will be querying to. One way is to find an archived host so the results are reproducible using the current mouse assembly. (February 2021: GRCm39 is current mouse assembly).


```r
listEnsemblArchives()
```

```
##              name     date                                url version
## 1  Ensembl GRCh37 Feb 2014          http://grch37.ensembl.org  GRCh37
## 2     Ensembl 103 Feb 2021 http://feb2021.archive.ensembl.org     103
## 3     Ensembl 102 Nov 2020 http://nov2020.archive.ensembl.org     102
## 4     Ensembl 101 Aug 2020 http://aug2020.archive.ensembl.org     101
## 5     Ensembl 100 Apr 2020 http://apr2020.archive.ensembl.org     100
## 6      Ensembl 99 Jan 2020 http://jan2020.archive.ensembl.org      99
## 7      Ensembl 98 Sep 2019 http://sep2019.archive.ensembl.org      98
## 8      Ensembl 97 Jul 2019 http://jul2019.archive.ensembl.org      97
## 9      Ensembl 96 Apr 2019 http://apr2019.archive.ensembl.org      96
## 10     Ensembl 95 Jan 2019 http://jan2019.archive.ensembl.org      95
## 11     Ensembl 94 Oct 2018 http://oct2018.archive.ensembl.org      94
## 12     Ensembl 93 Jul 2018 http://jul2018.archive.ensembl.org      93
## 13     Ensembl 92 Apr 2018 http://apr2018.archive.ensembl.org      92
## 14     Ensembl 91 Dec 2017 http://dec2017.archive.ensembl.org      91
## 15     Ensembl 90 Aug 2017 http://aug2017.archive.ensembl.org      90
## 16     Ensembl 89 May 2017 http://may2017.archive.ensembl.org      89
## 17     Ensembl 88 Mar 2017 http://mar2017.archive.ensembl.org      88
## 18     Ensembl 87 Dec 2016 http://dec2016.archive.ensembl.org      87
## 19     Ensembl 86 Oct 2016 http://oct2016.archive.ensembl.org      86
## 20     Ensembl 85 Jul 2016 http://jul2016.archive.ensembl.org      85
## 21     Ensembl 84 Mar 2016 http://mar2016.archive.ensembl.org      84
## 22     Ensembl 80 May 2015 http://may2015.archive.ensembl.org      80
## 23     Ensembl 77 Oct 2014 http://oct2014.archive.ensembl.org      77
## 24     Ensembl 75 Feb 2014 http://feb2014.archive.ensembl.org      75
## 25     Ensembl 67 May 2012 http://may2012.archive.ensembl.org      67
## 26     Ensembl 54 May 2009 http://may2009.archive.ensembl.org      54
##    current_release
## 1                 
## 2                *
## 3                 
## 4                 
## 5                 
## 6                 
## 7                 
## 8                 
## 9                 
## 10                
## 11                
## 12                
## 13                
## 14                
## 15                
## 16                
## 17                
## 18                
## 19                
## 20                
## 21                
## 22                
## 23                
## 24                
## 25                
## 26
```

```r
listMarts(host = "http://feb2021.archive.ensembl.org")
```

```
##                biomart                version
## 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 103
## 2   ENSEMBL_MART_MOUSE      Mouse strains 103
## 3     ENSEMBL_MART_SNP  Ensembl Variation 103
## 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 103
```

We will set the mart to be the newest archive of "Ensembl Genes 103" by specifying the biomart variable as "ENSEMBL_MART_ENSEMBL" in the `useMart()` function.


```r
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
```

To restrict our search to the **newest** mouse genome, let's search mouse-related data sets available and find the *R-formatted* name of most recent genome assembly (GRCm39).


```r
## [Mm]ouse matches both "Mouse" and "mouse"
searchDatasets(ensembl, pattern = "[Mm]ouse")
```

```
##                      dataset                                      description
## 95      mcaroli_gene_ensembl             Ryukyu mouse genes (CAROLI_EIJ_v1.1)
## 106    mmurinus_gene_ensembl                     Mouse Lemur genes (Mmur_3.0)
## 107   mmusculus_gene_ensembl                             Mouse genes (GRCm39)
## 110     mpahari_gene_ensembl              Shrew mouse genes (PAHARI_EIJ_v1.1)
## 112 mspicilegus_gene_ensembl                     Steppe mouse genes (MUSP714)
## 113    mspretus_gene_ensembl              Algerian mouse genes (SPRET_EiJ_v1)
## 150   pmbairdii_gene_ensembl Northern American deer mouse genes (HU_Pman_2.1)
##             version
## 95  CAROLI_EIJ_v1.1
## 106        Mmur_3.0
## 107          GRCm39
## 110 PAHARI_EIJ_v1.1
## 112         MUSP714
## 113    SPRET_EiJ_v1
## 150     HU_Pman_2.1
```

**Row 107** of the list has the assembly we want. Let's re-format our mart to include this information.


```r
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl",
                      mart = ensembl)
```

## Structuring the `getBM()` query

Now that our connection to the Ensembl mart is set up and structured as we want, we can actually query the data set for genes of interest. 

So, to structure we will need to know three things:

1. What **filters** are available,
2. What **attributes** are available, and
3. What **values** we want to specify 

(See [here](https://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#how-to-build-a-biomart-query) for more info). To find how those are formatted in R, we can use the `searchAttributes()` and `searchFilters()` function.

On the BioMart website, each data search has default attributes associated with the output:

* Ensembl gene stable ID
* Ensembl gene stable ID version
* Ensembl transcript stable ID
* Ensembl transcript stable ID version

We can search the attributes list for ones that match "some number of characters (transcript or gene), followed by ID" to narrow the results.


```r
## search attributes list for defaults
head(searchAttributes(ensembl, pattern="[a-z]+ stable ID"), n = 10)
```

```
##                              name                  description         page
## 1                 ensembl_gene_id               Gene stable ID feature_page
## 2         ensembl_gene_id_version       Gene stable ID version feature_page
## 3           ensembl_transcript_id         Transcript stable ID feature_page
## 4   ensembl_transcript_id_version Transcript stable ID version feature_page
## 5              ensembl_peptide_id            Protein stable ID feature_page
## 6      ensembl_peptide_id_version    Protein stable ID version feature_page
## 7                 ensembl_exon_id               Exon stable ID feature_page
## 177               ensembl_gene_id               Gene stable ID    structure
## 178       ensembl_gene_id_version       Gene stable ID version    structure
## 180         ensembl_transcript_id         Transcript stable ID    structure
```

Luckily enough, the first 4 are what we're searching for. Next, we search for attributes and filters we are personally interested in. For this question, we only want the RefSeq peptide ID, though there are many options available.


```r
## search attributes list for non-defaults
searchAttributes(ensembl, pattern="RefSeq peptide ID")
```

```
##              name       description         page
## 78 refseq_peptide RefSeq peptide ID feature_page
```

Finally, we want to build our gene filtering schema. We want to filter by the chromosome number (11), karyotype band E2 (start:end -- 110433446:122082543), the transcript count (>= 7), and the presence of a RefSeq peptide ID (TRUE).


```r
## search filter list
searchFilters(ensembl, pattern = "chr|start|end|greater_than|peptide")
```

```
##                              name
## 1                 chromosome_name
## 2                           start
## 3                             end
## 4                      band_start
## 5                        band_end
## 7              chromosomal_region
## 35            with_refseq_peptide
## 36  with_refseq_peptide_predicted
## 49             ensembl_peptide_id
## 50     ensembl_peptide_id_version
## 85                 refseq_peptide
## 86       refseq_peptide_predicted
## 145 transcript_count_greater_than
## 241     with_acchrysaetos_homolog
## 311     with_mochrogaster_homolog
## 325       with_bsplendens_homolog
## 330    with_apolyacanthus_homolog
##                                                       description
## 1                                        Chromosome/scaffold name
## 2                                                           Start
## 3                                                             End
## 4                                                      Band Start
## 5                                                        Band End
## 7                          e.g. 1:100:10000:-1, 1:100000:200000:1
## 35                                      With RefSeq peptide ID(s)
## 36                            With RefSeq peptide predicted ID(s)
## 49                 Protein stable ID(s) [e.g. ENSMUSP00000000001]
## 50  Protein stable ID(s) with version [e.g. ENSMUSP00000000001.5]
## 85                       RefSeq peptide ID(s) [e.g. NP_001001130]
## 86             RefSeq peptide predicted ID(s) [e.g. XP_001004117]
## 145                                           Transcript count >=
## 241                                Orthologous Golden eagle Genes
## 311                                Orthologous Prairie vole Genes
## 325                       Orthologous Siamese fighting fish Genes
## 330                               Orthologous Spiny chromis Genes
```

Now, we can finally structure the filters and attributes to restrict our biomaRt query to only those genes we're interested in.


```r
## these are used to search the database and filter out undesired genes
filters <- c("chromosome_name",
             "start",
             "end",
             "transcript_count_greater_than", 
             "with_refseq_peptide")

## these are the values for each filter we specified
values <- list("11",          ## chromosome_name
               "110433446",   ## start
               "122082543",   ## end
               "7",           ## transcript_count_greater_than
               TRUE)          ## with_refseq_peptide

## these are features used to format the output table
attributes <- c("ensembl_gene_id",                ## default
                "ensembl_gene_id_version",        ## default
                "ensembl_transcript_id",          ## default
                "ensembl_transcript_id_version",  ## default
                "refseq_peptide")                 ## non-default
```

Finally, we can query the database. The `values` variable lists that we want genes from Chromosome 11 with transcript counts greater than or equal to 7 that also contain a RefSeq Peptide ID.


```r
## save query to bm object
bm <- getBM(filters = filters,
             values = values,
             attributes = attributes,
             mart = ensembl)
```

Summary statistics on `bm` dataset.


```r
summary(bm)
```

```
##  ensembl_gene_id    ensembl_gene_id_version ensembl_transcript_id
##  Length:176         Length:176              Length:176           
##  Class :character   Class :character        Class :character     
##  Mode  :character   Mode  :character        Mode  :character     
##  ensembl_transcript_id_version refseq_peptide    
##  Length:176                    Length:176        
##  Class :character              Class :character  
##  Mode  :character              Mode  :character
```

```r
nrow(bm)
```

```
## [1] 176
```

```r
summary(unique(bm$ensembl_gene_id))
```

```
##    Length     Class      Mode 
##        65 character character
```

**There are 176 matches to our query, and 65 of them come from unique genes with multiple transcript matches.**

## Exporting results to file

Finally, let's export this to a TSV file. But first let's take a look at a segment of the results.


```r
head(bm)
```

```
##      ensembl_gene_id ensembl_gene_id_version ensembl_transcript_id
## 1 ENSMUSG00000041654   ENSMUSG00000041654.16    ENSMUST00000071539
## 2 ENSMUSG00000041654   ENSMUSG00000041654.16    ENSMUST00000071539
## 3 ENSMUSG00000041654   ENSMUSG00000041654.16    ENSMUST00000042657
## 4 ENSMUSG00000041654   ENSMUSG00000041654.16    ENSMUST00000106633
## 5 ENSMUSG00000041654   ENSMUSG00000041654.16    ENSMUST00000106633
## 6 ENSMUSG00000041654   ENSMUSG00000041654.16    ENSMUST00000106633
##   ensembl_transcript_id_version refseq_peptide
## 1         ENSMUST00000071539.10      NP_081492
## 2         ENSMUST00000071539.10   NP_001349870
## 3         ENSMUST00000042657.16   NP_001349867
## 4         ENSMUST00000106633.10   NP_001159975
## 5         ENSMUST00000106633.10   NP_001349869
## 6         ENSMUST00000106633.10   NP_001349868
```

```r
tail(bm)
```

```
##        ensembl_gene_id ensembl_gene_id_version ensembl_transcript_id
## 171 ENSMUSG00000039294   ENSMUSG00000039294.15    ENSMUST00000106115
## 172 ENSMUSG00000039294   ENSMUSG00000039294.15    ENSMUST00000106115
## 173 ENSMUSG00000039294   ENSMUSG00000039294.15    ENSMUST00000038709
## 174 ENSMUSG00000039294   ENSMUSG00000039294.15    ENSMUST00000038709
## 175 ENSMUSG00000039294   ENSMUSG00000039294.15    ENSMUST00000169393
## 176 ENSMUSG00000039230   ENSMUSG00000039230.15    ENSMUST00000103013
##     ensembl_transcript_id_version refseq_peptide
## 171          ENSMUST00000106115.8   NP_001239478
## 172          ENSMUST00000106115.8   NP_001241664
## 173         ENSMUST00000038709.14      NP_659081
## 174         ENSMUST00000038709.14   NP_001239477
## 175          ENSMUST00000169393.8   NP_001239479
## 176         ENSMUST00000103013.10      NP_084154
```


```r
write.table(bm, file = "mart_export.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```

## Validating our results with the Ensembl results

We can verify our query is correct by performing the same search on the [BioMart](https://www.ensembl.org/biomart/martview/) website. The results after searching with the above query are summarized in the image below.

![][1]

As we can see, the website-based search found 65 genes as well, hinting that our results from above are correct. To double-check, we can load the resulting file (I downloaded it as 'web_mart_export.txt') into our workspace and compare it to our biomaRt-generated results by using the `identical()` function. First, we need to 


```r
bm_web <- read.table("web_mart_export.txt", header = TRUE, sep = "\t")
head(bm_web)
```

```
##       Gene.stable.ID Gene.stable.ID.version Transcript.stable.ID
## 1 ENSMUSG00000041654  ENSMUSG00000041654.16   ENSMUST00000071539
## 2 ENSMUSG00000041654  ENSMUSG00000041654.16   ENSMUST00000071539
## 3 ENSMUSG00000041654  ENSMUSG00000041654.16   ENSMUST00000042657
## 4 ENSMUSG00000041654  ENSMUSG00000041654.16   ENSMUST00000106633
## 5 ENSMUSG00000041654  ENSMUSG00000041654.16   ENSMUST00000106633
## 6 ENSMUSG00000041654  ENSMUSG00000041654.16   ENSMUST00000106633
##   Transcript.stable.ID.version RefSeq.peptide.ID
## 1        ENSMUST00000071539.10         NP_081492
## 2        ENSMUST00000071539.10      NP_001349870
## 3        ENSMUST00000042657.16      NP_001349867
## 4        ENSMUST00000106633.10      NP_001159975
## 5        ENSMUST00000106633.10      NP_001349869
## 6        ENSMUST00000106633.10      NP_001349868
```

```r
## we have to format column names to be consistent with our output's
colnames(bm_web) <- colnames(bm)

## now we can check that they're identical
identical(bm, bm_web)
```

```
## [1] TRUE
```

The two data-frames are identical, so our search was accurate!

## Session info

```r
sessionInfo()
```

```
## R version 4.0.4 (2021-02-15)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] dplyr_1.0.4    biomaRt_2.44.4
## 
## loaded via a namespace (and not attached):
##  [1] progress_1.2.2       tidyselect_1.1.0     xfun_0.21           
##  [4] bslib_0.2.4          purrr_0.3.4          vctrs_0.3.6         
##  [7] generics_0.1.0       htmltools_0.5.1.1    stats4_4.0.4        
## [10] BiocFileCache_1.12.1 yaml_2.2.1           blob_1.2.1          
## [13] XML_3.99-0.5         rlang_0.4.10         jquerylib_0.1.3     
## [16] pillar_1.4.7         glue_1.4.2           DBI_1.1.1           
## [19] rappdirs_0.3.3       BiocGenerics_0.34.0  bit64_4.0.5         
## [22] dbplyr_2.1.0         lifecycle_1.0.0      stringr_1.4.0       
## [25] memoise_2.0.0        evaluate_0.14        Biobase_2.48.0      
## [28] knitr_1.31           IRanges_2.22.2       fastmap_1.1.0       
## [31] curl_4.3             parallel_4.0.4       AnnotationDbi_1.50.3
## [34] Rcpp_1.0.6           openssl_1.4.3        cachem_1.0.4        
## [37] S4Vectors_0.26.1     jsonlite_1.7.2       bit_4.0.4           
## [40] hms_1.0.0            askpass_1.1          digest_0.6.27       
## [43] stringi_1.5.3        tools_4.0.4          magrittr_2.0.1      
## [46] sass_0.3.1           RSQLite_2.2.3        tibble_3.0.6        
## [49] crayon_1.4.1         pkgconfig_2.0.3      ellipsis_0.3.1      
## [52] xml2_1.3.2           prettyunits_1.1.1    assertthat_0.2.1    
## [55] rmarkdown_2.7        httr_1.4.2           R6_2.5.0            
## [58] compiler_4.0.4
```

[1]: biomart_results.png "Biomart Results"

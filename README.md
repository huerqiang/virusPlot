
<!-- README.md is generated from README.Rmd. Please edit that file -->

# virusPlot

<!-- badges: start -->

<!-- badges: end -->

Visualization of virus insertion information

## :writing_hand: Authors

Erqiang Hu: Albert Einstein College of Medicine

Shanye Yin: Albert Einstein College of Medicine

## Availability and implementation

The R package version is available at:
<https://github.com/huerqiang/virusPlot>. The web version of virusPlot
is available at: <http://www.huerqiang.com:58211/app/virusplot>.

## :hammer: Installation

``` r
devtools::install_github("huerqiang/virusPlot")
```

## workflow

<img src="man/figures/overview.png" width="100%" />

## Example

Download virus genome

``` r
library(virusPlot)
#> 
#> 
genome <- get_virus_genom(accession_number = "NC_001526.2",
    email = "13766876214@163.com")
head(genome)
#>                                                                       V1
#> 1             >NC_001526.2 Human papillomavirus type 16, complete genome
#> 2 ACTACAATAATTCATGTATAAAACTAAGGGCGTAACCGAAATCGGTTGAACCGAAACCGGTTAGTATAAA
#> 3 AGCAGACATTTTATGCACCAAAAGAGAACTGCAATGTTTCAGGACCCACAGGAGCGACCCAGAAAGTTAC
#> 4 CACAGTTATGCACAGAGCTGCAAACAACTATACATGATATAATATTAGAATGTGTGTACTGCAAGCAACA
#> 5 GTTACTGCGACGTGAGGTATATGACTTTGCTTTTCGGGATTTATGCATAGTATATAGAGATGGGAATCCA
#> 6 TATGCTGTATGTGATAAATGTTTAAAGTTTTATTCTAAAATTAGTGAGTATAGACATTATTGTTATAGTT
```

Download virus annotation

``` r
gene_features <- get_virus_annotation(accession_number = "NC_001526.2",
    email = "13766876214@163.com")
virus_info <- deal_virus_annotation(gene_features)
virus_info
#>   gene start  end
#> 1   E6    83  559
#> 2   E7   562  858
#> 3   E1   865 2813
#> 4   E2  2755 3852
#> 5   E4  3332 3619
#> 6   E5  3849 4100
#> 7   L2  4236 5657
#> 8   L1  5560 7155
```

``` r
library(virusPlot)
virus_info <- data.frame(
         gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1"),
         start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560),
         end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155))

insert_num <- data.frame(start = seq(1, 7801, 100),
    end = seq(101, 7901, 100),
    num = sample(1:28, 79, replace = TRUE))
circle_virus(virus_info, insert_num)
#> Warning: Removed 8 rows containing missing values or values outside the scale range
#> (`geom_col()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_rect()`).
#> Removed 1 row containing missing values or values outside the scale range
#> (`geom_rect()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text()`).
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

Users need to prepare an input dataframe containing specific columns
that describe viral integration events. These columns include:  
1. Chromosome: The chromosome number where the viral integration
breakpoint is located in the host genome.  
2. Host Position: The position of the breakpoint on the host genome.  
3. Viral Position: The corresponding position on the viral genome.  
4. Read Count: The number of reads supporting the integration event.  
5. Sample ID (optional): The sample identifier, allowing for
visualization of integration events across multiple samples.  
Here we provide the insert_info object as an example:

``` r
data(insert_info)
head(insert_info)
#>    chr  host_loc hpv_loc reads   sample
#> 1 chr8 127886048    4667  1310  sample8
#> 2 chr8 127874880    5296   563  sample8
#> 3 chr8 127727167    2354   415 sample19
#> 4 chr8 127886058    5296   292  sample2
#> 5 chr8 127884965    3224   122 sample19
#> 6 chr8 127882447    2102    88  sample3
```

``` r

# virus_info <- data.frame(
#       gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
#       start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
#       end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))
# HPV16
gene_features <- get_virus_annotation(accession_number = "NC_001526.2",
    email = "13766876214@163.com")
virus_info_NC_001526 <- deal_virus_annotation(gene_features)
strudel_plot(virus_info = virus_info_NC_001526, insert_info, 
             hot_gene = 3)
#> Loading required package: org.Hs.eg.db
#> Loading required package: AnnotationDbi
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: IRanges
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> 
#> Attaching package: 'IRanges'
#> The following object is masked from 'package:grDevices':
#> 
#>     windows
#> 'select()' returned 1:many mapping between keys and columns
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r
strudel_plot(virus_info = virus_info_NC_001526, insert_info, 
             hot_gene = 3, sample_select = c("sample1", "sample2", "sample3", "sample4", "sample5"))
#> 'select()' returned 1:many mapping between keys and columns
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" /> hot
genes of virus

``` r
hot_gene <- get_hot_gene(virus_info = virus_info_NC_001526, insert_info)
insert_plot <- hot_gene_plot(hot_gene)
```

``` r
insert_plot[[2]]
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

``` r
strudel_plot(virus_info = virus_info_NC_001526, insert_info, 
             hot_gene = c( "SAV1", "ZPLD1", "CHMP6", "IRS4"))
#> 'select()' returned 1:many mapping between keys and columns
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

``` r
data(vcf_matrix)
data(col)
data(pdata)
data(cli_colors)
oncoplot(vcf_matrix, varis_color = col, 
    clinical = pdata[, c("ID", "gender", "race", "stage", "hpv16")], 
    clinical_color = cli_colors, na.value = "#F3F5F7")
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

Run shinyapp

``` r
run_virusPlot()
```

<img src="man/figures/shinyapp.png" width="100%" />

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 10 x64 (build 19045)
#> 
#> Matrix products: default
#>   LAPACK version 3.12.1
#> 
#> locale:
#> [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
#> [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
#> [3] LC_MONETARY=Chinese (Simplified)_China.utf8
#> [4] LC_NUMERIC=C                               
#> [5] LC_TIME=Chinese (Simplified)_China.utf8    
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] org.Hs.eg.db_3.21.0  AnnotationDbi_1.70.0 IRanges_2.42.0      
#> [4] S4Vectors_0.46.0     Biobase_2.68.0       BiocGenerics_0.54.0 
#> [7] generics_0.1.4       virusPlot_0.1.7     
#> 
#> loaded via a namespace (and not attached):
#>   [1] RColorBrewer_1.1-3                      
#>   [2] rstudioapi_0.17.1                       
#>   [3] jsonlite_2.0.0                          
#>   [4] magrittr_2.0.3                          
#>   [5] ggtangle_0.0.7                          
#>   [6] GenomicFeatures_1.60.0                  
#>   [7] farver_2.1.2                            
#>   [8] rmarkdown_2.29                          
#>   [9] fs_1.6.6                                
#>  [10] BiocIO_1.18.0                           
#>  [11] vctrs_0.6.5                             
#>  [12] config_0.3.2                            
#>  [13] memoise_2.0.1                           
#>  [14] Rsamtools_2.24.0                        
#>  [15] RCurl_1.98-1.17                         
#>  [16] ggtree_3.16.3                           
#>  [17] forcats_1.0.0                           
#>  [18] htmltools_0.5.8.1                       
#>  [19] S4Arrays_1.8.1                          
#>  [20] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 
#>  [21] plotrix_3.8-4                           
#>  [22] curl_7.0.0                              
#>  [23] SparseArray_1.8.1                       
#>  [24] gridGraphics_0.5-1                      
#>  [25] KernSmooth_2.23-26                      
#>  [26] aplotExtra_0.0.4                        
#>  [27] htmlwidgets_1.6.4                       
#>  [28] plyr_1.8.9                              
#>  [29] plotly_4.11.0                           
#>  [30] cachem_1.1.0                            
#>  [31] GenomicAlignments_1.44.0                
#>  [32] igraph_2.1.4                            
#>  [33] mime_0.13                               
#>  [34] lifecycle_1.0.4                         
#>  [35] pkgconfig_2.0.3                         
#>  [36] Matrix_1.7-3                            
#>  [37] R6_2.6.1                                
#>  [38] fastmap_1.2.0                           
#>  [39] GenomeInfoDbData_1.2.14                 
#>  [40] MatrixGenerics_1.20.0                   
#>  [41] shiny_1.11.1                            
#>  [42] digest_0.6.37                           
#>  [43] aplot_0.2.8                             
#>  [44] enrichplot_1.28.4                       
#>  [45] maftools_2.24.0                         
#>  [46] ps_1.9.1                                
#>  [47] patchwork_1.3.1                         
#>  [48] GenomicRanges_1.60.0                    
#>  [49] RSQLite_2.4.2                           
#>  [50] labeling_0.4.3                          
#>  [51] httr_1.4.7                              
#>  [52] abind_1.4-8                             
#>  [53] compiler_4.5.1                          
#>  [54] ggbreak_0.1.5                           
#>  [55] withr_3.0.2                             
#>  [56] attempt_0.3.1                           
#>  [57] bit64_4.6.0-1                           
#>  [58] S7_0.2.0                                
#>  [59] BiocParallel_1.42.1                     
#>  [60] DBI_1.2.3                               
#>  [61] gplots_3.2.0                            
#>  [62] R.utils_2.13.0                          
#>  [63] ChIPseeker_1.44.0                       
#>  [64] rappdirs_0.3.3                          
#>  [65] DelayedArray_0.34.1                     
#>  [66] rjson_0.2.23                            
#>  [67] DNAcopy_1.82.0                          
#>  [68] ggsci_3.2.0                             
#>  [69] gtools_3.9.5                            
#>  [70] caTools_1.18.3                          
#>  [71] tools_4.5.1                             
#>  [72] ape_5.8-1                               
#>  [73] rentrez_1.2.4                           
#>  [74] httpuv_1.6.16                           
#>  [75] R.oo_1.27.1                             
#>  [76] glue_1.8.0                              
#>  [77] callr_3.7.6                             
#>  [78] restfulr_0.0.16                         
#>  [79] nlme_3.1-168                            
#>  [80] GOSemSim_2.34.0                         
#>  [81] promises_1.3.3                          
#>  [82] shadowtext_0.1.5                        
#>  [83] grid_4.5.1                              
#>  [84] reshape2_1.4.4                          
#>  [85] fgsea_1.34.2                            
#>  [86] gtable_0.3.6                            
#>  [87] R.methodsS3_1.8.2                       
#>  [88] tidyr_1.3.1                             
#>  [89] data.table_1.17.8                       
#>  [90] XVector_0.48.0                          
#>  [91] ggrepel_0.9.6                           
#>  [92] pillar_1.11.0                           
#>  [93] stringr_1.5.1                           
#>  [94] yulab.utils_0.2.1                       
#>  [95] later_1.4.2                             
#>  [96] splines_4.5.1                           
#>  [97] dplyr_1.1.4                             
#>  [98] treeio_1.32.0                           
#>  [99] lattice_0.22-7                          
#> [100] survival_3.8-3                          
#> [101] rtracklayer_1.68.0                      
#> [102] bit_4.6.0                               
#> [103] tidyselect_1.2.1                        
#> [104] GO.db_3.21.0                            
#> [105] Biostrings_2.76.0                       
#> [106] knitr_1.50                              
#> [107] gridExtra_2.3                           
#> [108] SummarizedExperiment_1.38.1             
#> [109] xfun_0.52                               
#> [110] shinydashboard_0.7.3                    
#> [111] matrixStats_1.5.0                       
#> [112] stringi_1.8.7                           
#> [113] UCSC.utils_1.4.0                        
#> [114] lazyeval_0.2.2                          
#> [115] ggfun_0.2.0                             
#> [116] yaml_2.3.10                             
#> [117] boot_1.3-31                             
#> [118] TxDb.Hsapiens.UCSC.hg38.knownGene_3.21.0
#> [119] evaluate_1.0.4                          
#> [120] codetools_0.2-20                        
#> [121] tibble_3.3.0                            
#> [122] qvalue_2.40.0                           
#> [123] ggplotify_0.1.2                         
#> [124] cli_3.6.5                               
#> [125] xtable_1.8-4                            
#> [126] processx_3.8.6                          
#> [127] golem_0.5.1                             
#> [128] Rcpp_1.1.0                              
#> [129] GenomeInfoDb_1.44.1                     
#> [130] png_0.1-8                               
#> [131] XML_3.99-0.18                           
#> [132] parallel_4.5.1                          
#> [133] ggplot2_3.5.2                           
#> [134] blob_1.2.4                              
#> [135] DOSE_4.2.0                              
#> [136] bitops_1.0-9                            
#> [137] ggstar_1.0.4                            
#> [138] viridisLite_0.4.2                       
#> [139] tidytree_0.4.6                          
#> [140] scales_1.4.0                            
#> [141] purrr_1.1.0                             
#> [142] crayon_1.5.3                            
#> [143] rlang_1.1.6                             
#> [144] cowplot_1.2.0                           
#> [145] fastmatch_1.1-6                         
#> [146] KEGGREST_1.48.1
```

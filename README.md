
<!-- README.md is generated from README.Rmd. Please edit that file -->

# virusPlot

<!-- badges: start -->
<!-- badges: end -->

Visualization of virus insertion information

## :writing_hand: Authors

Erqiang Hu: Department of Bioinformatics, School of Basic Medical
Sciences, Southern Medical University.

Shanye Yin: Albert Einstein College of Medicine

## :hammer: Installation

``` r
devtools::install_github("huerqiang/virusPlot")
```

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
         gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
         start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
         end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))

insert_num <- data.frame(start = seq(1, 7801, 100),
    end = seq(101, 7901, 100),
    num = sample(1:28, 79, replace = TRUE))
circle_virus(virus_info, insert_num)
#> Warning: Removed 1 rows containing missing values (`geom_col()`).
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

``` r
data(insert_info)
# virus_info <- data.frame(
#       gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
#       start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
#       end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))
# HPV16
gene_features <- get_virus_annotation(accession_number = "NC_001526.2",
    email = "13766876214@163.com")
virus_info_NC_001526 <- deal_virus_annotation(gene_features)
strudel_plot(virus_info = virus_info_NC_001526, insert_info, hot_gene = 5)
#> >> done...                    2024-01-31 19:06:13
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

``` r
strudel_plot(virus_info = virus_info_NC_001526, insert_info, 
             hot_gene = c( "SAV1", "ZPLD1", "CHMP6", "IRS4"))
#> >> done...                    2024-01-31 19:06:42
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

``` r
data(vcf_matrix)
data(col)
data(pdata)
data(cli_colors)
oncoplot(vcf_matrix, varis_color = col, 
    clinical = pdata[, c("ID", "gender", "race", "stage", "hpv16")], 
    clinical_color = cli_colors, na.value = "#F3F5F7")
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

``` r
sessionInfo()
#> R version 4.3.2 (2023-10-31 ucrt)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 11 x64 (build 22621)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
#> [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
#> [3] LC_MONETARY=Chinese (Simplified)_China.utf8
#> [4] LC_NUMERIC=C                               
#> [5] LC_TIME=Chinese (Simplified)_China.utf8    
#> 
#> time zone: Asia/Shanghai
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] org.Hs.eg.db_3.18.0  AnnotationDbi_1.64.1 IRanges_2.36.0      
#> [4] S4Vectors_0.40.1     Biobase_2.62.0       BiocGenerics_0.48.1 
#> [7] virusPlot_0.1.1     
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.3.2                           
#>   [2] later_1.3.1                             
#>   [3] BiocIO_1.12.0                           
#>   [4] bitops_1.0-7                            
#>   [5] ggplotify_0.1.2                         
#>   [6] filelock_1.0.2                          
#>   [7] tibble_3.2.1                            
#>   [8] polyclip_1.10-6                         
#>   [9] XML_3.99-0.14                           
#>  [10] lifecycle_1.0.4                         
#>  [11] lattice_0.21-9                          
#>  [12] MASS_7.3-60                             
#>  [13] magrittr_2.0.3                          
#>  [14] rmarkdown_2.25                          
#>  [15] yaml_2.3.7                              
#>  [16] plotrix_3.8-4                           
#>  [17] httpuv_1.6.12                           
#>  [18] cowplot_1.1.1                           
#>  [19] DBI_1.1.3                               
#>  [20] RColorBrewer_1.1-3                      
#>  [21] abind_1.4-5                             
#>  [22] zlibbioc_1.48.0                         
#>  [23] GenomicRanges_1.54.0                    
#>  [24] purrr_1.0.2                             
#>  [25] ggraph_2.1.0                            
#>  [26] RCurl_1.98-1.12                         
#>  [27] yulab.utils_0.0.9                       
#>  [28] tweenr_2.0.2                            
#>  [29] rappdirs_0.3.3                          
#>  [30] GenomeInfoDbData_1.2.11                 
#>  [31] enrichplot_1.23.1.991                   
#>  [32] ggrepel_0.9.4                           
#>  [33] tidytree_0.4.5                          
#>  [34] rentrez_1.2.3                           
#>  [35] ChIPseeker_1.39.0                       
#>  [36] codetools_0.2-19                        
#>  [37] DelayedArray_0.28.0                     
#>  [38] DNAcopy_1.76.0                          
#>  [39] DOSE_3.29.1.991                         
#>  [40] xml2_1.3.5                              
#>  [41] ggforce_0.4.1                           
#>  [42] tidyselect_1.2.0                        
#>  [43] aplot_0.2.1                             
#>  [44] farver_2.1.1                            
#>  [45] viridis_0.6.4                           
#>  [46] matrixStats_1.0.0                       
#>  [47] BiocFileCache_2.10.1                    
#>  [48] GenomicAlignments_1.38.0                
#>  [49] jsonlite_1.8.8                          
#>  [50] ellipsis_0.3.2                          
#>  [51] tidygraph_1.2.3                         
#>  [52] survival_3.5-7                          
#>  [53] tools_4.3.2                             
#>  [54] progress_1.2.2                          
#>  [55] treeio_1.26.0                           
#>  [56] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 
#>  [57] ggstar_1.0.4                            
#>  [58] Rcpp_1.0.12                             
#>  [59] glue_1.7.0                              
#>  [60] gridExtra_2.3                           
#>  [61] SparseArray_1.2.0                       
#>  [62] xfun_0.40                               
#>  [63] qvalue_2.34.0                           
#>  [64] MatrixGenerics_1.14.0                   
#>  [65] GenomeInfoDb_1.38.1                     
#>  [66] dplyr_1.1.4                             
#>  [67] withr_3.0.0                             
#>  [68] BiocManager_1.30.22                     
#>  [69] fastmap_1.1.1                           
#>  [70] boot_1.3-28.1                           
#>  [71] fansi_1.0.6                             
#>  [72] caTools_1.18.2                          
#>  [73] digest_0.6.33                           
#>  [74] R6_2.5.1                                
#>  [75] mime_0.12                               
#>  [76] gridGraphics_0.5-1                      
#>  [77] colorspace_2.1-0                        
#>  [78] GO.db_3.18.0                            
#>  [79] gtools_3.9.4                            
#>  [80] biomaRt_2.58.0                          
#>  [81] RSQLite_2.3.1                           
#>  [82] utf8_1.2.4                              
#>  [83] tidyr_1.3.0                             
#>  [84] generics_0.1.3                          
#>  [85] ggsci_3.0.0                             
#>  [86] data.table_1.14.8                       
#>  [87] rtracklayer_1.62.0                      
#>  [88] prettyunits_1.2.0                       
#>  [89] graphlayouts_1.0.1                      
#>  [90] httr_1.4.7                              
#>  [91] S4Arrays_1.2.0                          
#>  [92] scatterpie_0.2.1                        
#>  [93] pkgconfig_2.0.3                         
#>  [94] gtable_0.3.4                            
#>  [95] blob_1.2.4                              
#>  [96] XVector_0.42.0                          
#>  [97] shadowtext_0.1.2                        
#>  [98] htmltools_0.5.7                         
#>  [99] fgsea_1.28.0                            
#> [100] scales_1.3.0                            
#> [101] TxDb.Hsapiens.UCSC.hg38.knownGene_3.18.0
#> [102] png_0.1-8                               
#> [103] ggfun_0.1.3                             
#> [104] knitr_1.45                              
#> [105] rstudioapi_0.15.0                       
#> [106] reshape2_1.4.4                          
#> [107] rjson_0.2.21                            
#> [108] nlme_3.1-163                            
#> [109] curl_5.1.0                              
#> [110] cachem_1.0.8                            
#> [111] stringr_1.5.1                           
#> [112] BiocVersion_3.18.1                      
#> [113] KernSmooth_2.23-22                      
#> [114] parallel_4.3.2                          
#> [115] HDO.db_0.99.1                           
#> [116] restfulr_0.0.15                         
#> [117] pillar_1.9.0                            
#> [118] grid_4.3.2                              
#> [119] vctrs_0.6.5                             
#> [120] gplots_3.1.3                            
#> [121] promises_1.2.1                          
#> [122] dbplyr_2.4.0                            
#> [123] xtable_1.8-4                            
#> [124] evaluate_0.23                           
#> [125] GenomicFeatures_1.54.1                  
#> [126] cli_3.6.2                               
#> [127] compiler_4.3.2                          
#> [128] Rsamtools_2.18.0                        
#> [129] rlang_1.1.3                             
#> [130] crayon_1.5.2                            
#> [131] labeling_0.4.3                          
#> [132] aplotExtra_0.0.2                        
#> [133] forcats_1.0.0                           
#> [134] maftools_2.18.0                         
#> [135] plyr_1.8.9                              
#> [136] stringi_1.7.12                          
#> [137] viridisLite_0.4.2                       
#> [138] BiocParallel_1.36.0                     
#> [139] munsell_0.5.0                           
#> [140] Biostrings_2.70.1                       
#> [141] lazyeval_0.2.2                          
#> [142] GOSemSim_2.28.0                         
#> [143] Matrix_1.6-4                            
#> [144] hms_1.1.3                               
#> [145] patchwork_1.1.3                         
#> [146] bit64_4.0.5                             
#> [147] ggplot2_3.4.4                           
#> [148] KEGGREST_1.42.0                         
#> [149] shiny_1.8.0                             
#> [150] SummarizedExperiment_1.32.0             
#> [151] interactiveDisplayBase_1.40.0           
#> [152] highr_0.10                              
#> [153] AnnotationHub_3.10.0                    
#> [154] igraph_1.5.1                            
#> [155] memoise_2.0.1                           
#> [156] ggtree_3.10.0                           
#> [157] fastmatch_1.1-4                         
#> [158] bit_4.0.5                               
#> [159] ape_5.7-1
```

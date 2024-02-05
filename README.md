
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
#> >> done...                    2024-02-05 17:45:21
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

``` r
strudel_plot(virus_info = virus_info_NC_001526, insert_info, 
             hot_gene = c( "SAV1", "ZPLD1", "CHMP6", "IRS4"))
#> >> done...                    2024-02-05 17:45:52
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

Run shinyapp

``` r
run_virusPlot()
```

<img src="man/figures/shinyapp.png" width="100%" />

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
#>  [14] plotly_4.10.3                           
#>  [15] rmarkdown_2.25                          
#>  [16] yaml_2.3.7                              
#>  [17] plotrix_3.8-4                           
#>  [18] httpuv_1.6.12                           
#>  [19] cowplot_1.1.1                           
#>  [20] DBI_1.1.3                               
#>  [21] RColorBrewer_1.1-3                      
#>  [22] golem_0.4.1                             
#>  [23] abind_1.4-5                             
#>  [24] zlibbioc_1.48.0                         
#>  [25] GenomicRanges_1.54.0                    
#>  [26] purrr_1.0.2                             
#>  [27] ggraph_2.1.0                            
#>  [28] RCurl_1.98-1.12                         
#>  [29] yulab.utils_0.0.9                       
#>  [30] tweenr_2.0.2                            
#>  [31] rappdirs_0.3.3                          
#>  [32] GenomeInfoDbData_1.2.11                 
#>  [33] enrichplot_1.23.1.991                   
#>  [34] ggrepel_0.9.4                           
#>  [35] rentrez_1.2.3                           
#>  [36] tidytree_0.4.5                          
#>  [37] ChIPseeker_1.39.0                       
#>  [38] codetools_0.2-19                        
#>  [39] DelayedArray_0.28.0                     
#>  [40] DNAcopy_1.76.0                          
#>  [41] DOSE_3.29.1.991                         
#>  [42] xml2_1.3.5                              
#>  [43] ggforce_0.4.1                           
#>  [44] tidyselect_1.2.0                        
#>  [45] aplot_0.2.1                             
#>  [46] farver_2.1.1                            
#>  [47] viridis_0.6.4                           
#>  [48] matrixStats_1.0.0                       
#>  [49] BiocFileCache_2.10.1                    
#>  [50] GenomicAlignments_1.38.0                
#>  [51] jsonlite_1.8.8                          
#>  [52] ellipsis_0.3.2                          
#>  [53] tidygraph_1.2.3                         
#>  [54] survival_3.5-7                          
#>  [55] tools_4.3.2                             
#>  [56] progress_1.2.2                          
#>  [57] treeio_1.26.0                           
#>  [58] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 
#>  [59] ggstar_1.0.4                            
#>  [60] Rcpp_1.0.12                             
#>  [61] glue_1.7.0                              
#>  [62] gridExtra_2.3                           
#>  [63] SparseArray_1.2.0                       
#>  [64] xfun_0.40                               
#>  [65] qvalue_2.34.0                           
#>  [66] MatrixGenerics_1.14.0                   
#>  [67] GenomeInfoDb_1.38.1                     
#>  [68] dplyr_1.1.4                             
#>  [69] shinydashboard_0.7.2                    
#>  [70] withr_3.0.0                             
#>  [71] BiocManager_1.30.22                     
#>  [72] fastmap_1.1.1                           
#>  [73] boot_1.3-28.1                           
#>  [74] fansi_1.0.6                             
#>  [75] caTools_1.18.2                          
#>  [76] digest_0.6.33                           
#>  [77] R6_2.5.1                                
#>  [78] mime_0.12                               
#>  [79] gridGraphics_0.5-1                      
#>  [80] colorspace_2.1-0                        
#>  [81] GO.db_3.18.0                            
#>  [82] gtools_3.9.4                            
#>  [83] biomaRt_2.58.0                          
#>  [84] RSQLite_2.3.1                           
#>  [85] config_0.3.2                            
#>  [86] ggsci_3.0.0                             
#>  [87] utf8_1.2.4                              
#>  [88] tidyr_1.3.0                             
#>  [89] generics_0.1.3                          
#>  [90] data.table_1.14.8                       
#>  [91] rtracklayer_1.62.0                      
#>  [92] htmlwidgets_1.6.3                       
#>  [93] prettyunits_1.2.0                       
#>  [94] graphlayouts_1.0.1                      
#>  [95] httr_1.4.7                              
#>  [96] S4Arrays_1.2.0                          
#>  [97] scatterpie_0.2.1                        
#>  [98] pkgconfig_2.0.3                         
#>  [99] gtable_0.3.4                            
#> [100] blob_1.2.4                              
#> [101] XVector_0.42.0                          
#> [102] shadowtext_0.1.2                        
#> [103] htmltools_0.5.7                         
#> [104] fgsea_1.28.0                            
#> [105] scales_1.3.0                            
#> [106] TxDb.Hsapiens.UCSC.hg38.knownGene_3.18.0
#> [107] png_0.1-8                               
#> [108] attempt_0.3.1                           
#> [109] ggfun_0.1.3                             
#> [110] knitr_1.45                              
#> [111] rstudioapi_0.15.0                       
#> [112] reshape2_1.4.4                          
#> [113] rjson_0.2.21                            
#> [114] nlme_3.1-163                            
#> [115] curl_5.1.0                              
#> [116] cachem_1.0.8                            
#> [117] stringr_1.5.1                           
#> [118] BiocVersion_3.18.1                      
#> [119] KernSmooth_2.23-22                      
#> [120] parallel_4.3.2                          
#> [121] HDO.db_0.99.1                           
#> [122] restfulr_0.0.15                         
#> [123] pillar_1.9.0                            
#> [124] grid_4.3.2                              
#> [125] vctrs_0.6.5                             
#> [126] gplots_3.1.3                            
#> [127] promises_1.2.1                          
#> [128] dbplyr_2.4.0                            
#> [129] xtable_1.8-4                            
#> [130] evaluate_0.23                           
#> [131] GenomicFeatures_1.54.1                  
#> [132] cli_3.6.2                               
#> [133] compiler_4.3.2                          
#> [134] Rsamtools_2.18.0                        
#> [135] rlang_1.1.3                             
#> [136] crayon_1.5.2                            
#> [137] labeling_0.4.3                          
#> [138] aplotExtra_0.0.2                        
#> [139] forcats_1.0.0                           
#> [140] maftools_2.18.0                         
#> [141] plyr_1.8.9                              
#> [142] stringi_1.7.12                          
#> [143] viridisLite_0.4.2                       
#> [144] BiocParallel_1.36.0                     
#> [145] munsell_0.5.0                           
#> [146] Biostrings_2.70.1                       
#> [147] lazyeval_0.2.2                          
#> [148] GOSemSim_2.28.0                         
#> [149] Matrix_1.6-4                            
#> [150] hms_1.1.3                               
#> [151] patchwork_1.1.3                         
#> [152] bit64_4.0.5                             
#> [153] ggplot2_3.4.4                           
#> [154] KEGGREST_1.42.0                         
#> [155] shiny_1.8.0                             
#> [156] highr_0.10                              
#> [157] SummarizedExperiment_1.32.0             
#> [158] interactiveDisplayBase_1.40.0           
#> [159] AnnotationHub_3.10.0                    
#> [160] igraph_1.5.1                            
#> [161] memoise_2.0.1                           
#> [162] ggtree_3.10.0                           
#> [163] fastmatch_1.1-4                         
#> [164] bit_4.0.5                               
#> [165] ape_5.7-1
```

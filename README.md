
<!-- README.md is generated from README.Rmd. Please edit that file -->

# virusPlot

<!-- badges: start -->
<!-- badges: end -->

Visualization of virus insertion information

## :writing_hand: Authors

Erqiang Hu

Department of Bioinformatics, School of Basic Medical Sciences, Southern
Medical University.

## :hammer: Installation

``` r
devtools::install_github("huerqiang/virusPlot")
```

## Example

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

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
data(insert_info)
virus_info <- data.frame(
      gene = c("E6", "E7", "E1", "E2", "E4", "E5", "L2", "L1", "LCR"),
      start = c(83, 562, 865, 2755, 3332, 3849, 4236, 5560, 7200),
      end = c(559, 858, 2813, 3852, 3619, 4100, 5657, 7155, 7904))
strudel_plot(virus_info, insert_info)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

``` r
data(vcf_matrix)
data(col)
data(pdata)
data(cli_colors)
oncoplot(vcf_matrix, varis_color = col, 
    clinical = pdata[, c("ID", "gender", "race", "stage", "hpv16")], 
    clinical_color = cli_colors, na.value = "#F3F5F7")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

``` r
sessionInfo()
#> R version 4.3.1 (2023-06-16 ucrt)
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
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] virusPlot_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] yulab.utils_0.1.0  tidyr_1.3.0        utf8_1.2.4         generics_0.1.3    
#>  [5] ggplotify_0.1.2    maftools_2.18.0    lattice_0.22-5     digest_0.6.33     
#>  [9] magrittr_2.0.3     shadowtext_0.1.2   evaluate_0.22      grid_4.3.1        
#> [13] RColorBrewer_1.1-3 fastmap_1.1.1      Matrix_1.6-1.1     survival_3.5-7    
#> [17] gridExtra_2.3      purrr_1.0.2        fansi_1.0.5        aplot_0.2.2       
#> [21] scales_1.2.1       ggstar_1.0.4       cli_3.6.1          rlang_1.1.1       
#> [25] munsell_0.5.0      splines_4.3.1      withr_2.5.1        cachem_1.0.8      
#> [29] yaml_2.3.7         DNAcopy_1.76.0     tools_4.3.1        memoise_2.0.1     
#> [33] dplyr_1.1.3        colorspace_2.1-0   ggplot2_3.4.4      forcats_1.0.0     
#> [37] vctrs_0.6.4        R6_2.5.1           gridGraphics_0.5-1 lifecycle_1.0.3   
#> [41] fs_1.6.3           ggfun_0.1.3        pkgconfig_2.0.3    pillar_1.9.0      
#> [45] gtable_0.3.4       glue_1.6.2         data.table_1.14.8  aplotExtra_0.0.2  
#> [49] xfun_0.40          tibble_3.2.1       tidyselect_1.2.0   rstudioapi_0.15.0 
#> [53] knitr_1.44         farver_2.1.1       htmltools_0.5.6.1  patchwork_1.1.3   
#> [57] rmarkdown_2.25     ggsci_3.0.0        labeling_0.4.3     compiler_4.3.1
```

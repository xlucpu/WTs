# WTs
Delineating the Interplay Between Oncogenic Pathways and Immunity in Anaplastic Wilms Tumors: Implications for Prognosis and Therapeutic Vulnerability

The uploaded R script was the entire analytical pipelines/codes for this  project, and were run in the evironment below:
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

Matrix products: default

Random number generation:
 RNG:     L'Ecuyer-CMRG 
 Normal:  Inversion 
 Sample:  Rejection 
 
locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] prodlim_2019.11.13  readxl_1.4.1        survminer_0.4.9     ggpubr_0.4.0        ggplot2_3.4.0       survival_3.4-0     
[7] Biobase_2.58.0      BiocGenerics_0.44.0

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3              rtracklayer_1.58.0          wateRmelon_2.4.0            pkgmaker_0.32.2            
  [5] tidyr_1.2.1                 bumphunter_1.40.0           minfi_1.44.0                bit64_4.0.5                
  [9] knitr_1.40                  irlba_2.3.5.1               DelayedArray_0.23.2         data.table_1.14.4          
 [13] KEGGREST_1.38.0             RCurl_1.98-1.9              GEOquery_2.66.0             doParallel_1.0.17          
 [17] generics_0.1.3              GenomicFeatures_1.50.2      preprocessCore_1.60.0       ScaledMatrix_1.6.0         
 [21] callr_3.7.3                 cowplot_1.1.1               usethis_2.1.6               RSQLite_2.2.18             
 [25] shadowtext_0.1.2            future_1.29.0               bit_4.0.4                   tzdb_0.3.0                 
 [29] enrichplot_1.18.0           xml2_1.3.3                  httpuv_1.6.6                SummarizedExperiment_1.28.0
 [33] assertthat_0.2.1            viridis_0.6.2               xfun_0.34                   hms_1.1.2                  
 [37] promises_1.2.0.1            fansi_1.0.3                 restfulr_0.0.15             scrime_1.3.5               
 [41] progress_1.2.2              dbplyr_2.3.0                km.ci_0.5-6                 igraph_1.3.5               
 [45] DBI_1.1.3                   geneplotter_1.76.0          htmlwidgets_1.5.4           reshape_0.8.9              
 [49] stats4_4.2.2                purrr_1.0.1                 ROC_1.74.0                  ellipsis_0.3.2             
 [53] backports_1.4.1             dplyr_1.0.10                annotate_1.76.0             gridBase_0.4-7             
 [57] biomaRt_2.54.0              sparseMatrixStats_1.10.0    MatrixGenerics_1.10.0       vctrs_0.5.0                
 [61] SingleCellExperiment_1.20.0 remotes_2.4.2               abind_1.4-5                 cachem_1.0.6               
 [65] withr_2.5.0                 ggforce_0.4.1               HDO.db_0.99.1               BSgenome_1.66.1            
 [69] treeio_1.22.0               GenomicAlignments_1.34.0    prettyunits_1.1.1           mclust_6.0.0               
 [73] cluster_2.1.4               DOSE_3.25.0.001             lazyeval_0.2.2              ape_5.6-2                  
 [77] crayon_1.5.2                genefilter_1.80.0           pkgconfig_2.0.3             tweenr_2.0.2               
 [81] GenomeInfoDb_1.34.3         nlme_3.1-160                pkgload_1.3.1               devtools_2.4.5             
 [85] globals_0.16.1              rlang_1.0.6                 lifecycle_1.0.3             miniUI_0.1.1.1             
 [89] nleqslv_3.3.3               downloader_0.4              registry_0.5-1              filelock_1.0.2             
 [93] affyio_1.68.0               BiocFileCache_2.6.0         rsvd_1.0.5                  cellranger_1.1.0           
 [97] polyclip_1.10-4             GSVA_1.46.0                 matrixStats_0.62.0          graph_1.76.0               
[101] rngtools_1.5.2              aplot_0.1.8                 base64_2.0.1                Matrix_1.5-3               
[105] KMsurv_0.1-5                carData_3.0-5               zoo_1.8-11                  Rhdf5lib_1.20.0            
[109] GlobalOptions_0.1.2         processx_3.8.0              png_0.1-7                   viridisLite_0.4.1          
[113] rjson_0.2.21                bitops_1.0-7                gson_0.0.9                  KernSmooth_2.23-20         
[117] rhdf5filters_1.10.0         Biostrings_2.66.0           blob_1.2.3                  DelayedMatrixStats_1.20.0  
[121] doRNG_1.8.2                 shape_1.4.6                 stringr_1.4.1               qvalue_2.30.0              
[125] nor1mix_1.3-0               parallelly_1.32.1           rstatix_0.7.1               gridGraphics_0.5-1         
[129] readr_2.1.3                 ggsignif_0.6.4              S4Vectors_0.36.0            beachmat_2.14.0            
[133] scales_1.2.1                memoise_2.0.1               GSEABase_1.60.0             magrittr_2.0.3             
[137] plyr_1.8.8                  zlibbioc_1.44.0             scatterpie_0.1.8            compiler_4.2.2             
[141] BiocIO_1.8.0                RColorBrewer_1.1-3          illuminaio_0.40.0           clue_0.3-62                
[145] DESeq2_1.38.0               Rsamtools_2.14.0            cli_3.4.1                   affy_1.76.0                
[149] XVector_0.38.0              urlchecker_1.0.1            listenv_0.8.0               patchwork_1.1.2            
[153] ps_1.7.2                    MASS_7.3-58.1               mgcv_1.8-41                 tidyselect_1.2.0           
[157] stringi_1.7.8               forcats_0.5.2               yaml_2.3.6                  GOSemSim_2.24.0            
[161] BiocSingular_1.14.0         askpass_1.1                 locfit_1.5-9.6              ggrepel_0.9.2              
[165] survMisc_0.5.6              grid_4.2.2                  fastmatch_1.1-3             tools_4.2.2                
[169] future.apply_1.10.0         parallel_4.2.2              circlize_0.4.15             rstudioapi_0.14            
[173] foreach_1.5.2               gridExtra_2.3               farver_2.1.1                ggraph_2.1.0               
[177] digest_0.6.30               BiocManager_1.30.19         lava_1.7.1                  shiny_1.7.3                
[181] quadprog_1.5-8              Rcpp_1.0.9                  car_3.1-1                   broom_1.0.1                
[185] GenomicRanges_1.49.0        siggenes_1.72.0             later_1.3.0                 httr_1.4.4                 
[189] AnnotationDbi_1.60.0        ComplexHeatmap_2.13.4       lumi_2.50.0                 colorspace_2.0-3           
[193] XML_3.99-0.12               fs_1.5.2                    IRanges_2.32.0              splines_4.2.2              
[197] yulab.utils_0.0.5           tidytree_0.4.1              graphlayouts_0.8.3          multtest_2.54.0            
[201] ggplotify_0.1.0             sessioninfo_1.2.2           xtable_1.8-4                jsonlite_1.8.3             
[205] ggtree_3.6.2                tidygraph_1.2.2             ggfun_0.0.8                 R6_2.5.1                   
[209] profvis_0.3.7               pillar_1.8.1                htmltools_0.5.3             mime_0.12                  
[213] NMF_0.24.0                  glue_1.6.2                  fastmap_1.1.0               clusterProfiler_4.6.0      
[217] BiocParallel_1.32.1         beanplot_1.3.1              codetools_0.2-18            fgsea_1.24.0               
[221] pkgbuild_1.3.1              utf8_1.2.2                  lattice_0.20-45             tibble_3.1.8               
[225] curl_4.3.3                  GO.db_3.16.0                openssl_2.0.4               limma_3.54.0               
[229] methylumi_2.44.0            munsell_0.5.0               GetoptLong_1.0.5            rhdf5_2.42.0               
[233] GenomeInfoDbData_1.2.9      iterators_1.0.14            HDF5Array_1.26.0            reshape2_1.4.4             
[237] gtable_0.3.1  

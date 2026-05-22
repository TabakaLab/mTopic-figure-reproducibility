# Second part of this figure was generated in python

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(dplyr)
library(purrr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)
library(this.path)

setwd(this.dir())

dir.create("../figures/Figure_E8", recursive = T)


res <- readRDS("../data/PBMCs_GO_MSigDB.RDS")
top_n_per_set <- 20

select_top_per_set_ordered <- function(mat, top_n = 20, remove_zero_terms = TRUE) {
  selected_rows <- c()
  row_source <- c()
  
  for (set_name in colnames(mat)) {
    vals <- mat[, set_name]
    if (remove_zero_terms) {
      vals <- vals[vals > 0]
    }
    vals <- sort(vals, decreasing = TRUE)
    if (length(vals) == 0) {
      next
    }
    top_terms <- names(vals)[seq_len(min(top_n, length(vals)))]
    selected_rows <- c(selected_rows, top_terms)
    row_source <- c(row_source, rep(set_name, length(top_terms)))
  }
  
  if (length(selected_rows) == 0) {
    mat_out <- matrix(0, nrow = 1, ncol = ncol(mat),
      dimnames = list("no_enriched_terms", colnames(mat))
    )
    
    return(list(mat = mat_out, row_source = factor("no_enriched_terms")))
  }
  
  mat_out <- mat[selected_rows, , drop = FALSE]
  rownames(mat_out) <- make.unique(rownames(mat_out), sep = " | dup")
  list(mat = mat_out, row_source = factor(row_source, levels = unique(row_source)))
}

go_top <- select_top_per_set_ordered(
  mat = res$GO_matrix,
  top_n = top_n_per_set,
  remove_zero_terms = TRUE
)

msigdb_top <- select_top_per_set_ordered(
  mat = res$MSigDB_matrix,
  top_n = top_n_per_set,
  remove_zero_terms = TRUE
)

go_mat_top_per_set <- go_top$mat
msigdb_mat_top_per_set <- msigdb_top$mat

max_score <- max(c(go_mat_top_per_set, msigdb_mat_top_per_set), na.rm = TRUE)

if (!is.finite(max_score) || max_score == 0) {max_score <- 1}

col_fun <- colorRamp2(
  c(0, 2,10, 15,25,40,max_score),
  c("white",brewer.pal(n = 5, name = 'Reds'),"black")
)

go_terms_to_mark <- c(
  "MF: ubiquitin-protein transferase activity [GO:0004842]",
  "MF: protein folding chaperone [GO:0044183]",
  "BP: defense response to virus [GO:0051607]",
  "BP: response to type I interferon [GO:0034340]",
  "BP: regulation of T cell tolerance induction [GO:0002664]",
  "BP: regulation of immune effector process [GO:0002697]",
  "MF: MHC class I protein binding [GO:0042288]",
  "BP: leukocyte mediated cytotoxicity [GO:0001909]",
  "BP: T cell activation [GO:0042110]",
  "CC: cytolytic granule [GO:0044194]",
  "BP: natural killer cell activation [GO:0030101] | dup2",
  "BP: cell killing [GO:0001906] | dup3",
  "BP: B cell activation [GO:0042113]",
  "MF: MHC class II receptor activity [GO:0032395] | dup1",
  "BP: regulation of cell migration [GO:0030334]",
  "BP: kynurenine metabolic process [GO:0070189]"
)

msigdb_terms_to_mark <- c(
  "GSE11057_NAIVE_VS_EFF_MEMORY_CD4_TCELL_UP | Genes up-regulated in comparison of naive T cells versus effector memory T cells.",
  "GSE11057_NAIVE_VS_CENT_MEMORY_CD4_TCELL_UP | Genes up-regulated in comparison of naive T cells versus central memory T cells.", 
  "GSE13738_RESTING_VS_BYSTANDER_ACTIVATED_CD4_TCELL_UP | Genes up-regulated in comparison of resting CD4 [GeneID=920] T cells versus bystander activated CD4 [GeneID=920] T cells. | dup1",
  "GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_UP | Genes up-regulated in comparison of naive T cells versus memory T cells.",
  "HAY_BONE_MARROW_NAIVE_T_CELL |  | dup3",
  "GAVISH_3CA_METAPROGRAM_CD4_T_CELLS_STRESS_HSP | Genes upregulated in subsets of cells of a given type within various tumors",
  "GAVISH_3CA_METAPROGRAM_CD4_T_CELLS_INTERFERON | Genes upregulated in subsets of cells of a given type within various tumors",
  "GSE11057_NAIVE_VS_CENT_MEMORY_CD4_TCELL_DN | Genes down-regulated in comparison of naive T cells versus central memory T cells.",  
  "GSE25087_TREG_VS_TCONV_ADULT_UP | Genes up-regulated in comparison of adult regulatory T cell (Treg) versus adult conventional T cells.",  
  "GSE7852_TREG_VS_TCONV_LN_UP | Genes up-regulated in comparison of lymph node regulatory T cells versus lymph node conventional T cells.",
  "GSE40068_CXCR5NEG_BCL6NEG_CD4_TCELL_VS_CXCR5POS_BCL6NEG_TFH_DN | Genes down-regulated in CXCR5- BCL6- CD4+ [GeneID=643;604;920] T cells versus CXCR5+ BCL6- [GeneID=643;604] follicular helper T cells.",
  "GSE11057_CD4_EFF_MEM_VS_PBMC_UP | Genes up-regulated in comparison of effector memory T cells versus peripheral blood mononuclear cells (PBMC). | dup1",
  "GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_DN | Genes down-regulated in comparison of naive T cells versus memory T cells. | dup4",
  "GSE26495_NAIVE_VS_PD1HIGH_CD8_TCELL_UP | Genes up-regulated in comparison of naive CD8 T cells versus PD-1 high CD8 T cells. | dup4",
  "GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_UP | Genes up-regulated in comparison of naive CD8 T cells versus PD-1 low CD8 T cells. | dup3",
  "GSE26495_NAIVE_VS_PD1LOW_CD8_TCELL_DN | Genes down-regulated in comparison of naive CD8 T cells versus PD-1 low CD8 T cells. | dup2",
  "JIANG_MELANOMA_TRM9_CD8 |",
  "HAY_BONE_MARROW_NK_CELLS |  | dup1",
  "FAN_OVARY_CL4_T_LYMPHOCYTE_NK_CELL_1 | The ovaries analyzed showed a pronounced population of CD53high/CXCR4high immune cells (Fig. 2e), including separate clusters for adaptive T lymphocytes and Natural Killer (NK) cells (CL4 and CL12), B lymphocytes (CL18), and innate immune system, such as monocytes and macrophages (CL13) | dup2",
  "AIZARANI_LIVER_C3_NK_NKT_CELLS_2 |  | dup1",
  "GSE10325_CD4_TCELL_VS_BCELL_DN | Genes down-regulated in comparison of healthy CD4 [GeneID=920] T cells versus healthy CD19 [GeneID=920] B cells.",
  "GAVISH_3CA_METAPROGRAM_B_CELLS_MHC_II | Genes upregulated in subsets of cells of a given type within various tumors",
  "GSE29618_BCELL_VS_PDC_DAY7_FLU_VACCINE_DN | Genes down-regulated in comparison of B cells from influenza vaccinee at day 7 post-vaccination versus plasmacytoid dendritic cells (pDC) at day 7 post-vaccination.",
  "TRAVAGLINI_LUNG_EREG_DENDRITIC_CELL",
  "GSE2706_UNSTIM_VS_8H_R848_DC_DN | Genes down-regulated in comparison of unstimulated dendritic cells (DC) at 0 h versus DCs stimulated with R848 for 8 h."
)

find_terms_to_mark <- function(mat, terms_to_mark, exact = FALSE, one_hit_per_term = TRUE) {
  row_ids <- rownames(mat)
  # remove make.unique suffixes like " | dup1", " | dup2"
  clean_row_ids <- sub(" \\| dup[0-9]+$", "", row_ids)
  selected_indices <- integer()
  selected_labels <- character()
  selected_query <- character()
  selected_terms <- character()
  
  for (i in seq_along(terms_to_mark)) {
    query <- terms_to_mark[i]
    if (exact) {
      hits <- which(clean_row_ids == query | row_ids == query)
    } else {
      hits <- grep(query, clean_row_ids, ignore.case = TRUE, fixed = TRUE)
      if (length(hits) == 0) {
        hits <- grep(query, row_ids, ignore.case = TRUE, fixed = TRUE)
      }
    }
    if (length(hits) == 0) {
      message("Not found: ", query)
      next
    }
    if (one_hit_per_term) {hits <- hits[1]}
    selected_indices <- c(selected_indices, hits)
    selected_labels <- c(selected_labels, rep(as.character(i), length(hits)))
    selected_query <- c(selected_query, rep(query, length(hits)))
    selected_terms <- c(selected_terms, row_ids[hits])
  }
  
  keep <- !duplicated(selected_indices)
  data.frame(
    row_index = selected_indices[keep],
    label = selected_labels[keep],
    query = selected_query[keep],
    term = selected_terms[keep],
    stringsAsFactors = FALSE
  )
}

go_marks <- find_terms_to_mark(
  mat = go_mat_top_per_set,
  terms_to_mark = go_terms_to_mark,
  exact = FALSE,
  one_hit_per_term = TRUE
)

msigdb_marks <- find_terms_to_mark(
  mat = msigdb_mat_top_per_set,
  terms_to_mark = msigdb_terms_to_mark,
  exact = FALSE,
  one_hit_per_term = TRUE
)

go_marks
msigdb_marks

go_right_anno <- rowAnnotation(
  mark = anno_mark(
    at = go_marks$row_index,
    labels = go_marks$label,
    labels_gp = gpar(fontsize = 10, fontface = "bold"),
    link_gp = gpar(lwd = 1)
  )
)

msigdb_right_anno <- rowAnnotation(
  mark = anno_mark(
    at = msigdb_marks$row_index,
    labels = msigdb_marks$label,
    labels_gp = gpar(fontsize = 10, fontface = "bold"),
    link_gp = gpar(lwd = 1)
  )
)

ht_go <- Heatmap(
  go_mat_top_per_set,
  name = "-log10(p)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_split = go_top$row_source,
  column_title = paste0("GO enrichment — top ", top_n_per_set, " terms per gene set"),
  row_title = "GO terms",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  border = FALSE,
  rect_gp = gpar(col = NA),
  right_annotation = go_right_anno
)

ht_msigdb <- Heatmap(
  msigdb_mat_top_per_set,
  name = "-log10(p)",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_split = msigdb_top$row_source,
  column_title = paste0("MSigDB enrichment — top ", top_n_per_set, " terms per gene set"),
  row_title = "MSigDB terms",
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  border = FALSE,
  rect_gp = gpar(col = NA),
  right_annotation = msigdb_right_anno
)

draw(ht_go %v% ht_msigdb, heatmap_legend_side = "right", merge_legends = TRUE)

go_number_legend <- go_marks |> dplyr::select(label, query, term)

msigdb_number_legend <- msigdb_marks |> dplyr::select(label, query, term)

go_number_legend
msigdb_number_legend

pdf("../figures/Figure_E8/Figure_E8a.pdf", width = 8, height = 8)
draw(ht_go %v% ht_msigdb, heatmap_legend_side = "right", merge_legends = TRUE)
dev.off()


# > sessionInfo()
# R version 4.5.2 (2025-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 24.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Warsaw
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-3     circlize_0.4.17        ComplexHeatmap_2.26.1  tidyr_1.3.2            purrr_1.2.1            dplyr_1.2.1           
# [7] msigdbr_26.1.0         org.Hs.eg.db_3.22.0    AnnotationDbi_1.72.0   IRanges_2.44.0         S4Vectors_0.48.0       Biobase_2.70.0        
# [13] BiocGenerics_0.56.0    generics_0.1.4         clusterProfiler_4.18.1
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3               gson_0.1.0              rlang_1.2.0             magrittr_2.0.5          clue_0.3-66            
# [6] GetoptLong_1.1.0        DOSE_4.4.0              matrixStats_1.5.0       compiler_4.5.2          RSQLite_2.4.6          
# [11] png_0.1-9               systemfonts_1.3.2       vctrs_0.7.2             reshape2_1.4.5          stringr_1.6.0          
# [16] shape_1.4.6.1           pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0           XVector_0.50.0         
# [21] enrichplot_1.30.3       bit_4.6.0               cachem_1.1.0            aplot_0.2.9             jsonlite_2.0.0         
# [26] blob_1.3.0              BiocParallel_1.44.0     cluster_2.1.8.1         parallel_4.5.2          R6_2.6.1               
# [31] stringi_1.8.7           GOSemSim_2.36.0         Rcpp_1.1.1              Seqinfo_1.0.0           assertthat_0.2.1       
# [36] iterators_1.0.14        ggtangle_0.1.1          R.utils_2.13.0          Matrix_1.7-4            splines_4.5.2          
# [41] igraph_2.2.3            tidyselect_1.2.1        qvalue_2.42.0           rstudioapi_0.18.0       doParallel_1.0.17      
# [46] codetools_0.2-20        curl_7.0.0              lattice_0.22-7          tibble_3.3.1            plyr_1.8.9             
# [51] withr_3.0.2             treeio_1.34.0           KEGGREST_1.50.0         S7_0.2.1                gridGraphics_0.5-1     
# [56] Biostrings_2.78.0       pillar_1.11.1           ggtree_4.0.4            foreach_1.5.2           ggfun_0.2.0            
# [61] ggplot2_4.0.2           scales_1.4.0            tidytree_0.4.7          glue_1.8.0              gdtools_0.5.0          
# [66] lazyeval_0.2.3          tools_4.5.2             data.table_1.18.2.1     fgsea_1.36.2            babelgene_22.9         
# [71] ggiraph_0.9.4           fs_2.0.1                Cairo_1.7-0             fastmatch_1.1-8         cowplot_1.2.0          
# [76] ape_5.8-1               colorspace_2.1-2        nlme_3.1-168            patchwork_1.3.2         cli_3.6.6              
# [81] rappdirs_0.3.4          fontBitstreamVera_0.1.1 gtable_0.3.6            R.methodsS3_1.8.2       yulab.utils_0.2.4      
# [86] digest_0.6.39           fontquiver_0.2.1        ggrepel_0.9.8           ggplotify_0.1.3         rjson_0.2.23           
# [91] htmlwidgets_1.6.4       farver_2.1.2            memoise_2.0.1           htmltools_0.5.9         R.oo_1.27.1            
# [96] lifecycle_1.0.5         httr_1.4.8              GlobalOptions_0.1.3     GO.db_3.22.0            fontLiberation_0.1.0   
# [101] bit64_4.6.0-1 
#   [1] data.table_1.18.2.1 RColorBrewer_1.1-3  pheatmap_1.0.13  
#   [1] gridExtra_2.3       ggrepel_0.9.8       ggplot2_4.0.2       data.table_1.18.2.1
#   [1] colorRamps_2.3.4    viridis_0.6.5       viridisLite_0.4.3   dplyr_1.2.1         data.table_1.18.2.1 rgexf_0.16.3        ggraph_2.2.2       
# [8] ggplot2_4.0.2       tidygraph_1.3.1     igraph_2.2.3       
#   [1] RColorBrewer_1.1-3      patchwork_1.3.2         Signac_1.16.0           Seurat_5.4.0            SeuratObject_5.3.0     
# [6] sp_2.2-1                gridExtra_2.3           viridis_0.6.5           viridisLite_0.4.3       MuDataSeurat_0.0.0.9000
# [11] SeuratData_0.2.2.9002   tidyr_1.3.2             dplyr_1.2.1             data.table_1.18.2.1     ggpubr_0.6.2           
# [16] ggplot2_4.0.2    
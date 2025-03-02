library(targets)
library(tarchetypes)
source("../../script/workflow_functions.R")

# Set target-specific options such as packages:
# tar_option_set(packages = "utils") # nolint
tar_option_set(
  packages = c(
    "tidyverse",
    "rlang",
    "ggrepel",
    "edgeR",
    "ComplexHeatmap",
    "circlize",
    "svglite"
  )
)
# End this file with a list of target objects.
list(
  ## Input
  tar_target(countFile, "./count.csv", format = "file"),
  #tar_target(countFilePath, c("../data/iLabData.csv", "../data/salmonCount.csv")),
  #tar_target(iLabCountFile, "../data/iLabData.csv", format = "file"),
  #tar_target(salmonCountFile, "../data/salmon.merged.gene_tpm.tsv", format = "file"),
  #tar_target(countFile, countFilePath, format = "file", pattern = map(countFilePath)),
  tar_target(count, read_csv(countFile)),
  tar_target(
    metaFile,
    #"../data/group_by_Sample-Speed-Biorep-Techrep.csv",
    "../../meta/group_by_Sample-Speed-Biorep-Techrep.csv",
    format = "file"
  ),
  tar_target(RINFile, "../../meta/RIN.csv", format = "file"),
  #tar_target(iLabCount, read_csv(iLabCountFile, col_names = TRUE)),
  #tar_target(
  #  salmonCount,
  #  read_tsv(salmonCountFile, col_names = TRUE) %>%
  #    select(-gene_id) %>%
  #    rename(Name = gene_name)
  #),
  tar_target(
    meta,
    read_csv(metaFile, col_names = TRUE) %>%
      mutate(Speed = str_sub(Speed, start = 2))
  ),
  tar_target(RIN, read_csv(RINFile, col_names = TRUE)),
  # tar_target(meta_RIN, meta %>% full_join(RIN, by = "Sample")),

  ## Remove low-expression genes
  tar_target(tpm_hist01, count %>% ggplot(aes(x=S01)) + geom_histogram() + scale_x_log10()),
  tar_target(tpm_hist01_file, ggsave("tpm_hist01_log.svg", plot = tpm_hist01), format = "file"),
  tar_target(removedLowExpression, removeLowExpression(count = count)),
  tar_target(
    removedLowExpressionFile,
    tar_saveTSV(
      removedLowExpression,
      filePath = "./gene_expression_table_remove_low.tsv"
    ),
    format = "file"
  ),
  tar_target(plotList, combn(meta$Sample, 2) %>% t() %>% as.data.frame()),

  ## MDS
  tar_target(correlation, cor(removedLowExpression[, -1])),
  tar_target(
    correlationFile,
    tar_saveTSV(
      as_tibble(correlation),
      filePath = "./correlation/gene_expression_correlation.tsv"
    ),
    format = "file"
  ),
  tar_target(
    scatterPlots,
    {
      for(i in seq_along(plotList[,1])) {
        plotScatter_export(removedLowExpression, plotList[i, 1], plotList[i, 2], "./correlation/scatter_log/")
      }
      return("./correlation/scatter_log/")
    }
  ),
  tar_target(MDS, runMDS(correlation = correlation)),
  tar_target(
    MDSFile,
    tar_saveTSV(MDS, filePath = "./MDS/gene_expression_MDS_pearson.tsv"),
    format = "file"
  ),
  tar_target(MDSScatter, plotMDS(MDS = MDS, meta = meta)),
  tar_target(MDSScatterSVG, ggsave("./MDS/gene_expression_MDS_pearson_scatter.svg", plot = MDSScatter), format = "file"),
  # tar_target(
  #   col_RIN,
  #   colorRamp2(
  #     c(min(meta_RIN$RIN), max(meta_RIN$RIN)),
  #     c("#ffffff", "#fe88f1")
  #   )
  # ),
  # tar_target(
  #   ha,
  #   HeatmapAnnotation(
  #     Speed = meta_RIN$Speed,
  #     Biorep = meta_RIN$Biorep,
  #     Techrep = meta_RIN$Techrep,
  #     RIN = meta_RIN$RIN,
  #     col = list(
  #       Speed = c(
  #         "Speed050" = "#afeeee",
  #         "Speed130" = "#00bfff",
  #         "Speed210" = "#0000ff",
  #         "Speed290" = "#191970"
  #       ),
  #       Biorep = c(
  #         "Biorep1" = "#ffd700",
  #         "Biorep2" = "#32cd32",
  #         "Biorep3" = "#ffb6c1"
  #       ),
  #       Techrep = c(
  #         "Techrep1" = "#ff7f50",
  #         "Techrep2" = "#48d1cc"
  #       ),
  #       RIN = col_RIN
  #     )
  #   )
  # ),
  # tar_target(
  #   ra,
  #   rowAnnotation(
  #     Speed = meta_RIN$Speed,
  #     Biorep = meta_RIN$Biorep,
  #     Techrep = meta_RIN$Techrep,
  #     RIN = meta_RIN$RIN,
  #     col = list(
  #       Speed = c(
  #         "Speed050" = "#afeeee",
  #         "Speed130" = "#00bfff",
  #         "Speed210" = "#0000ff",
  #         "Speed290" = "#191970"
  #       ),
  #       Biorep = c(
  #         "Biorep1" = "#ffd700",
  #         "Biorep2" = "#32cd32",
  #         "Biorep3" = "#ffb6c1"
  #       ),
  #       Techrep = c(
  #         "Techrep1" = "#ff7f50",
  #         "Techrep2" = "#48d1cc"
  #       ),
  #       RIN = col_RIN
  #     ),
  #     show_legend = FALSE
  #   )
  # ),
  # tar_target(
  #   cor_heatmap,
  #   Heatmap(
  #     correlation,
  #     name = "PCC",
  #     top_annotation = ha,
  #     left_annotation = ra,
  #     width = unit(15, "cm"),
  #     height = unit(15, "cm")
  #   )
  # ),
  tar_target(cor_heatmap, plotHeatMap(correlation, meta, RIN)),
  tar_target(
    cor_heatmap_file,
    {
      svg(filename = "./correlation/gene_expression_correlation_heatmap.svg", width = 10, height = 10)
      print(cor_heatmap)
      dev.off()
      return("./correlation/gene_expression_correlation_heatmap.svg")
    },
    format = "file"
  ),
  # tar_target(x, removedLowExpression %>% select(-Name) %>% as.matrix()),
  # tar_target(y, DGEList(counts = x)),
  tar_target(dgeList, makeDGEList(removedLowExpression)),
  # tar_target(design, model.matrix(~Speed+Biorep, data = meta)),
  tar_target(AveLogCPM, aveLogCPM(dgeList)),

  ## GLM

  # tar_target(y_disp, estimateDisp(dgeList, design)),
  # tar_target(fit, glmQLFit(y_disp, design)),
  # tar_target(res, glmQLFTest(fit, coef = 2:ncol(design))),
  # tar_target(
  #   df_res,
  #   topTags(res, n = Inf, adjust.method = "BH") %>%
  #     {tibble(gene_id = rownames(.), as_tibble(.$table))}
  # ),
  tar_target(result_GLM, runGLM(removedLowExpression, meta)),
  tar_target(result_GLM_file, tar_saveTSV(result_GLM, "./DEG/edger_result.tsv"), format = "file"),
  tar_target(
    FDR_hist_GLM,
    result_GLM %>%
      ggplot(aes(x=FDR)) +
      geom_histogram() +
      theme_bw()
  ),
  tar_target(FDR_hist_GLM_file, ggsave("./DEG/edger_result_fdr_hist.svg", plot = FDR_hist_GLM), format = "file"),
  tar_target(filter_fdr01, result_GLM %>% filter(FDR < 0.1) %>% tar_saveTSV("./DEG/edger_result_fdr01.tsv"), format = "file"),
  ## 2-way ANOVA
  tar_target(
    dftmp,
    as_tibble(count) %>%
      pivot_longer(-Name, names_to = "Sample", values_to = "TPM") %>%
      inner_join(meta, by = "Sample")
  ),
  tar_target(
    dfcv1,
    dftmp %>%
      group_by(Name, Speed) %>%
      summarize(
        CV_Biorep_Techrep = sd(TPM, na.rm = TRUE) / mean(TPM, na.rm = TRUE)
      )
  ),
  # tar_target(Y_forAnova, dftmp %>% filter(Name == "YAL001C") %>% pull(TPM)),
  # tar_target(
  #   Speed_forAnova,
  #   dftmp %>%
  #     filter(Name == "YAL001C") %>%
  #     pull(Speed)
  #   ),
  #   tar_target(
  #     Biorep_forAnova,
  #     dftmp %>% filter(Name == "YAL001C") %>% pull(Biorep)
  #   ),
  # tar_target(res_anova, aov(Y_forAnova ~ Speed_forAnova + Biorep_forAnova)),
  # tar_target(
  #   list_anova,
  #   dftmp %>%
  #     split(.$Name) %>%
  #     map(~ aov(TPM ~ Speed + Biorep, data = .))
  # ),
  # tar_target(
  #   df_anova,
  #   lapply(
  #     # 1 : length(list_anova),
  #     seq_along(list_anova),
  #     function(i) {
  #       tibble(
  #         Name = names(list_anova)[i],
  #         pval_Speed = summary(list_anova[[i]])[[1]]$"Pr(>F)"[1],
  #         pval_Biorep = summary(list_anova[[i]])[[1]]$"Pr(>F)"[2]
  #       )
  #     }
  #   ) %>%
  #     bind_rows() %>%
  #     mutate(
  #       FDR_Speed = p.adjust(pval_Speed, method = "BH"),
  #       FDR_Biorep = p.adjust(pval_Biorep, method = "BH")
  #     )
  # ),
  tar_target(result_2wayANOVA, run2wayANOVA(removedLowExpression, meta)),
  tar_target(result_2wayANOVA_file, tar_saveTSV(result_2wayANOVA, "./ANOVA/2way_anova_result.tsv"), format = "file"),
  tar_target(
    FDR_hist_2wayANOVA,
    tar_read(result_2wayANOVA) %>%
      pivot_longer(
        cols = c(FDR_Speed, FDR_Biorep),
        names_to = "FDR_type",
        values_to = "FDR"
      ) %>%
      ggplot(aes(x=FDR)) +
      geom_histogram() +
      theme_bw() +
      theme(aspect.ratio = 1) +
      facet_wrap(~FDR_type)
  ),
  tar_target(fdr_hist_anova, ggsave("./ANOVA/2way_anova_result_fdr_hist.svg", plot = FDR_hist_2wayANOVA), format = "file"),
  tar_target(anova_fdrSpeed01, result_2wayANOVA %>% filter(FDR_Speed < 0.1) %>% tar_saveTSV("./ANOVA/2way_anova_result_fdrSpeed01.tsv"), format = "file"),
  tar_target(anova_fdrBiorep01, result_2wayANOVA %>% filter(FDR_Biorep < 0.1) %>% tar_saveTSV("./ANOVA/2way_anova_result_fdrBiorep01.tsv"), format = "file"),
  tar_quarto(report, path = "../../report_template/report.qmd")
)

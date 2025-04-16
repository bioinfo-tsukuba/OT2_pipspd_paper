library(targets)
library(tarchetypes)
source("../../script/workflow_functions.R")

# Set target-specific options such as packages:
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
  tar_target(count, read_csv(countFile)),
  tar_target(
    metaFile,
    "../../meta/group_by_Sample-Speed-Biorep-Techrep.csv",
    format = "file"
  ),
  tar_target(RINFile, "../../meta/RIN.csv", format = "file"),
  tar_target(
    meta,
    read_csv(metaFile, col_names = TRUE) %>%
      mutate(Speed = str_sub(Speed, start = 2))
  ),
  tar_target(RIN, read_csv(RINFile, col_names = TRUE)),

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
  tar_target(dgeList, makeDGEList(removedLowExpression)),
  tar_target(AveLogCPM, aveLogCPM(dgeList)),

  ## GLM
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

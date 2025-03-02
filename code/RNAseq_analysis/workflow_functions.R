
# 関数定義
tar_saveTSV <- function(x, filePath) {
  write_tsv(x, file = filePath)
  return(filePath)
}
removeLowExpression <- function(count) {
  isLowExpression <- apply(count[, -1], 1, mean) <= 1
  return(count[isLowExpression == FALSE, ])
}
runMDS <- function(correlation) {
  MDS <- cmdscale(1 - correlation, eig = TRUE, k = 2)
  return(
    tibble(
      Sample = rownames(MDS$points),
      xPos = MDS$points[, 1],
      yPos = MDS$points[, 2]
    )
  )
}
plotMDS <- function(MDS, meta) {
  MDSPlotData <- full_join(meta, MDS, by = c("Sample"))
  MDSScatter <- ggplot(
    MDSPlotData,
    aes(x=xPos, y=yPos, color = Speed, shape = Biorep, size = Techrep)
  ) +
    geom_point(size=4) +
    geom_text_repel(aes(label = Sample)) +
    scale_color_manual(
      values = c(
        "#afeeee",
        "#00bfff",
        "#0000ff",
        "#191970"
      )
    ) +
    theme_bw() +
    theme(aspect.ratio = 1)
  return(MDSScatter)
}

plotHeatMap <- function(correlation, meta, RIN) {
    meta_RIN <- meta %>% full_join(RIN, by = "Sample")
    col_RIN <- colorRamp2(
        c(min(meta_RIN$RIN), max(meta_RIN$RIN)),
        c("#ffffff", "#fe88f1")
    )
    ha <- HeatmapAnnotation(
        Speed = meta_RIN$Speed,
        Biorep = meta_RIN$Biorep,
        Techrep = meta_RIN$Techrep,
        RIN = meta_RIN$RIN,
        col = list(
            Speed = c(
            "Speed050" = "#afeeee",
            "Speed130" = "#00bfff",
            "Speed210" = "#0000ff",
            "Speed290" = "#191970"
            ),
            Biorep = c(
            "Biorep1" = "#ffd700",
            "Biorep2" = "#32cd32",
            "Biorep3" = "#ffb6c1"
            ),
            Techrep = c(
            "Techrep1" = "#ff7f50",
            "Techrep2" = "#48d1cc"
            ),
            RIN = col_RIN
        )
    )
    ra <- rowAnnotation(
        Speed = meta_RIN$Speed,
        Biorep = meta_RIN$Biorep,
        Techrep = meta_RIN$Techrep,
        RIN = meta_RIN$RIN,
        col = list(
            Speed = c(
            "Speed050" = "#afeeee",
            "Speed130" = "#00bfff",
            "Speed210" = "#0000ff",
            "Speed290" = "#191970"
            ),
            Biorep = c(
            "Biorep1" = "#ffd700",
            "Biorep2" = "#32cd32",
            "Biorep3" = "#ffb6c1"
            ),
            Techrep = c(
            "Techrep1" = "#ff7f50",
            "Techrep2" = "#48d1cc"
            ),
            RIN = col_RIN
        ),
        show_legend = FALSE
    )
    return(
        Heatmap(
            correlation,
            name = "PCC",
            top_annotation = ha,
            left_annotation = ra,
            width = unit(15, "cm"),
            height = unit(15, "cm")
        )
    )
}

makeDGEList <- function(count) {
    countMatrix <- count %>% select(-Name) %>% as.matrix()
    DGEList(counts = countMatrix) %>% return()
}

runGLM <- function(count, meta) {
    dgeList <- makeDGEList(count)
    design <- model.matrix(~Speed + Biorep, data = meta)
    disp <- estimateDisp(dgeList, design)
    fit <- glmQLFit(disp, design)
    res <- glmQLFTest(fit, coef=2:ncol(design))
    df_res <- topTags(res, n = Inf, adjust.method = "BH") %>%
        {tibble(gene_id = rownames(.), as_tibble(.$table))}
    return(df_res)
}

run2wayANOVA <- function(count, meta) {
    dftmp <- as_tibble(count) %>%
        pivot_longer(-Name, names_to = "Sample", values_to = "TPM") %>%
        inner_join(meta, by = "Sample")
    list_anova <- dftmp %>%
        split(.$Name) %>%
        map(~ aov(TPM ~ Speed + Biorep, data = .))
    df_anova <- lapply(
        # 1 : length(list_anova),
        seq_along(list_anova),
        function(i) {
            tibble(
                Name = names(list_anova)[i],
                pval_Speed = summary(list_anova[[i]])[[1]]$"Pr(>F)"[1],
                pval_Biorep = summary(list_anova[[i]])[[1]]$"Pr(>F)"[2]
            )
        }
    ) %>%
        bind_rows() %>%
        mutate(
            FDR_Speed = p.adjust(pval_Speed, method = "BH"),
            FDR_Biorep = p.adjust(pval_Biorep, method = "BH")
        )
    return(df_anova)
}

plotScatter_export <- function(df, x, y, dir_name) {
  df_rm0 <- df[(df[x] != 0 & df[y] != 0),]
  correlation <- cor(df[x], df[y])
  #correlation <- cor(df_rm0[x], df_rm0[y])
  #correlation <- cor(log10(df_rm0[x]), log10(df_rm0[y]))
  p <- ggplot(df_rm0, aes(x=!!sym(x), y=!!sym(y))) + geom_point(size=1) + ggtitle(str_c(x, " vs. ", y,  ": r = ", correlation)) + scale_x_log10() + scale_y_log10() + theme(aspect.ratio = 1)
  ggsave(str_c(dir_name, x, '_', y, '.png'), p)
  #print(p)
}

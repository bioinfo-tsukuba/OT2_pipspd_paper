# 下準備

寒天培地があったDeck

    agarPlateDeck <- "6"

パッケージ・関数の準備

    library(tidyverse)
    library(ggsignif)

    # extract last element of a vector after splitting it by a pattern
    str_split_last <- function(string, pattern) {
      str_split(string, pattern) %>%
        map_chr(last)
    }

    extractRoute <- function(timeTable, startNode, resultTimeTable=tibble(), rowNumber=Inf) {
      nextEdge <- timeTable %>%
        filter(id < rowNumber) %>%
        filter(to == startNode) %>%
        tail(n=1)
      if (nrow(nextEdge) == 0) {
        return(resultTimeTable %>% arrange(elapsedTime_source))
      } else {
        resultTimeTable <- bind_rows(resultTimeTable, nextEdge)
        extractRoute(timeTable, nextEdge %>% pull(from), resultTimeTable, nextEdge %>% pull(id))
      }
    }

    formatTime <- function(agarPlateWell, timeTable) {
      extractRoute(timeTable, str_c(agarPlateDeck, "_", agarPlateWell)) %>%
        group_by(to) %>%
        summarise(time = max(elapsedTime_source)) %>%
        arrange(time) %>%
        mutate(
          to=case_when(
            startsWith(to, "7")  ~ "source",
            startsWith(to, "8") ~ "RNAseq",
            startsWith(to, agarPlateDeck) ~ "agarPlate",
            startsWith(to, "11") ~ "96well"
          )
        ) %>%
        pivot_wider(names_from = to, values_from = time) %>%
        mutate(agarPlateWell = agarPlateWell)
    }

データ読込

    firstLog <- read_csv("pipetting_record_first.csv", show_col_types = FALSE)
    secondLog <- read_csv("pipetting_record_second.csv", show_col_types = FALSE)
    wellAnnotation <- read_csv("./agar_plate_annotation.csv", show_col_types = FALSE)

ログの統合

    mergedLog <- bind_rows(firstLog, secondLog)

経過時間を追加

    startTime <- min(mergedLog$unix_time)
    mergedLog <- mergedLog %>% mutate(elapsedTime = unix_time - startTime)

ウェル番号・deck番号を抽出→必要な列を抽出

    mergedLog <- mergedLog %>%
      mutate(wellID = str_split_i(well, " of ", 1),
             deck = str_split_last(well, " on ")) %>%
      select(deck, wellID, operation, volume, elapsedTime)

タイムテーブルの作成

-   「well Aでmix」→ 「well Aから吸ってwell Aに注ぐ」に書き換え
-   列
    -   id: 管理用の連番
    -   from: どこから吸ったか
    -   to: どこへ注いだか
    -   elapsedTime\_source: 吸い始めた時間

<!-- -->

    timeTable <- mergedLog %>%
      bind_rows(
        mergedLog %>% filter(operation=="mix") %>% mutate(operation="dest")
      ) %>%
      arrange(elapsedTime) %>%
      mutate(operation=recode(operation, mix="source")) %>%
      rowid_to_column("id") %>%
      mutate(id=(id+1) %/% 2) %>%
      pivot_wider(
        id_cols = id,
        names_from = operation,
        values_from = c(deck, wellID, elapsedTime, volume)
      ) %>%
      rename(volume=volume_source) %>%
      select(-volume_dest) %>%
      mutate(
        from=str_c(deck_source, "_", wellID_source),
        to=str_c(deck_dest, "_", wellID_dest)
      ) %>%
      select(id, from, to, elapsedTime_source)

スポッティングに至るまでのタイムスタンプを抽出

    agarPlateWell <- timeTable %>% filter(str_starts(to, agarPlateDeck)) %>% pull(to)
    well_vec <- wellAnnotation %>% pull(name) %>% unique()
    timeStampToSpotting <- tibble()
    for(well in well_vec) {timeStampToSpotting <- timeStampToSpotting %>% bind_rows(formatTime(well, timeTable))}

# qfaの解析

## データ

    well_name_annotation <- tibble(
      row = rep(2:5, 6),
      col = ((1:24 - 1) %/% 4)+2,
      name = well_vec
    )

    with_qfa <- read_csv("../Gmp_fitness.csv", show_col_types = FALSE) %>%
      select(Row, Column, maxslp, maxslp_t, nr, nr_t, g, b) %>%
      rename(row=Row, col=Column) %>%
      full_join(well_name_annotation) %>%
      full_join(timeStampToSpotting, join_by(name == agarPlateWell)) %>%
      mutate(elapsed_time = agarPlate - source) %>%
      full_join(wellAnnotation, join_by(name == name)) %>%
      filter(!is.na(name))

    ## New names:
    ## Joining with `by = join_by(row, col)`
    ## • `` -> `...1`

    write_csv(with_qfa, file = "with_qfa.csv")

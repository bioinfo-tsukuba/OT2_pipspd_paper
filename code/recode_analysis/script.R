library(tidyverse)

# extract last element of a vector after splitting it by a pattern
str_split_last <- function(string, pattern) {
  str_split(string, pattern) %>%
    map_chr(last)
}

first_log <- read_csv("sample_pipetting_record_first.csv")
second_log <- read_csv("sample_pipetting_record_second.csv")

# add column that indicates elapsed time from unix_time
start_time <- min(first_log$unix_time)
first_log <- first_log %>%
  mutate(elapsed_time = round(unix_time - start_time, 2))
second_log <- second_log %>%
  mutate(elapsed_time = round(unix_time - start_time, 2))

# merge two logs
merged_log <- bind_rows(first_log, second_log)

# extract well & deck number from "well" column
merged_log <- merged_log %>%
  mutate(well_id = str_split_i(well, " of ", 1),
         deck = str_split_last(well, " on ")) %>%
  select(deck, well_id, operation, volume, elapsed_time)

# make timetable
timeTable <- merged_log %>%
  bind_rows(
    merged_log %>% filter(operation=="mix") %>% mutate(operation="dest")
  ) %>%
  arrange(elapsed_time) %>%
  mutate(operation=recode(operation, mix="source")) %>%
  rowid_to_column("id") %>%
  mutate(id=(id+1) %/% 2) %>%
  pivot_wider(
    id_cols = id,
    names_from = operation,
    values_from = c(deck, well_id, elapsed_time, volume)
  ) %>%
  rename(volume=volume_source) %>%
  select(-volume_dest)

# timeTable <- merged_log %>%
#   filter(operation != "mix") %>%
#   rowid_to_column("id") %>% mutate(id=(id+1) %/% 2) %>%
#   pivot_wider(
#     id_cols = id,
#     names_from = operation,
#     values_from = c(deck, well_id, elapsed_time, volume)
#   ) %>%
#   rename(volume=volume_source) %>%
#   select(-volume_dest)

timeTableForGraph <- timeTable %>%
  mutate(
    from=str_c(deck_source, "_", well_id_source),
    to=str_c(deck_dest, "_", well_id_dest)
  ) %>%
  select(id, from, to, elapsed_time_source)# %>%
  #tidygraph::as_tbl_graph(directed = TRUE, weighted=elapsed_time_source)

# extract last action of each operation for each well
last_operation <- merged_log %>%
  group_by(deck, well_id, operation) %>%
  summarise(last_elapsed_time = max(elapsed_time)) %>%
  ungroup()

wellAnnotation <- read_csv("./agar_plate_annotation.csv")


extractRoute <- function(timeTable, startNode, resultTimeTable=tibble(), rowNumber=Inf) {
  nextEdge <- timeTable %>%
    filter(id < rowNumber) %>%
    filter(to == startNode) %>%
    tail(n=1)
  if (nrow(nextEdge) == 0) {
    return(resultTimeTable %>% arrange(elapsed_time_source))
  } else {
    resultTimeTable <- bind_rows(resultTimeTable, nextEdge)
    extractRoute(timeTable, nextEdge %>% pull(from), resultTimeTable, nextEdge %>% pull(id))
  }
}

agarPlateWell <- timeTableForGraph %>% filter(str_starts(to, "9")) %>% pull(to)

formatTime <- function(agarPlateWell) {
  extractRoute(timeTableForGraph, str_c("9_", agarPlateWell)) %>%
    group_by(to) %>%
    summarise(time = max(elapsed_time_source)) %>%
    arrange(time) %>%
    mutate(
      to=case_when(
        startsWith(to, "7")  ~ "source",
        startsWith(to, "8") ~ "RNAseq",
        startsWith(to, "9") ~ "agarPlate",
        startsWith(to, "11") ~ "96well"
      )
    ) %>%
    pivot_wider(names_from = to, values_from = time) %>%
    mutate(agarPlateWell = agarPlateWell)
}

#well_vec <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "C3", "C4", "D1", "D2", "D3", "D4", "E1", "E2", "E3", "E4")
well_vec <- wellAnnotation %>% pull(name) %>% unique()
rettib <- tibble()
for(well in well_vec) {rettib <- rettib %>% bind_rows(formatTime(well))}

rettib %>%
  pivot_longer(cols=c("source", "RNAseq", "96well", "agarPlate")) %>%
  ggplot(aes(y=agarPlateWell, x=value, color=name)) +
  geom_point() +
  xlab("time") +
  ylab("Well name of agar plate")

well_name_annotation <- tibble(
  row = rep(1:6, 4),
  col = ((1:24 - 1) %/% 6)+1,
  name = well_vec
)


with_qfa <- read_tsv("GmpFit.tsv") %>%
  select(Row, Column, maxslp, maxslp_t, nr, nr_t, g, b) %>%
  rename(row=Row, col=Column) %>%
  full_join(well_name_annotation) %>%
  full_join(rettib, join_by(name == agarPlateWell)) %>%
  mutate(elapsed_time = agarPlate - source) %>%
  full_join(wellAnnotation, join_by(name == name))


ggplot(with_qfa, aes(x=agarPlate, y=nr, color=-1 * speed, shape=as_factor(biorep))) +
  geom_point()# +
  #xlim(0, 150)

ggplot(with_qfa, aes(x=as_factor(biorep), y=maxslp_t)) + geom_boxplot() +
  geom_signif(comparisons = list(c("1", "2")))

ggplot(with_qfa, aes(x=as_factor(biorep), y=maxslp_t)) + geom_boxplot() +
  geom_signif(comparisons = list(c("1", "2")), map_signif_level = TRUE, textsize = 6)

ggplot(with_qfa, aes(x=as_factor(speed), y=maxslp_t)) + geom_boxplot()# +
    #geom_signif(comparisons = list(c("1", "2"), c("1", "3")), map_signif_level = TRUE, textsize = 6)

ggplot(with_qfa %>% filter(biorep == 1), aes(x=as_factor(speed), y=maxslp)) + geom_boxplot()
ggplot(with_qfa %>% filter(biorep == 2), aes(x=as_factor(speed), y=maxslp)) + geom_boxplot()

for_anova <- with_qfa %>% mutate(biorep = as_factor(biorep), speed=as_factor(speed)) %>% select(biorep, speed, agarPlate, maxslp, maxslp_t, nr, nr_t, g, b, neighbour)
anova_model <- lm(b ~ biorep + speed + agarPlate + biorep * speed + speed * agarPlate + agarPlate * biorep + biorep * speed * agarPlate, data=for_anova)
anova_model_nointeract <- lm(b ~ biorep + speed + agarPlate + neighbour, data=for_anova)
anova_result <- anova(anova_model)
anova_result_nointeract <- anova(anova_model_nointeract)
ggplot(with_qfa, aes(x=speed, y=g, color=-1 * speed, shape=as_factor(biorep))) + geom_point() + theme_bw()
ggplot(with_qfa, aes(x=col, y=row, color=b, shape=as_factor(biorep), size=I(5))) + geom_point() + theme_bw()
ggplot(with_qfa, aes(x=col, y=row, fill=b)) + geom_tile() + theme_bw()

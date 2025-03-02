library(tidyverse)
data <- read_csv("../qfa_table.csv") %>%
  select(Row, Column, Expt.Time, Growth)

annotation <- read_csv("agar_plate_annotation.csv")
well_name_annotation <- tibble(
  row = rep(2:5, 6),
  col = ((1:24 - 1) %/% 4)+2,
  name = unique(annotation$name)
)
annotation <- annotation %>% full_join(well_name_annotation, by=join_by("name" == "name"))
data <- data %>%
  right_join(annotation, join_by("Row" == "row", "Column" == "col"))

data %>%
  ggplot(aes(x=Expt.Time, y=Growth, color=as.factor(speed))) +
  geom_point(size=0.1) +
  facet_grid(Row ~ Column)

data %>%
  ggplot(aes(x=Expt.Time, y=Growth, linetype = as.factor(techrep), color=as.factor(speed))) +
  geom_line() +
  scale_color_manual(values = c("#87ceeb", "#00bfff", "#0000ff", "#000080")) +
  facet_grid(. ~ biorep) +
  theme_bw()

data %>%
  ggplot(aes(x=Expt.Time, y=Growth, color=as.factor(speed))) +
  geom_line() +
  scale_color_manual(values = c("#87ceeb", "#00bfff", "#0000ff", "#000080")) +
  facet_grid(techrep ~ biorep) +
  theme_bw()


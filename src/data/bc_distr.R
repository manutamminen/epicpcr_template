
library(tidyverse)


bact_bc <- read_tsv("../../data/processed/final/16S_bc_tax.txt")
euk_bc <- read_tsv("../../data/processed/final/18S_bc_tax.txt")


bact_bc_counts <-
  bact_bc %>%
  group_by(Sample, BC) %>%
  summarise(n = sum(Count)) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Ribo", "Junk"), sep="_")


euk_bc_counts <-
  euk_bc %>%
  group_by(Sample, BC) %>%
  summarise(n = sum(Count)) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Ribo", "Junk"), sep="_")


bc_counts <-
  bind_rows(bact_bc_counts,
            euk_bc_counts) %>%
  select(-Junk)


bc_counts %>%
  ggplot(aes(x = id, y = n)) +
  geom_density(stat = "identity") +
  facet_grid(Sample ~ Ribo, scales = "free") +
  scale_x_log10() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))


bact_bc_taxa <-
  bact_bc %>%
  select(-Count) %>%
  count(Sample, BC) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Ribo", "Junk"), sep="_")


euk_bc_taxa <-
  euk_bc %>%
  select(-Count) %>%
  count(Sample, BC) %>%
  arrange(Sample, desc(n)) %>%
  group_by(Sample) %>%
  mutate(id = row_number()) %>%
  separate(Sample, c("Sample", "Ribo", "Junk"), sep="_")


bc_taxa <-
  bind_rows(bact_bc_taxa,
            euk_bc_taxa) %>%
  select(-Junk)

bc_taxa %>%
  ggplot(aes(x = id, y = n)) +
  geom_density(stat = "identity") +
  facet_grid(Sample ~ Ribo, scales = "free") +
  scale_x_log10() +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5))







  bact_bc %>%
    filter(str_detect(Sample, "Mock")) %>%
    select(-Count) %>%
    mutate(Mock1 = str_detect(Taxonomy, "Mock1"),
           Mock2 = str_detect(Taxonomy, "Mock2"),
           Mock3 = str_detect(Taxonomy, "Mock3"),
           Mock4 = str_detect(Taxonomy, "Mock4"),
           Mock5 = str_detect(Taxonomy, "Mock5"),
           Mock6 = str_detect(Taxonomy, "Mock6")) %>%
    select(-Taxonomy) %>%
    unique %>%
    count(Sample, Mock1, Mock2, Mock3, Mock4, Mock5, Mock6) %>%
    arrange(desc(n))


  euk_bc %>%
    filter(str_detect(Sample, "Mock")) %>%
    select(-Count) %>%
    mutate(Mock1 = str_detect(Taxonomy, "Mock1"),
           Mock2 = str_detect(Taxonomy, "Mock2"),
           Mock3 = str_detect(Taxonomy, "Mock3"),
           Mock4 = str_detect(Taxonomy, "Mock4"),
           Mock5 = str_detect(Taxonomy, "Mock5"),
           Mock6 = str_detect(Taxonomy, "Mock6")) %>%
    select(-Taxonomy) %>%
    unique %>%
    count(Sample, Mock1, Mock2, Mock3, Mock4, Mock5, Mock6) %>%
    arrange(desc(n))

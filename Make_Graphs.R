library(tidyverse)
setwd("~/Documents/BC_Files")
data_general = read_tsv("GSE72308_RAW_Methylation_Analysis_Data.txt")
data_confirm = read_tsv("GSE84207_RAW_Methylation_Analysis_Data.txt")
data_subtypes = read_tsv("GSE141338_RAW_Methylation_Analysis_Data.txt")
data_normal = read_tsv("GSE89278_RAW_Control_Methylation_Analysis_Data.txt")
all_output = read_tsv("joined_final_results_confirmation_test_all.tsv")

data_general = data_general %>%
  mutate(Dataset = "General")
data_confirm = data_confirm %>%
  mutate(Dataset = "Confirm")
data_subtypes = data_subtypes %>%
  mutate(Dataset = "Subtypes")
data_normal = data_normal %>%
  mutate(Dataset = "Normal")

data_all = list(data_general, data_confirm, data_subtypes, data_normal) %>%
  reduce(rbind)

setwd("~/Google Drive (alyssaclaire31@gmail.com)/Winter 2021/BIO 463/out/organized_out_files/")
results_confirmation = read_tsv("test_1_results/significant_genes/significant_genes/results_confirmation_test.tsv")
results_gwas = read_tsv("test_2_results/sumarized_data/significant_genes/significant_gwas_genes_summary.tsv")
results_gwas_subtypes = read_tsv("test_2_results/raw_data/subtypes_gwas_genes.tsv")
raw_subtypes = read_tsv("test_3_results/raw_data/common_genes_all_vs_normal.tsv")

gwas_genes = read_tsv("../../genes.tsv")

dir.create("../../Plots")
setwd("../../Plots")

# Histogram of all diffs for all genes
hist_data = all_output %>%
  pivot_longer(c(general_diff, confirm_diff), names_to = "category")
(hist_plot = ggplot(hist_data, aes(x = value)) +
    geom_histogram() +
    facet_wrap(~category) +
    ylab("Number of Genes") +
    xlab("Difference from Normal") +
    theme_bw())
ggsave("hist_plot.png", hist_plot)

# Scatter plot of all diffs for all genes
all_scatter_data = hist_data %>%
  rename(Difference = value) %>%
  pivot_longer(c(p_value_test_1, p_value_test_2), names_to = "P_Value_Category")
(all_scatter_plot = ggplot(all_scatter_data, aes(x = Gene, y = Difference, color = Gene %in% results_confirmation$Gene)) +
    scale_colour_manual(values = setNames(c('red','black'), c(T, F))) +
    geom_jitter(size = .5) +
    theme(axis.text.x = element_blank()) +
    labs(color = "Statistically Significant"))
ggsave("all_scatter_plot.png", all_scatter_plot)

# Density plot of all methylation values faceted by dataset
(density_plot = ggplot(data_all, aes(x = Value)) +
    geom_density(mapping = aes(y = ..scaled..)) +
    facet_wrap(~Dataset) +
    theme_bw() +
    ylab("Scaled Proportion"))
ggsave("density_plot.png", density_plot)

# Scatterplot of pvals from test 1
pval_scatter_data = results_confirmation %>%
  pivot_longer(c(general_diff, confirm_diff), names_to = "category")
(scatter_plot_1 = ggplot(pval_scatter_data %>% filter(value > 0), aes(x = Gene, y = value)) +
    geom_point() +
    facet_wrap(~category) +
    ylab("Value") +
    theme(axis.text = element_blank()))
ggsave("scatter_plot_1.png", scatter_plot_1)
(scatter_plot_2 = ggplot(pval_scatter_data %>% filter(value < 0), aes(x = Gene, y = value)) +
    geom_point() +
    facet_wrap(~category) +
    ylab("Value") +
    theme(axis.text = element_blank()))
ggsave("scatter_plot_2.png", scatter_plot_2)

# Line graph of genes with shared changes across datasets
set.seed(0)
greatest_diffs = data_all %>%
  filter(Dataset != "Normal") %>%
  filter(Gene %in% gwas_genes$Gene) %>%
  group_by(Gene) %>%
  summarize(Median = median(Value)) %>%
  left_join(data_normal %>% filter(Gene %in% gwas_genes$Gene) %>% group_by(Gene) %>% summarize(Normal_Median = median(Value))) %>%
  mutate(Difference = Median - Normal_Median) %>%
  arrange(desc(Difference)) %>%
  head(10) %>%
  pull(Gene)
line_data = data_all %>%
  filter(Gene %in% greatest_diffs) %>%
  filter(Dataset %in% c("Normal", "General", "Confirm")) %>%
  group_by(Gene, Dataset) %>%
  summarize(Median = median(Value))
(line_plot = ggplot(line_data, aes(x = Gene, y = Median, group = Dataset, color = Dataset)) +
    geom_line() +
    theme_bw() +
    theme(axis.text = element_text(size = 8, angle = 25, hjust = 1)))
ggsave("line_plot.png", line_plot)

# Boxplot of genes with greatest differences between datasets

# Bar graph of change values faceted by dataset
determine_category = function(general_positive, confirm_positive) {
  if (general_positive & confirm_positive) {
    return("both_positive")
  } else if (general_positive | confirm_positive) {
    return("one_positive")
  } else {
    return("both_negative")
  }
}
bar_data = all_output %>%
  mutate(general_positive = general_diff > 0) %>%
  mutate(confirm_positive = confirm_diff > 0) %>%
  rowwise() %>%
  mutate(category = determine_category(general_positive, confirm_positive))
(bar_plot = ggplot(bar_data, aes(x = category)) +
    geom_bar() +
    theme_bw() +
    ylab("Number of Genes") +
    xlab("Category"))
ggsave("bar_plot.png", bar_plot)


# Box plot of gwas genes with expected differences
set.seed(3)
box_data = data_all %>%
  filter(Gene %in% sample(gwas_genes$Gene, 10))
(box_plot = ggplot(box_data, aes(x = Dataset, y = Value)) +
    geom_boxplot() +
    facet_wrap(~Gene) +
    theme_bw() +
    theme(axis.text = element_text(size = 8, angle = 25)))
ggsave("box_plot.png", box_plot)


# Box plot of all significant gwas genes for subtypes
set.seed(6)
subtypes = raw_subtypes %>%
  mutate(Patient_ID = paste(Accession, Patient_ID, sep = "_")) %>%
  select(Patient_ID, Category)
subtype_box_data = data_subtypes %>%
  filter(Gene %in% sample(data_subtypes$Gene, 6)) %>%
  left_join(subtypes) %>%
  filter(!is.na(Category)) %>%
  mutate(Category = factor(Category, levels = c("TN", "HER2", "Luminal A", "Luminal B", "Luminal-HER2", "Normal")))
(subtypes_box_plot = ggplot(subtype_box_data, aes(x = Category, y = Value, fill = Category)) +
    geom_boxplot() +
    facet_wrap(~Gene) +
    theme_bw() +
    theme(axis.text.x = element_blank()))#element_text(size = 8, angle = 25, hjust = 1)))
ggsave("subtypes_box_plot.png", subtypes_box_plot)

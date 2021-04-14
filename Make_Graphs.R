#BiocManager::install("EnhancedVolcano")

library(tidyverse)
library(EnhancedVolcano)

setwd("~/Documents/BC_Files")
data_general = read_tsv("GSE72308_RAW_Methylation_Analysis_Data.txt")
data_confirm = read_tsv("GSE84207_RAW_Methylation_Analysis_Data.txt")
data_subtypes = read_tsv("GSE141338_RAW_Methylation_Analysis_Data.txt")
data_normal = read_tsv("GSE89278_RAW_Control_Methylation_Analysis_Data.txt")
all_output = read_tsv("joined_final_results_confirmation_test_all.tsv")

data_general = data_general %>%
  mutate(Dataset = "General1")
data_confirm = data_confirm %>%
  mutate(Dataset = "General2")
data_subtypes = data_subtypes %>%
  mutate(Dataset = "Subtypes")
data_normal = data_normal %>%
  mutate(Dataset = "Normal")

data_all = list(data_general, data_confirm, data_subtypes, data_normal) %>%
  reduce(rbind)

setwd("~/Google Drive (alyssaclaire31@gmail.com)/Winter 2021/BIO 463/out/organized_out_files/")
results_confirmation = read_tsv("test_1_results/significant_genes/significant_genes/results_confirmation_test.tsv")
raw_subtypes = read_tsv("test_3_results/raw_data/common_genes_all_vs_normal.tsv")
common_genes_all = read_tsv("test_3_results/raw_data/common_genes_not_her2.tsv")

gwas_genes = read_tsv("../../genes.tsv")

subtype_normal_patients = raw_subtypes %>%
  filter(Category == "Normal") %>%
  mutate(Patient_ID = str_c(Accession, Patient_ID, sep = "_")) %>%
  pull(Patient_ID)

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

png("volcano_general.png")
EnhancedVolcano(all_output, all_output$Gene, x = "general_diff", y = "p_value_test_1", title = "General1 Comparison", 
                subtitle = element_blank(), pCutoff = 9.425e-5, FCcutoff = 0.3,
                xlab = "Methylation beta change", pointSize = 1, caption = "Total = 22,579 genes",
                legendPosition = "bottom", legendLabels = c("not significant", "beta change", "p-value", "p-value and beta change"))
dev.off()

png("volcano_confirm.png")
EnhancedVolcano(all_output, all_output$Gene, x = "confirm_diff", y = "p_value_test_2", title = "General2 Comparison", 
                subtitle = element_blank(), pCutoff = 8.012e-6, FCcutoff = 0.3,
                xlab = "Methylation beta change", pointSize = 1, caption = "Total = 22,579 genes",
                legendPosition = "bottom", legendLabels = c("not significant", "beta change", "p-value", "p-value and beta change"))
dev.off()

num_inc_sig_1 = all_output %>%
  filter(general_diff > 0.3 & p_value_test_1 < 9.425e-5) %>%
  nrow()

num_dec_sig_1 = all_output %>%
  filter(general_diff < -0.3 & p_value_test_1 < 9.425e-5) %>%
  nrow()

num_inc_sig_2 = all_output %>%
  filter(confirm_diff > 0.3 & p_value_test_2 < 8.012e-6) %>%
  nrow()

num_dec_sig_2 = all_output %>%
  filter(confirm_diff < -0.3 & p_value_test_2 < 8.012e-6) %>%
  nrow()

(all_scatter_plot = ggplot(all_scatter_data, aes(x = Gene, y = Difference, color = Gene %in% results_confirmation$Gene)) +
    scale_colour_manual(values = setNames(c('red','black'), c(T, F))) +
    geom_jitter(size = .5) +
    theme(axis.text.x = element_blank()) +
    labs(color = "Statistically Significant"))
ggsave("all_scatter_plot.png", all_scatter_plot)

# Density plot of all methylation values faceted by dataset
(density_plot = ggplot(data_all %>% filter(!(Patient_ID %in% subtype_normal_patients)), aes(x = Value)) +
    geom_density(mapping = aes(y = ..scaled..)) +
    facet_wrap(~Dataset) +
    theme_bw() +
    ylab("Scaled proportion of genes") +
    xlab("Methylation beta value"))
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
greatest_diffs = data_all %>%
  filter(Dataset != "Normal") %>%
  filter(Gene %in% gwas_genes$Gene) %>%
  group_by(Gene) %>%
  summarize(Median = median(Value)) %>%
  left_join(data_normal %>% filter(Gene %in% gwas_genes$Gene) %>% group_by(Gene) %>% summarize(Normal_Median = median(Value))) %>%
  mutate(Abs_Difference = abs(Median - Normal_Median)) %>%
  arrange(desc(Abs_Difference)) %>%
  head(5) %>%
  pull(Gene)

line_data = data_all %>%
  filter(Gene %in% greatest_diffs) %>%
  mutate(Dataset = factor(Dataset, levels = c("General1", "General2", "Subtypes", "Normal")))
(violin_plot = ggplot(line_data, aes(x = Gene, y = Value, fill = Gene)) +
    geom_violin() +
    geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  facet_wrap(~ Dataset))
ggsave("violin_plot.png", violin_plot)


(line_plot = ggplot(line_data, aes(x = Gene, y = Median, group = Dataset, color = Dataset)) +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 25, hjust = 1)))
ggsave("line_plot.png", line_plot)

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
    theme(axis.text = element_text(size = 8, angle = 25)) +
    ylab("Methylation beta value"))
ggsave("box_plot.png", box_plot)


# Box plot of all significant genes for subtypes
subtype_box_data = common_genes_all %>%
  mutate(Category = factor(Category, levels = c("TN", "HER2", "Luminal A", "Luminal B", "Luminal-HER2", "Normal")))
(subtypes_box_plot = ggplot(subtype_box_data, aes(x = Category, y = Value)) + 
    geom_boxplot(aes(fill = Category)) + 
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    facet_wrap(~Gene) +
    ylab("Methylation beta value"))
ggsave("subtypes_box_plot.png", subtypes_box_plot)

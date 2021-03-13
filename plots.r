library('tidyverse')
common_genes_all = read_tsv('out/test_3_results/raw_data/common_genes_all_vs_normal.tsv')
print(common_genes_all)
"PLOTS"

plot <- ggplot(common_genes_all, aes(x=Category, y=Value)) + 
  geom_boxplot(aes(fill=Category)) + facet_grid(. ~ Gene) + theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("Methylation Level") + 
  ggtitle("Plot of Methylation Levels of Interesting Genes") + 
  labs(fill = "Cancer Subtypes")
plot

common_genes_all = read_tsv('out/test_3_results/raw_data/common_genes_not_her2.tsv')
print(common_genes_all)
"PLOTS"

plot2 <- ggplot(common_genes_all, aes(x=Category, y=Value)) + 
  geom_boxplot(aes(fill=Category)) + facet_grid(. ~ Gene) + theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  ylab("Methylation Level") + 
  ggtitle("Plot of Methylation Levels of Interesting Genes") + 
  labs(fill = "Cancer Subtypes")
plot2

# Minimal overlap analysis for the DEG outputs.

library(readr)
library(dplyr)
library(ggplot2)
library(ggVennDiagram)

read_up_genes <- function(path, gene_col = '...1', padj_col = 'p_val_adj', log2fc_col = 'avg_log2FC', padj_cutoff = 0.05, log2fc_cutoff = 0.25) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(.data[[padj_col]] < padj_cutoff, .data[[log2fc_col]] >= log2fc_cutoff) %>%
    pull(.data[[gene_col]]) %>%
    unique()
}

file1 <- '/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_gem_x/Chimera/Long/Chimera_long_proNeu_Veh_WTvsKO_DEG.csv'
file2 <- '/data/Amblerwg_Schaughency/will/2025_Mu_Integrate/Integrate_gem_x/Chimera/Long/Chimera_long_proNeu_E_WTvsKO_DEG.csv'

up1 <- read_up_genes(file1)
up2 <- read_up_genes(file2)
shared_up <- intersect(up1, up2)

write.csv(data.frame(shared_upregulated_gene = shared_up), 'shared_upregulated_genes_Chimera_E&Veh_WTvsKO.csv', row.names = FALSE)

p <- ggVennDiagram(list(Dataset_1_Up = up1, Dataset_2_Up = up2), label_alpha = 0) +
  scale_fill_gradient(low = '#F7FBFF', high = '#08519C') +
  theme_void(base_size = 14) +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5, face = 'bold')) +
  ggtitle('Overlap of Upregulated Differentially Expressed Genes')

print(p)
ggsave('upregulated_venn.pdf', p, width = 7, height = 6, device = cairo_pdf)
ggsave('upregulated_venn.png', p, width = 7, height = 6, dpi = 600)

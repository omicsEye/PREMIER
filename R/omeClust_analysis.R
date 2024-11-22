library(plyr)
library(ggplot2)

# installation: https://github.com/omicsEye/omicsArt
library(omicsArt)


setwd("~/Downloads")
metabolites_filtered = metabolites[, head(order(sapply(metabolites, var), decreasing = TRUE), 50)]

fakepcl <- list(
  meta = metadata,
  x = as.matrix(metabolites_filtered),
  ns = dim(metabolites_filtered)[2],
  nf = dim(metabolites_filtered)[1]
)
metabolites_heat_plot <- omicsArt:::pcl.heatmap(
  fakepcl,
  sqrtspace = T,
  gamma = 4,
  meta = T,
  show_colnames = F,
  show_rownames = T,
  treeheight_row = 0,
  treeheight_col = 5
)
ggsave(
  filename = 'metabolites_heat_plot.png',
  plot = metabolites_heat_plot,
  width = 7.2,
  height = 5.5,
  units = "in",
  dpi = 300
)
ggsave(
  filename = 'metabolites_heat_plot.pdf',
  plot = metabolites_heat_plot,
  width = 7.2,
  height = 5.5,
  units = "in",
  dpi = 300
)

metaCorr_results <- omicsArt::metadataCorr(
  metadata,
  entropy_threshold = 0.5,
  p_threshold = 0.001,
  cluster = T
)


if (do_write) {
  ggsave(
    filename = 'SFig1a.png',
    plot = metaCorr_results$pval_hetamap,
    width = 7.2,
    height = 6,
    units = "in",
    dpi = 300
  )
  ggsave(
    filename = 'SFig1b.png',
    plot = metaCorr_results$stat_pval_hetamap,
    width = 7.2,
    height = 6,
    units = "in",
    dpi = 300
  )
  write.table(
    metaCorr_results$P_perm,
    "SFig1a_metadataCorrelation_Pvalue.txt",
    sep = "\t",
    eol = "\n",
    quote = F,
    col.names = NA,
    row.names = T
  )
  #
  write.table(
    metadata,
    "metadata_processed.txt",
    sep = "\t",
    eol = "\n",
    quote = F,
    col.names = NA,
    row.names = T
  )
  write.table(
    proteins,
    "proteins_processed.txt",
    sep = "\t",
    eol = "\n",
    quote = F,
    col.names = NA,
    row.names = T
  )
}

# calculate similarity between samples based on omics measurements
library(vegan)
veg_dist <- as.matrix(vegdist(proteins, method = "bray", na.rm = T))

#  write the  in you computer as a tab-delimited file
write.table(
  veg_dist,
  'proteins_disatnce.txt',
  sep = "\t",
  eol = "\n",
  na = "",
  col.names = NA,
  quote = F,
  row.names = T
)


#### omeClust ######
# run omeClust from command line
# omeClust -i distance.txt --metadata metadata.txt -o omeClust_output_meta"

# omeClustviz omeClust_output_meta/adist.txt omeClust_output_meta/clusters.txt --metadata metadata_processed.txt --shapeby INTERVENTION -o omeClust_output_meta/
# omeClustviz omeClust_output_meta/adist.txt omeClust_output_meta/clusters.txt --metadata metadata_processed.txt --shapeby Visit -o omeClust_output_meta/
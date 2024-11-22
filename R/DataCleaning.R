
library(plyr)
library(ggplot2)

# installation: https://github.com/omicsEye/omicsArt
library(omicsArt)


setwd("~/Downloads")

metabolites <- read.delim(
  "PREMIER_Nightingale_Laboratory_Lipoproteins.csv",
  sep = ',',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = 1
)
colnames(metabolites)[which(names(metabolites) == "Total_BCAA")] <- "BCAA"
colnames(metabolites)[which(names(metabolites) == "Clinical_LDL_C")] <- "LDL_C"

metadata <- read.delim(
  "PREMIER_basic.csv",
  sep = ',',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = 1
)
raw_metadata <- metadata
rownames(raw_metadata) <- raw_metadata$STUDY_ID
metadata$Visit <- 1
metadata1 <- metadata
metadata$Visit <- 2
metadata2 <- metadata
metadata  <- rbind(metadata1, metadata2)
metabolites$ID <- paste0(metabolites$STUDY_ID, "_", metabolites$VISIT)
metadata$ID <- paste0(metadata$STUDY_ID, "_", metadata$Visit)
metadata$Sample_id <- metabolites[match(metadata$ID, metabolites$ID), "Sample id"]
metadata <- metadata[!is.na(metadata$Sample_id), ]
rownames(metadata) <- metadata$Sample_id
rownames(metabolites) <- metabolites$`Sample id`
metadata = subset(metadata, select = -c(Sample_id, ID))

metabolites_temp <- metabolites
metabolites = subset(metabolites, select = -c(`Sample id`, ID, STUDY_ID, VISIT))
metadata$AGE_REL[metadata$AGE_REL == ">= 71"] <- "71 = <"
metadata <- metadata[order(row.names(metadata)), ]
metabolites <- metabolites[order(row.names(metabolites)), ]
df <- cbind(metadata, metabolites)
Ile_ApoA1 <- ggplot(df, aes(
  x = Ile,
  y = ApoA1,
  color = as.character(Visit),
  shape = SEX
)) +
  geom_point(alpha = .7) +
  labs(x = "Ile", y = "ApoA1", color = "Visit") +
  theme_omicsEye()
ggsave(
  filename = 'Ile_ApoA1.pdf',
  plot = Ile_ApoA1,
  width = 3.6,
  height = 3,
  units = "in",
  dpi = 350
)


BCAA_cols <- c("BCAA", "Ile", "Leu", "Val")

BCAAs <- metabolites[, BCAA_cols]
lipoproteins <- metabolites[, !colnames(metabolites) %in% BCAA_cols]
dim(lipoproteins)

Y_var <- read.delim(
  "Y_variables.txt",
  sep = ',',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  #row.names = 1
)
lipoproteins <-  lipoproteins[, Y_var$Y]
dim(lipoproteins)

write.table(
  t(BCAAs),
  "BCAA.tsv",
  sep = "\t",
  eol = "\n",
  quote = F,
  col.names = NA,
  row.names = T
)
write.table(
  t(lipoproteins),
  "lipoproteins.tsv",
  sep = "\t",
  eol = "\n",
  quote = F,
  col.names = NA,
  row.names = T
)




cor.test(BCAAs$Ile, lipoproteins$ApoA1, method = "spearman")



## calculate diff

metabolites_base <- metabolites_temp[metabolites_temp$VISIT == 1, ]
rownames(metabolites_base) <- metabolites_base$STUDY_ID

metabolites_followup <- metabolites_temp[metabolites_temp$VISIT == 2, ]
rownames(metabolites_followup) <- metabolites_followup$STUDY_ID

comp_subjects <- intersect(metabolites_base$STUDY_ID, metabolites_followup$STUDY_ID)
metabolites_base <- metabolites_base[comp_subjects, ]
metabolites_base = subset(metabolites_base, select = -c(`Sample id`, ID, STUDY_ID, VISIT))

metabolites_followup <- metabolites_followup[comp_subjects, ]
metabolites_followup = subset(metabolites_followup,
                              select = -c(`Sample id`, ID, STUDY_ID, VISIT))

diff_df <- metabolites_followup - metabolites_base


df <- cbind(raw_metadata[rownames(diff_df), ], diff_df)
Ile_ApoA1 <- ggplot(df, aes(
  x = Ile,
  y = ApoA1,
  color = as.character(INTERVENTION),
  shape = SEX
)) +
  geom_point(alpha = .7) +
  labs(x = "Ile", y = "ApoA1", color = "INTERVENTION") +
  theme_omicsEye()
Ile_ApoA1
ggsave(
  filename = 'Ile_ApoA1_diff.pdf',
  plot = Ile_ApoA1,
  width = 3.6,
  height = 3,
  units = "in",
  dpi = 350
)


BCAA_cols <- c("BCAA", "Ile", "Leu", "Val")

diff_BCAA <- diff_df[, BCAA_cols]
diff_lipoproteins_diff <- diff_df[, !colnames(diff_df) %in% BCAA_cols]
write.table(
  t(diff_BCAA),
  "diff_BCAA.tsv",
  sep = "\t",
  eol = "\n",
  quote = F,
  col.names = NA,
  row.names = T
)
write.table(
  t(diff_lipoproteins_diff),
  "diff_lipoproteins.tsv",
  sep = "\t",
  eol = "\n",
  quote = F,
  col.names = NA,
  row.names = T
)



# btest BCAAs vs. Lipoproteins
# btest -X BCAA.tsv -Y lipoproteins.tsv -o btest_BCAA_Lipoproteins --fdr 0.05 --diagnostics-plot
# blockplot btest_BCAA_Lipoproteins/simtable.tsv btest_BCAA_Lipoproteins/X_Y.tsv --strongest 100 --similarity Spearman --axlabels "BCAAs" "Lipoproteins" --outfile btest_BCAA_Lipoproteins/heatmap_100.pdf
# blockplot btest_BCAA_Lipoproteins/simtable.tsv btest_BCAA_Lipoproteins/X_Y.tsv --strongest 180 --similarity Spearman --axlabels "BCAAs" "Lipoproteins" --outfile btest_BCAA_Lipoproteins/heatmap_180.pdf
# b_scatter --datax BCAA.tsv --datay lipoproteins.tsv --b_test btest_BCAA_Lipoproteins/X_Y.tsv --ind 0-179 --out btest_BCAA_Lipoproteins/scatters
#
# btest -X diff_BCAA.tsv -Y diff_lipoproteins.tsv -o btest_diff_BCAAs_lipoproteins --fdr 0.05 --diagnostics-plot
# blockplot btest_diff_BCAAs_lipoproteins/simtable.tsv btest_diff_BCAAs_lipoproteins/X_Y.tsv --strongest 100 --similarity Spearman --axlabels "BCAAs" "Lipoproteins" --outfile btest_diff_BCAAs_lipoproteins/heatmap_100.pdf
# blockplot btest_diff_BCAAs_lipoproteins/simtable.tsv btest_diff_BCAAs_lipoproteins/X_Y.tsv --strongest 368 --similarity Spearman --axlabels "BCAAs" "Lipoproteins" --outfile btest_diff_BCAAs_lipoproteins/heatmap_368.pdf
# b_scatter --datax diff_BCAA.tsv --datay diff_lipoproteins.tsv --b_test btest_diff_BCAAs_lipoproteins/X_Y.tsv --ind 0-367 --out btest_diff_BCAAs_lipoproteins/scatters

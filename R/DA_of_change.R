
###### lipoproteins changes in relation to intervention #####
do_write = T
metadata_temp2 <- raw_metadata[rownames(diff_lipoproteins),]
KW_results_lipoproteins <- lapply(diff_lipoproteins, function(x) {
  test <- kruskal.test(x ~ metadata_temp2$INTERVENTION) # Perform the test
  data.frame(
    Statistic = test$statistic, 
    P_value = test$p.value, 
    DF = test$parameter
  )
})

# Combine results into a single data frame
KW_results_lipoproteins_df <- do.call(rbind, KW_results_lipoproteins)
KW_results_lipoproteins_df$Omics <- "Lipoproteins"


# Function to perform Shapiro-Wilk test for each column
normality_test <- function(df) {
  results <- sapply(df, function(x) {
    if (is.numeric(x)) {
      test_result <- shapiro.test(x)
      return(test_result$p.value)
    } else {
      return(NA)  # For non-numeric columns, return NA
    }
  })
  return(results)
}

# Apply normality test to each column in the data frame
p_values <- normality_test(diff_lipoproteins)
KW_results_lipoproteins_df$Shapiro_Wilk_Test_P_value <- p_values






###### BCAA changes in relation to intervention #####

metadata_temp <- raw_metadata[rownames(diff_BCAA),]
KW_results_BCAA <- lapply(diff_BCAA, function(x) {
  test <- kruskal.test(x ~ metadata_temp$INTERVENTION) # Perform the test
  data.frame(
    Statistic = test$statistic, 
    P_value = test$p.value, 
    DF = test$parameter
  )
})
# Combine results into a single data frame
KW_results_BCAA_df <- do.call(rbind, KW_results_BCAA)
KW_results_BCAA_df$Omics <- "BCAA"

p_values <- normality_test(diff_BCAA)
KW_results_BCAA_df$Shapiro_Wilk_Test_P_value <- p_values

KW_results <- rbind(KW_results_BCAA_df, KW_results_lipoproteins_df)

KW_Leu_Treatment <- ggplot(diff_BCAA, aes(x = metadata_temp$INTERVENTION, y = diff_BCAA$Leu)) +
  geom_boxplot(fill = "#033C5A", color = "#033C5A", alpha = 0.4, outlier.colour = NA) +
  geom_jitter(width = 0.2, color = "#033C5A", alpha=.6, stroke=.1, size=1) +
  labs(title = "Changes of Leu (followup vs. baseline)", x = "Intervention", y = "Change (followup vs. baseline)") +
  omicsArt::theme_omicsEye()+
  ggplot2::annotate(
    geom = "text",
    x = Inf,
    y = Inf,
    hjust = 1,
    vjust = 1,
    label = sprintf(
      "p-value: %s\nKruskal-Wallis chi-squared: %s\ndf: %s",
      formatC(
        KW_results_BCAA$Leu$P_value,
        format = "f",
        digits = 4
      ),
      formatC(
        KW_results_BCAA$Leu$Statistic,
        format = 'f',
        digits = 2
      ),
      formatC(
        KW_results_BCAA$Leu$DF,
        format = 'f',
        digits = 0
      )
    ) ,
    color = "black",
    size = 2,
    fontface = "italic"
  )
KW_Leu_Treatment

KW_results$Significant <- KW_results$P_value < 0.05

# Add a column for -log10(p-value)
KW_results$neg_log_pval <- -log10(KW_results$P_value)

# Sort data frame by -log10(p-value) in descending order
KW_results <- KW_results[order(-KW_results$neg_log_pval), ]

# Convert Name to a factor with levels in sorted order
KW_plot <- ggplot(KW_results, aes(x = Name, y = neg_log_pval, fill = Omics)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.7, linewidth=.1) +  # Bar plot with border
  scale_fill_manual(values = c("#033C5A", "#AA9868")) +  # Custom colors
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", size = .5) +  # Threshold line
  labs(
    #title = " Significant changes among treatments measured by Kruskal Wallis test",
    x="",
    y = "-Log10(p-value)"
  ) +
  omicsArt::theme_omicsEye() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5, family = "Arial"),  # Rotate x-axis labels and set size
    axis.text.y = element_text(size = 5, family = "Arial"),  # Y-axis tick size
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.95, 0.95),                  # Legend in top-right corner inside the plot
    legend.justification = c(1, 1),                   # Adjust alignment to top-right
    axis.title.y = element_text(size = 7, family = "Arial")   # Customize y-axis label
  )


if (do_write) {
  ggsave(
    filename = 'KW_plot.png',
    plot = KW_plot,
    width = 7.2,
    height = 2,
    units = "in",
    dpi = 350
  )
  ggsave(
    filename = 'KW_plot.pdf',
    plot = KW_plot,
    width = 7.2,
    height = 2,
    units = "in",
    dpi = 350
  )
  ggsave(
    filename = 'KW_plot.png',
    plot = KW_plot,
    width = 7.2,
    height = 2,
    units = "in",
    dpi = 350
  )
  ggsave(
    filename = 'KW_plot.pdf',
    plot = KW_plot,
    width = 7.2,
    height = 2,
    units = "in",
    dpi = 350
  )
  write.table(
    KW_results,
    "KW_results.tsv",
    sep = "\t",
    eol = "\n",
    quote = F,
    col.names = NA,
    row.names = T
  )
  ggsave(
    filename = 'KW_Leu_Treatment.png',
    plot = KW_Leu_Treatment,
    width = 2.4,
    height = 2,
    units = "in",
    dpi = 350
  )
  ggsave(
    filename = 'KW_Leu_Treatment.pdf',
    plot = KW_Leu_Treatment,
    width = 2.4,
    height = 2,
    units = "in",
    dpi = 350
  )
}

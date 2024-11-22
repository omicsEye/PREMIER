# https://github.com/himelmallick/Tweedieverse
library(Tweedieverse)

# https://github.com/biobakery/Maaslin2
library(Maaslin2)

tweedieverse_results <- Tweedieverse::Tweedieverse(
  metabolites,
  metadata,
  output = "~/Downloads/Tweedieverse_output",
  fixed_effects = c(
    "INTERVENTION",
    "Visit",
    "SEX",
    "AGE_REL",
    "RACE_REL",
    "BMI_BASE"
  ),
  random_effects = c("CENTER", "STUDY_ID"),
  reference = "INTERVENTION,A;AGE_REL,=< 30",
  max_significance = 0.05,
  plot_scatter = T,
  plot_heatmap = T
)

Maaslin2_results <- Maaslin2::Maaslin2(
  metabolites,
  metadata,
  output = "~/Downloads/Maaslin2_output",
  fixed_effects = c(
    "INTERVENTION",
    "Visit",
    "SEX",
    "AGE_REL",
    "RACE_REL",
    "BMI_BASE"
  ),
  random_effects = c("CENTER", "STUDY_ID"),
  max_significance = 0.05,
  reference = "INTERVENTION,A;AGE_REL,=< 30"
)
for (INT in unique(metadata$INTERVENTION))
  tweedieverse_results <- Tweedieverse::Tweedieverse(
    metabolites,
    metadata[metadata$INTERVENTION == INT, ],
    output = paste0("~/Downloads/Tweedieverse_outpu_", INT),
    fixed_effects = c("Visit", "SEX", "AGE_REL", "RACE_REL", "BMI_BASE"),
    random_effects = c("CENTER", "STUDY_ID"),
    max_significance = 0.05,
    reference = "AGE_REL,=< 30",
    plot_scatter = T,
    plot_heatmap = T
  )

for (INT in unique(metadata$INTERVENTION))
  Maaslin2_results <- Maaslin2::Maaslin2(
    metabolites,
    metadata[metadata$INTERVENTION == INT, ],
    output = paste0("~/Downloads/Maaslin2_output_", INT),
    fixed_effects = c("Visit", "SEX", "AGE_REL", "RACE_REL", "BMI_BASE"),
    random_effects = c("CENTER", "STUDY_ID"),
    max_significance = 0.05,
    reference = "AGE_REL,=< 30"
  )

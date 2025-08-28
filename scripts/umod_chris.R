#=========================================#
#   UMOD role in eGFR decline in CHRIS
#=========================================#

# CKDGen lead SNP from Cristian (hg38)
# rs13335818, chr16:20348509

# COJO independent SNPs at UMOD (Stanzik, 2021)
# rsid, GRCh37 (study build), GRCh38
# rs66487676,  16:21089905, 16:21078584
# rs112709995, 16:20390734, 16:20379412
# rs34882080,  16:20361441, 16:20350119
# rs77924615,  16:20392332, 16:20381010  # CHRIS index

# Locus, lead SNP, GRCh37, GRCh38 (Liu, 2025)
# chr16:20180730-21200153, rs9928936, 16:20353049, 16:20341727

# PhD/papers/susztak/suppl. data in zip
# Credible sets from fine-mapping in EUA (Liu, 2025)
# locus: chr16_20080730_21412460
# rs7498751,  16:20354280
# rs28362063,	16:20365012
# rs4997081,	16:20365234
# rs13334589,	16:20366459
# rs13329952,	16:20366507
# rs77924615,	16:20392332 ***
# rs113009812, 16:20392890
# rs12598584,	16:20488107
# rs163277,	  16:20671260
# rs9926773,	16:21088023
# rs7186298,	16:21088031


# inputs
path_umod_chris <- "/scratch/dariush.ghasemi/projects/UMOD/data/umod_hg38.dosage"
path_clinicals <- "~/projects/shroom3/data/chris_q-norm.csv"
path_pcs <- "~/projects/shroom3/data/CHRIS13K.GT.evecs"

# outputs
rds_summary <- "/scratch/dariush.ghasemi/projects/UMOD/results/27-Aug-25_models_summaries.RDS"

#----------------------#
# inverse-normal transformation
do_INT = function(values) {
  qnorm((rank(values, na.last = "keep", ties.method = "random") - 0.5) / sum(!is.na(values)))
}


#----------------------#
# read dosage of UMOD in CHRIS
headers <- c("AID", "SNP", "DS")
umod_dosage <- read.delim(path_umod_chris, col.names = headers, sep = "\t", stringsAsFactors = FALSE)

# convert dosage file to wide format
dosage <- umod_dosage %>% pivot_wider(names_from = SNP, values_from = DS)

snps <- dosage %>% dplyr::select(starts_with("chr16:")) %>% colnames()

# load chris phenotypes
chris <- read.csv(path_clinicals, stringsAsFactors = FALSE)
pcs <- read.delim(path_pcs, stringsAsFactors = FALSE)

#----------------------#
# prepare data for model/visualization
umod_combined <- chris %>%
  dplyr::select(
    AID, Sex, Age, eGFRw.log.Res #, eGFR, UACR 
    ) %>%
  dplyr::filter(!is.na(eGFRw.log.Res)) %>% # remove rows with missing phenotype
  dplyr::mutate(across(Age, as.numeric)) %>% # convert age to numeric
  inner_join(
    umod_dosage, # merge phenotype with genotype 
    join_by(AID)
    ) %>%
  dplyr::rename(dosage = DS) %>%
  left_join(pcs, join_by(AID == X.IND_ID))


#group_by(Age_cat) %>%
#summarise(n(), M = mean(DS), SD = sd(DS))

#tab <- as.matrix(prop.table(table(UMOD$DS_Level, UMOD$Age_cat), margin = 2))

#res_cor <- cor.test(umod_combined$Age, umod_combined$dosage, method = "spearman")


#----------------------#
# plots

# age vs. SNP
umod_combined %>%
  # categorize numeric variables for model/visualization
  dplyr::mutate(
    across(
      all_of(snps), 
      function(ds) cut(ds, breaks=c(-Inf, 0.500, 1.500, Inf),
                           labels=c("0", "1", "2"))),
    age_cat = cut(
      Age, 
      breaks=c(18, 30, 40, 50, 60, 70, 94),
      labels=c("18-30","30-40","40-50","50-60","60-70","70-94"))
  ) %>%
  ggplot(aes(age_cat, snps[1])) +
  geom_violin(aes(fill = age_cat), trim = FALSE, alpha=0.9) +
  labs(
    x = "Categorized age of CHRIS participants",
    y = "Dosage of UMOD lead SNP\n(rs13335818)",
    #y = "CHRIS lead SNP dosage\n(rs77924615)",
    fill = "Categorized\nAge"
       ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.ticks.length = unit(5, "points"),
    axis.title = element_text(size = 12, face = 2)
  )


ggsave(filename = "19-Aug-25_violin_plot_chris_umod_vs_age.png", height = 5.5, width = 9, dpi = 150)

#----------------------#
# linear model 

# different models
base_model <- "do_INT(eGFRw.log.Res) ~ dosage"
inclusive_model <- "do_INT(eGFRw.log.Res) ~ dosage + Age"
interaction_model <- "do_INT(eGFRw.log.Res) ~ dosage * Age"


my_lm <- function(df, formool) {
  lm(formula = formool, data = df)
}

# to include other covariates
#+ Age + Sex + PC1 + PC2 + PC3 + PC4 + PC5

# fit models
umod_model <- lm(interaction_model, data = umod_combined)


# illustrate the model summary
broom::tidy(umod_model) %>%
  dplyr::mutate(across(c(estimate:statistic), ~round(.x, 3))) %>%
  DT::datatable()


# iterate models for all SNPs
results <- umod_combined %>%
  dplyr::select(- Sex, - starts_with("PC")) %>%
  group_by(SNP) %>%
  nest() %>%
  dplyr::mutate(
    model_base = map(data, function(df) lm(base_model, data = df)),
    model_incl = map(data, function(df) lm(inclusive_model, data = df)),
    model_intr = map(data, function(df) lm(interaction_model, data = df)),
    fitnes_base = map(model_base, broom::tidy),
    fitnes_incl = map(model_incl, broom::tidy),
    fitnes_intr = map(model_intr, broom::tidy)
  )


results %>%
  select(SNP, fitnes_intr) %>%
  unnest(fitnes_intr) %>%
  write.csv()

saveRDS(results, rds_summary)

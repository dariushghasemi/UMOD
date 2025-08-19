#=========================================#
#   UMOD role in eGFR decline in CHRIS
#=========================================#


# inputs
path_umod_chris <- "/scratch/dariush.ghasemi/projects/UMOD/data/chr16_20381010_rs77924615.txt"
path_clinicals <- "~/projects/shroom3/data/chris_q-norm.csv"


# outputs


#----------------------#
# read dosage of UMOD in CHRIS
headers <- c("AID", "ID", "DS")
umod_dosage <- read.delim(path_umod_chris, sep = "\t", col.names = headers, stringsAsFactors = FALSE)

# load chris phenotypes
chris <- read.csv(path_clinicals, stringsAsFactors = FALSE)


#----------------------#
# prepare data for model/visualization
umod_combined <- chris %>%
  dplyr::select(
    AID, Sex, Age, UACR, eGFRw.log.Res, eGFR
    ) %>%
  dplyr::filter(!is.na(eGFRw.log.Res)) %>% # remove rows with missing phenotype
  dplyr::mutate(across(Age, as.numeric)) %>% # convert age to numeric
  inner_join(
    umod_dosage %>% dplyr::select(AID, DS), # merge phenotype with genotype 
    join_by(AID)
    ) %>%
  # categorize numeric variables for model/visualization
  dplyr::mutate(
    ds_cat  = cut(
      DS,
      breaks=c(-Inf, 0.500, 1.500, Inf),
      labels=c("0", "1", "2")),
    age_cat = cut(
      Age, 
      breaks=c(18, 30, 40, 50, 60, 70, 94),
      labels=c("18-30","30-40","40-50","50-60","60-70","70-94"))
  )


#group_by(Age_cat) %>%
#summarise(n(), M = mean(DS), SD = sd(DS))

#tab <- as.matrix(prop.table(table(UMOD$DS_Level, UMOD$Age_cat), margin = 2))


#----------------------#
# plots

# age vs. SNP
umod_combined %>%
  ggplot(aes(age_cat, DS)) +
  geom_violin(aes(color = age_cat), trim = FALSE, alpha=0.9) +
  labs(
    x = "Categorized age of CHRIS participants",
    y = "Dosage of UMOD lead SNP\n(rs13335818)",
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

umod_model <- lm(eGFRw.log.Res ~ DS, data = umod_combined)

# illustrate the model summary
broom::tidy(umod_model) %>%
  dplyr::mutate(across(c(estimate:statistic), ~round(.x, 3))) %>%
  DT::datatable()



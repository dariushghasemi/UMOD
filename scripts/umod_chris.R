#=========================================#
#   UMOD role in eGFR decline in CHRIS
#=========================================#


#-----------------------------------------------------#
#---------------------- UMOD -------------------------
#-----------------------------------------------------#
UMODvcf <- read.delim("D:\\Dariush\\PhD\\Analysis\\Data\\chr16-20348509-rs13335818.txt", 
                      sep = "\t", stringsAsFactors = FALSE)

UMOD <-  merge( chris[c("AID","Sex","Age","eGFRw.log.Res","eGFR","UACR")],
                UMODvcf[c("AID", "DS")], 
                by="AID", all=FALSE)

UMOD <- UMOD %>% 
  mutate(DS_Level = cut(DS, breaks=c(-Inf, 0.500, 1.500, Inf), labels=c("0", "1", "2"))) %>% 
  mutate(Age_cat = cut(Age, breaks=c(18, 30, 40, 50, 60, 70, 94),
                       labels=c("18-30","30-40","40-50","50-60","60-70","70-94"))) #%>%
#group_by(Age_cat) %>%
#summarise(n(), M = mean(DS), SD = sd(DS))
ggplot(aes(Age_cat, DS)) +
  geom_violin(aes(fill = Age_cat), trim = FALSE, alpha=0.9) + theme_bw()

tab <- as.matrix(prop.table(table(UMOD$DS_Level, UMOD$Age_cat), margin = 2))



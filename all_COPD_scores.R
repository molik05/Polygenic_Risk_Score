library(data.table)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)

#ziskani vsech skore pro dany trait do tabulky pro kazdeho pacienta
soubory <- list.files(path = paste0(getwd(),"/final_COPD"), pattern = "\\.tsv$", full.names = TRUE)
seznam_dt <- lapply(soubory, fread)
jmena_souboru <- basename(soubory)
vysledek <- rbindlist(seznam_dt, idcol = "Soubor")
vysledek[, Soubor := rep(jmena_souboru, times = sapply(seznam_dt, nrow))]

#write_xlsx(export,"COPD_results.xlsx")
samples_with_at_least_1 <- unique(vysledek$Soubor)
samples_with_at_least_1_from_selected <- unique(data_t$Soubor)
chybi <- setdiff(samples_with_at_least_1,samples_with_at_least_1_from_selected)

#nacteni dat a rozdeleni na dilci casti pro prehlednost
data_t <- vysledek %>% 
  filter(`Included SNPs` >= 10)
data_trait <- data_t[1:344,]
data_trait <- data_t[345:680,]
data_trait <- data_t[681:997,]
data_trait <- data_t[998:1319,]
data_trait <- data_t[1320:1471,]
data_trait <- data_t

#zjisteni pro ktera skore jsou a nejsou skore u kazdeho vzorku
df <- data_t[,c("Study ID","Soubor")]
df[, value := 1]
df_wide <- dcast(  df,  Soubor ~ `Study ID`,  value.var = "value",  fill = 0)
df_wide[, Total := rowSums(.SD), .SDcols = !c("Soubor")]
total_row <- data.table(
  Soubor = "TOTAL",
  GCST004147 = sum(df_wide$GCST004147),
  GCST007692 = sum(df_wide$GCST007692),
  GCST011766 = sum(df_wide$GCST011766),
  GCST90016586 = sum(df_wide$GCST90016586),
  GCST90016588 = sum(df_wide$GCST90016588),
  GCST90016589 = sum(df_wide$GCST90016589),
  GCST90016591 = sum(df_wide$GCST90016591),
  GCST90016594 = sum(df_wide$GCST90016594),
  Total = sum(df_wide$Total)
)
df_final <- rbind(df_wide, total_row, fill = TRUE)

percentile_box_plots(data_trait)
PRS_box_plots(data_trait)

zpracuj_skore <- function(s) {
  if (is.na(s)) return(NA_real_)
  s <- as.character(s)
  if (grepl("-", s)) {
    cisla <- as.numeric(unlist(strsplit(gsub(" ", "", s), "-")))
    return(mean(cisla))
  } else {
    return(as.numeric(s))
  }
}

percentile_box_plots <- function(data_trait){
  data_trait[, Skore := sapply(Percentile, zpracuj_skore)]
  mediany <- data_trait[, .(MedSkore = median(Skore, na.rm = TRUE)), by = Sample]
  mediany[, Barva := fifelse(MedSkore <= 20, "0-20",
                             fifelse(MedSkore <= 50, "21-50", fifelse(MedSkore <= 80, "51-80", "81-100")))]
  data_trait <- merge(data_trait, mediany[, .(Sample, Barva)], by = "Sample", all.x = TRUE)
  ggplot(data_trait, aes(x = factor(Sample), y = Skore, fill = Barva)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.7, color = "black") +
    scale_fill_manual(values = c("0-20" = "lightblue", "21-50"= "forestgreen", "51-80" = "gold", "81-100" = "firebrick")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 58, hjust = 1),
          legend.position = "right") +
    labs(x = "Sample",
         y = "Percentile",
         fill = "",
         title = "COPD Risk Percentile 241... Samples"
    )
}

PRS_box_plots <- function(data_trait){
  mediany <- data_trait[, .(MedSkore = median(`Polygenic Risk Score`, na.rm = TRUE)), by = Sample]
  mediany[, Barva := fifelse(MedSkore <= 1, "lower",
                             fifelse(MedSkore <= 1.05, "slightly higher", "higher"))]
  data_trait <- merge(data_trait, mediany[, .(Sample, Barva)], by = "Sample", all.x = TRUE)
  set.seed(123) # pro reprodukovatelnost
  data_trait[, `Study ID` := factor(`Study ID`, 
                                    levels = sample(levels(factor(`Study ID`))))]
  ggplot(data_trait, aes(x = factor(Sample), y = `Polygenic Risk Score`, fill=Barva)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = `Study ID`),   # Barva bodů podle Study ID
                width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = c("lower" = "lightgreen", "slightly higher"= "lightyellow", "higher" = "#F08080")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 58, hjust = 1),
          legend.position = "right") +
    labs(x = "Sample",
         y = "Score",
         color = "PRS Study ID",
         title = "COPD Risk Score 241... Samples"
    )
}

#pro percentile bez semaforu
# data_trait[, Skore := sapply(Percentile, zpracuj_skore)]
# mediany <- data_trait[, .(MedSkore = median(Skore, na.rm = TRUE)), by = Sample]
# mediany[, Barva := fifelse(MedSkore <= 1, "lower",
#                            fifelse(MedSkore <= 1.05, "slightly higher", "higher"))]
# data_trait <- merge(data_trait, mediany[, .(Sample, Barva)], by = "Sample", all.x = TRUE)
# set.seed(123) # pro reprodukovatelnost
# data_trait[, `Study ID` := factor(`Study ID`, 
#                                   levels = sample(levels(factor(`Study ID`))))]
# ggplot(data_trait, aes(x = factor(Sample), y = Skore)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(aes(color = `Study ID`),   # Barva bodů podle Study ID
#               width = 0.2, alpha = 0.7) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 58, hjust = 1),
#         legend.position = "right") +
#   labs(x = "Sample",
#        y = "Score",
#        color = "PRS Study ID",
#        title = "COPD Risk Score 1-60 Samples"
#   )
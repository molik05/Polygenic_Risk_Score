library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)

#nacteni vsech dat ze slozky
#pacienti plicni kohorty
soubory <- list.files(path = paste0(getwd(),"/final_COPD"), pattern = "\\.tsv$", full.names = TRUE)
type = "l"
#pacienti zdravi
soubory <- list.files(path = paste0(getwd(),"/healthy_results"), pattern = "\\.tsv$", full.names = TRUE)
type = "h"
seznam_dt <- lapply(soubory, fread)
vysledek <- rbindlist(seznam_dt, idcol = "soubor")
data_trait <- vysledek %>% 
  filter(`Included SNPs` >= 10)

#spusteni funkci a tvorba grafu #ls,lp,hs,hp
hs <- PRS_distribution_plot(data_trait,type)
hp <- Percentile_distribution_plot(data_trait,type)
ls <- PRS_distribution_plot(data_trait,type)
lp <- Percentile_distribution_plot(data_trait,type)
grid.arrange(ls,lp,hs,hp,ncol=2)

sample_counts <- table(data_trait$Sample)
samples_with_3_or_more <- names(sample_counts[sample_counts >= 3])
filtered_data <- data_trait[data_trait$Sample %in% samples_with_3_or_more, ]

filtered_data[, `Polygenic Risk Score` := as.numeric(`Polygenic Risk Score`)]
mediany_PRS <- filtered_data[, .(MedSkore = mean(`Polygenic Risk Score`, na.rm = TRUE)), by = Sample]

filtered_data[, Skore := sapply(Percentile, zpracuj_skore)]
mediany_Per <- filtered_data[, .(MedSkore = mean(Skore, na.rm = TRUE)), by = Sample]

mediany_all <- merge(mediany_PRS, mediany_Per, by = "Sample", suffixes = c("_PRS", "_Per"))

ggplot(mediany_all, aes(y = MedSkore_PRS, x = MedSkore_Per)) +
  geom_point(size = 3) +
  theme_minimal()


#zpracovani percentilu v pripade rozsahu
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

#vytvoreni rolzozeni PRS skore, kde jsou vyrazeni pacienti s mene jak 3
PRS_distribution_plot <- function(data_trait,type){
  sample_counts <- table(data_trait$Sample)
  samples_with_3_or_more <- names(sample_counts[sample_counts >= 3])
  filtered_data <- data_trait[data_trait$Sample %in% samples_with_3_or_more, ]
  filtered_data[, `Polygenic Risk Score` := as.numeric(`Polygenic Risk Score`)]
  mediany <- filtered_data[, .(MedSkore = mean(`Polygenic Risk Score`, na.rm = TRUE)), by = Sample]
  mu <- mean(mediany$MedSkore)
  sigma <- sd(mediany$MedSkore)
  
  shapiro_result <- shapiro.test(mediany$MedSkore)
  p_val <- shapiro_result$p.value
  normality_text <- if (p_val < 0.05) {
    paste0("Shapiro-Wilk p=", signif(p_val, 3), " -> NOT normal distribution")
  } else {
    paste0("Shapiro-Wilk p=", signif(p_val, 3), " -> normal distribution")
  }
  
  if (type == "h"){
    title <- "Distribution of PRS of normal cohort patients"
  } else {
    title <- "Distribution of PRS of lung screening cohort patients"
  }
  ggplot(mediany, aes(x = MedSkore)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 60,
                   fill = "lightblue",
                   color = "black") +
    geom_density(aes(color = "PRS distribution"), linewidth = 0.65) +
    stat_function(fun = dnorm, args = list(mean = mu, sd = sigma),
                  aes(color = "Normal distribution"), linewidth = 0.6) +
    scale_color_manual(name = NULL,
                       values = c("PRS distribution" = "blue",
                                  "Normal distribution" = "black")) +
    annotate("text", x = Inf, y = Inf, label = normality_text,
             hjust = 1.1, vjust = 2, size = 4, color = "black") +
    labs(
      title = title,
      x = "Polygenic Risk Score",
      y = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_blank(), 
      legend.key = element_blank()          
    )
}

#vytvoreni rozlozeni percentilu
Percentile_distribution_plot <- function(data_trait,type){
  sample_counts <- table(data_trait$Sample)
  samples_with_3_or_more <- names(sample_counts[sample_counts >= 3])
  filtered_data <- data_trait[data_trait$Sample %in% samples_with_3_or_more, ]
  filtered_data[, Skore := sapply(Percentile, zpracuj_skore)]
  mediany <- filtered_data[, .(MedSkore = mean(Skore, na.rm = TRUE)), by = Sample]
  mu <- mean(mediany$MedSkore)
  sigma <- sd(mediany$MedSkore)
  
  
  shapiro_result <- shapiro.test(mediany$MedSkore)
  p_val <- shapiro_result$p.value
  normality_text <- if (p_val < 0.05) {
    paste0("Shapiro-Wilk p=", signif(p_val, 3), " -> NOT normal distribution")
  } else {
    paste0("Shapiro-Wilk p=", signif(p_val, 3), " -> normal distribution")
  }
  
  if (type == "h"){
    title <- "Distribution of Risk Percentile of normal cohort patients"
  } else {
    title <- "Distribution of Risk Percentile of lung screening cohort patients"
  }
  
  ggplot(mediany, aes(x = MedSkore)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 60,
                   fill = "lightblue",
                   color = "black") +
    geom_density(aes(color = "Risk Percentile distribution"), linewidth = 0.65) +
    stat_function(fun = dnorm, args = list(mean = mu, sd = sigma),
                  aes(color = "Normal distribution"), linewidth = 0.6) +
    scale_color_manual(name = NULL,
                       values = c("Risk Percentile distribution" = "blue",
                                  "Normal distribution" = "black")) +
    annotate("text", x = Inf, y = Inf, label = normality_text,
             hjust = 1.1, vjust = 2, size = 4, color = "black") +
    labs(
      title = title,
      x = "Risk Percentile",
      y = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = c(0.02, 0.98),   
      legend.justification = c(0, 1),
      legend.background = element_blank(),  
      legend.key = element_blank()          
    )
}

#qq plot
qqnorm(mediany$MedSkore)
qqline(mediany$MedSkore, col = "red", lwd = 2)


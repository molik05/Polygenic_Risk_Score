library(data.table)
library(ggplot2)
library(writexl)

#získání všech skore pro daný traint do 1 tabulky se všemi pacienty
soubory1 <- list.files(path = getwd(), pattern = "\\.tsv$", full.names = TRUE)
soubory2 <- list.files(path = "C:/Users/Admin/Documents/PolygenicRiskScore/411_results", pattern = "\\.tsv$", full.names = TRUE)
soubory <- c(soubory1,soubory2)
seznam_dt <- lapply(soubory, fread)
traits_zajem <- c(
  "Chronic Obstructive Pulmonary Disease",
  "Chronic Obstructive Pulmonary Disease (Moderate To Severe)",
  "Chronic Obstructive Pulmonary Disease (Severe)"
)

seznam_dt <- lapply(soubory, function(soubor) {
  dt <- fread(soubor)
  dt[Trait %in% traits_zajem]
})
names(seznam_dt) <- tools::file_path_sans_ext(basename(soubory))
vysledek <- rbindlist(seznam_dt, idcol = "soubor")

#ziskani jedne souhrne tabulky medianu pro každého pacient
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
vysledky <- lapply(soubory, function(soubor) {
  dt <- fread(soubor)
  
  # Ověření přítomnosti nutných sloupců
  if (!all(c("Trait", "Polygenic Risk Score") %in% names(dt))) {
    message("Přeskakuji soubor: ", basename(soubor), " (chybí sloupce)")
    return(NULL)
  }
  
  # Filtrování podle požadovaných traitů
  dt_filtr <- dt[Trait %in% traits_zajem]
  if (nrow(dt_filtr) == 0) return(NULL) 
  
  # Zpracování percentilů
  skore_cisla <- sapply(dt_filtr$Percentile, zpracuj_skore)
  prumer_perc <- median(skore_cisla, na.rm = TRUE)
  
  # Medián PRS (převedeme na čísla pro jistotu)
  dt_filtr[, `Polygenic Risk Score` := as.numeric(`Polygenic Risk Score`)]
  median_prs <- median(dt_filtr$`Polygenic Risk Score`, na.rm = TRUE)
  
  # Výstup pro jednoho pacienta
  data.table(
    Pacient = tools::file_path_sans_ext(basename(soubor)),
    MedianScore = prumer_perc,
    MedianPRS = median_prs
  )
})

# Spojení výsledků do jedné tabulky
vysledky <- rbindlist(Filter(Negate(is.null), vysledky))
# vysledky <- lapply(soubory, function(soubor) {
#   dt <- fread(soubor)
#   if (!all(c("Trait", "Polygenic Risk Score") %in% names(dt))) {
#     message("Přeskakuji soubor: ", basename(soubor), " (chybí sloupce)")
#     return(NULL)
#   }
#   dt_filtr <- dt[Trait %in% traits_zajem]
#   if (nrow(dt_filtr) == 0) return(NULL) 
#   skore_cisla <- sapply(dt_filtr$Percentile, zpracuj_skore)
#   prumer_perc <- median(skore_cisla, na.rm = TRUE)
#   data.table(
#     Pacient = tools::file_path_sans_ext(basename(soubor)),
#     MedianScore = prumer
#   )
# })
# vysledky <- rbindlist(Filter(Negate(is.null), vysledky))

write_xlsx(vysledek, path = "411_souhrna_tabulka_skore_pro_COPD.xlsx")


create_box_plots_pgs(getwd())
# make box plots for each trait and all patient in one run
create_box_plots_pgs <- function(path){
  slozka_vstup <- getwd()
  slozka_vystup <- file.path(slozka_vstup, "grafy_traity_skore")
  dir.create(slozka_vystup, showWarnings = FALSE)
  soubory <- list.files(path = slozka_vstup, pattern = "\\.tsv$", full.names = TRUE)
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
  vsechna_data <- lapply(soubory, function(soubor) {
    dt <- fread(soubor)
    if (!all(c("Trait", "Percentile") %in% names(dt))) {
      message("Přeskakuji soubor (chybí sloupce): ", basename(soubor))
      return(NULL)
    }
    dt[, Skore := sapply(`Polygenic Risk Score`, zpracuj_skore)]
    dt[, Pacient := tools::file_path_sans_ext(basename(soubor))]
    dt[, .(Pacient, Trait, Skore)]
  })
  vsechna_data <- rbindlist(Filter(Negate(is.null), vsechna_data))
  unique_traits <- unique(vsechna_data$Trait)
  for (trait in unique_traits) {
    data_trait <- vsechna_data[Trait == trait]
    mediany <- data_trait[, .(MedSkore = median(Skore, na.rm = TRUE)), by = Pacient]
    mediany[, Barva := fifelse(MedSkore <= 50, "0-50",
                               fifelse(MedSkore <= 80, "51-80", "81-100"))]
    data_trait <- merge(data_trait, mediany[, .(Pacient, Barva)], by = "Pacient", all.x = TRUE)
    p <- ggplot(data_trait, aes(x = Pacient, y = Skore, fill = Barva)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.7, color = "black") +
      #scale_fill_manual(values = c("0-50" = "forestgreen", "51-80" = "gold", "81-100" = "firebrick")) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(
        title = paste(trait),
        x = "Pacient",
        y = "Risk Score",
        fill = ""
      )
    nazev_souboru <- paste0(gsub("[^a-zA-Z0-9]", "_", trait), ".png")
    ggsave(filename = file.path(slozka_vystup, nazev_souboru), plot = p, width = 10, height = 6)
  }
}


median_list <- vysledky[, median(`Polygenic Risk Score`, na.rm = TRUE), by = Pacient]
setnames(median_list, "V1", "median_PRS")

# Sloučení s původní tabulkou podle Sample
v3 <- merge(vysledky, median_list, by = "Pacient", all.x = TRUE)

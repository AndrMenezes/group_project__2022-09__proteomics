library(dplyr)
library(readxl)
library(ggrepel)
library(cowplot)

path_data <- "./"

data_raw <- read_xlsx(
  file.path(path_data, "/Supplementary Data Stat Sigs.xlsx"),
  sheet = "All protein IDs",
  range = "A2:Y947")
head(data_raw)
glimpse(data_raw)
data_raw

data_raw <- data_raw |> 
  rename(
    protein__name = `Protein names`,
    protein__id = `Majority protein IDs`,
    gene__name = `Gene names`,
    peptides = Peptides,
    sequence_coverage = `Sequence coverage [%]`,
    molecular_weight = `Mol. weight [kDa]`,
    score = Score,
    intensity = Intensity,
    count = `MS/MS count`,
    sequence_length = `Sequence length`,
    lqf__Ampicillin__1 = `LFQ intensity Ampicillin Rep 1`,
    lqf__Ampicillin__2 = `LFQ intensity Ampicillin Rep 2`,
    lqf__Ampicillin__3 = `LFQ intensity Ampicillin Rep 3`,
    lqf__Cefotaxime__1 = `LFQ intensity Cefotaxime Rep 1`,
    lqf__Cefotaxime__2 = `LFQ intensity Cefotaxime Rep 2`,
    lqf__Cefotaxime__3 = `LFQ intensity Cefotaxime Rep 3`,
    lqf__Impipenem__1 = `LFQ intensity Impipenem Rep 1`,
    lqf__Impipenem__2 = `LFQ intensity Impipenem Rep 2`,
    lqf__Impipenem__3 = `LFQ intensity Impipenem Rep 3`,
    lqf__Ciprofloxacin__1 = `LFQ intensity Ciprofloxacin Rep 1`,
    lqf__Ciprofloxacin__2 = `LFQ intensity Ciprofloxacin Rep 2`,
    lqf__Ciprofloxacin__3 = `LFQ intensity Ciprofloxacin Rep 3`,
    lqf__Control__1 = `LFQ intensity Control Rep 1`,
    lqf__Control__2 = `LFQ intensity Control Rep 2`,
    lqf__Control__3 = `LFQ intensity Control Rep 3`,
  )


# Pivotting ---------------------------------------------------------------
data_pivotted <- data_raw |> 
  select(-c(peptides, protein__name, gene__name, sequence_coverage,
            molecular_weight, score, intensity, count, sequence_length)) |> 
  tidyr::pivot_longer(cols = -protein__id) |> 
  tidyr::separate(col = name, into = c("variable", "group", "replicate"),
                  sep = "__") |> 
  select(protein__id, group, replicate, variable, value)

# Dictionary --------------------------------------------------------------
data_dict <- data_raw |> 
  select(protein__id, protein__name, gene__name,
         peptides, sequence_coverage, molecular_weight, score, intensity,
         count, sequence_length)


# Exporting ---------------------------------------------------------------
write.csv(
  x = data_dict, file = file.path(path_data, "dictionary.csv"),
  row.names = FALSE)

write.csv(
  x = data_pivotted,
  file = file.path(path_data, "all_protein__pivotted.csv"),
  row.names = FALSE)

data <- read.csv("all_protein__pivotted.csv")

VolcanoPlot <- function(Treatment) {
data_summarized <- data |>
  filter(group %in% c("Control",Treatment)) |>
  select(protein__id, replicate, group, value)

### Control vs Ampicillin
if (Treatment == "Ampicillin") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(protein__id) |>
  summarise(
    logFC = mean(Control - Ampicillin),
    pvalue = t.test(Control, Ampicillin)$p.value,
    .groups = "drop"
 )
p1 <- ggplot(check, aes(logFC, -log(pvalue,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  ggtitle("Control Vs Ampicillin")
save_plot(file.path("./", "VolcanoPlotAMP.png"), p1)
return(p1)}

### Control vs Imipenem
if (Treatment == "Impipenem") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(protein__id) |>
  summarise(
    logFC = mean(Control - Impipenem),
    pvalue = t.test(Control, Impipenem)$p.value,
    .groups = "drop"
  )
p1 <- ggplot(check, aes(logFC, -log(pvalue,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  ggtitle("Control Vs Impipenem")
save_plot(file.path("./", "VolcanoPlotImi.png"), p1)
return(p1)
}

### Control vs Cefotaxime
if (Treatment == "Cefotaxime") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(protein__id) |>
  summarise(
    logFC = mean(Control - Cefotaxime),
    pvalue = t.test(Control, Cefotaxime)$p.value,
    .groups = "drop"
  )
p1 <- ggplot(check, aes(logFC, -log(pvalue,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  ggtitle("Control Vs Cefotaxime")
save_plot(file.path("./", "VolcanoPlotCefo.png"), p1)
return(p1)}

### Control vs Ciprofloxacin
if (Treatment == "Ciprofloxacin") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(protein__id) |>
  summarise(
    logFC = mean(Control - Ciprofloxacin),
    pvalue = t.test(Control, Ciprofloxacin)$p.value,
    .groups = "drop"
  )
p1 <- ggplot(check, aes(logFC, -log(pvalue,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  ggtitle("Control Vs Ciprofloxacin")
save_plot(file.path("./", "VolcanoPlotCipro.png"), p1)
return(p1)}}

## obtaining VolcanoPlots
VolcanoPlot("Ampicillin")
VolcanoPlot("Impipenem")
VolcanoPlot("Ciprofloxacin")
VolcanoPlot("Cefotaxime")
library(dplyr)
library(readxl)
library(ggrepel)
library(cowplot)

# Reading Data ------------------------------------------------------------
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
  select(-c(peptides, protein__name, sequence_coverage,
            molecular_weight, score, intensity, count, sequence_length)) |> 
  tidyr::pivot_longer(cols = -protein__id & -gene__name) |> 
  tidyr::separate(col = name, into = c("variable", "group", "replicate"),
                  sep = "__") |> 
  select(protein__id, gene__name, group, replicate, variable, value)

# Exporting ---------------------------------------------------------------
write.csv(
  x = data_pivotted,
  file = file.path(path_data, "all_protein__pivotted.csv"),
  row.names = FALSE)

data <- read.csv("all_protein__pivotted.csv")


# Volcano Plot -------------------------------------------------------------
VolcanoPlot <- function(Treatment) {
data_summarized <- data |>
  filter(group %in% c("Control",Treatment)) |>
  select(protein__id, gene__name, replicate, group, value)

## Control vs Ampicillin --------------------------------------------------
### logFC and pvalue 
if (Treatment == "Ampicillin") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(gene__name) |>
  summarise(
    logFC = mean(Ampicillin - Control),
    pvalue = t.test(Ampicillin, Control, var.equal = TRUE)$p.value,
    .groups = "drop"
 )
### highlight significantly regulated genes 
df <-  check %>%
  filter(abs(logFC) >= (1) & (pvalue) <= (0.05)) %>%
  head(20)

### Plotting 
p1 <- ggplot(check, aes((logFC), -log(pvalue,10)), options(ggrepel.max.overlaps = Inf)) + # -log10 conversion
  geom_point(size = 0.005, col="gray") +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  ggtitle("Control Vs Ampicillin") +
  geom_vline(xintercept=c(-0.5, 0.5), linetype = "dashed", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") +
  geom_point(data = df,
             shape = 21,
             size = 0.5,
             fill = "red",
             colour = "red") +
  geom_text_repel(data=df,
                  aes(x=logFC,
                      y=-log(pvalue,10),
                      label=gene__name),
                  size=1.5,
                  color="black") +
  theme_bw()

### Exporting
save_plot(file.path("./", "VolcanoPlot_Ampicillin.png"), p1)
return(p1)}

## Control vs Imipenem -------------------------------------------------------
### LogFC and pvalue
if (Treatment == "Impipenem") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(gene__name) |>
  summarise(
    logFC = mean(Impipenem - Control),
    pvalue = t.test(Impipenem, Control, var.equal = TRUE)$p.value,
    .groups = "drop"
  )
### highlight significantly regulated genes 
df <-  check %>%
  filter((abs(logFC) >= (1) & (pvalue) <= (0.05))) %>%
  head(20)

### plotting
p1 <- ggplot(check, aes((logFC), -log(pvalue,10)), options(ggrepel.max.overlaps = Inf)) + # -log10 conversion
  geom_point(size = 0.005, col="gray") +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  ggtitle("Control Vs Imipenem") +
  geom_vline(xintercept=c(-0.5, 0.5), linetype = "dashed", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") +
  geom_point(data = df,
             shape = 21,
             size = 0.5,
             fill = "red",
             colour = "red") +
  geom_text_repel(data=df,
                  aes(x=logFC,
                      y=-log(pvalue,10),
                      label=gene__name),
                  size=1.5,
                  color="black") +
  theme_bw()

### Exporting
save_plot(file.path("./", "VolcanoPlot_Imipenem.png"), p1)
return(p1)
}

## Control vs Cefotaxime -----------------------------------------------------
### LogFC and pvalue
if (Treatment == "Cefotaxime") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(gene__name) |>
  summarise(
    logFC = mean(Cefotaxime - Control),
    pvalue = t.test(Cefotaxime, Control, var.equal = TRUE)$p.value,
    .groups = "drop"
  )
### highlight significantly regulated genes 
df <-  check %>%
  filter((abs(logFC) >= (1) & (pvalue) <= (0.05))) %>%
  head(20)

### Plotting
p1 <- ggplot(check, aes((logFC), -log(pvalue,10)), options(ggrepel.max.overlaps = Inf)) + # -log10 conversion
  geom_point(size = 0.005, col="gray") +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  ggtitle("Control Vs Cefotaxime") +
  geom_vline(xintercept=c(-0.5, 0.5), linetype = "dashed", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") +
  geom_point(data = df,
             shape = 21,
             size = 0.5,
             fill = "red",
             colour = "red") +
  geom_text_repel(data=df,
                  aes(x=logFC,
                      y=-log(pvalue,10),
                      label=gene__name),
                  size=1.5,
                  color="black") +
  theme_bw()
### Exporting
save_plot(file.path("./", "VolcanoPlot_Cefotaxime.png"), p1)
return(p1)}

## Control vs Ciprofloxacin --------------------------------------------------
### LogFC and pvalue
if (Treatment == "Ciprofloxacin") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(gene__name) |>
  summarise(
    logFC = mean(Ciprofloxacin - Control),
    pvalue = t.test(Ciprofloxacin, Control, var.equal = TRUE)$p.value,
    .groups = "drop"
  )
### highlight significantly regulated genes 
df <-  check %>%
  filter((abs(logFC) >= (1) & (pvalue) <= (0.05))) %>%
  head(20)

### Plotting
p1 <- ggplot(check, aes((logFC), -log(pvalue,10)), options(ggrepel.max.overlaps = Inf)) + # -log10 conversion
  geom_point(size = 0.005, col="gray") +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"pvalue")) +
  ggtitle("Control Vs Ciprofloxacin") +
  geom_vline(xintercept=c(-0.5, 0.5), linetype = "dashed", col="black") +
  geom_hline(yintercept=-log10(0.05), linetype = "dashed", col="black") +
  geom_point(data = df,
             shape = 21,
             size = 0.5,
             fill = "red",
             colour = "red") +
  geom_text_repel(data=df,
                  aes(x=logFC,
                      y=-log(pvalue,10),
                      label=gene__name),
                  size=1.5,
                  color="black") +
  theme_bw()

### Exporting
save_plot(file.path("./", "VolcanoPlot_Ciprofloxacin.png"), p1)
return(p1)}}

# Generate and save Volcano Plots -------------------------------------------
VolcanoPlot("Ampicillin")
VolcanoPlot("Impipenem")
VolcanoPlot("Cefotaxime")
VolcanoPlot("Ciprofloxacin")


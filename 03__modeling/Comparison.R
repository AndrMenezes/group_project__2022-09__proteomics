library(dplyr)
library(readxl)

# Reading Data ------------------------------------------------------------
path_data <- "./"
data <- read.csv("all_protein__pivotted.csv")

## create a summary table --------------------------------------------------
tableSummary <- function(Treatment) {
  data_summarized <- data |>
    filter(group %in% c("Control",Treatment)) |>
    select(protein__id, gene__name, replicate, group, value)
  
  ## Control vs Ampicillin --------------------------------------------------
  ### logFC and pvalue 
  if (Treatment == "Ampicillin") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(protein__id, gene__name) |>
  summarise(
    group1="Control",
    group2=Treatment,
    logFC = mean(Ampicillin - Control),
    method="t.test",
    statistic = t.test(Ampicillin, Control, var.equal = TRUE)$statistic,
    pvalue = t.test(Ampicillin, Control, var.equal = TRUE)$p.value,
    .groups = "drop"
  )
FDR <- p.adjust(check$pvalue, method = "BH")
check <- cbind(check, FDR)
colnames(check)[1] <- "Protein"
check <- check[, c("group1", "group2", "Protein", "gene__name", "logFC", "method", "statistic", "pvalue", "FDR")]

## Exporting summary table -----------------------------------------------
write.csv(
  x = check,
  file = file.path(path_data, "Summary_Ampicillin.csv"),
  row.names = FALSE)
return(check)
}

## Control vs Imipenem --------------------------------------------------

### logFC and pvalue
if (Treatment == "Impipenem") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(protein__id, gene__name) |>
  summarise(
    group1="Control",
    group2=Treatment,
    logFC = mean(Impipenem - Control),
    method="t.test",
    statistic = t.test(Impipenem, Control, var.equal = TRUE)$statistic,
    pvalue = t.test(Impipenem, Control, var.equal = TRUE)$p.value,
    .groups = "drop"
  )
FDR <- p.adjust(check$pvalue, method = "BH")
check <- cbind(check, FDR)
colnames(check)[1] <- "Protein"
check <- check[, c("group1", "group2", "Protein", "gene__name", "logFC", "method", "statistic", "pvalue", "FDR")]
write.csv(
  x = check,
  file = file.path(path_data, "Summary_Imipenem.csv"),
  row.names = FALSE)
return(check)
}

## Control vs Cefotaxime --------------------------------------------------

### logFC and pvalue
if (Treatment == "Cefotaxime") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(protein__id, gene__name) |>
  summarise(
    group1="Control",
    group2=Treatment,
    logFC = mean(Cefotaxime - Control),
    method="t.test",
    statistic = t.test(Cefotaxime, Control, var.equal = TRUE)$statistic,
    pvalue = t.test(Cefotaxime, Control, var.equal = TRUE)$p.value,
    .groups = "drop"
  )
FDR <- p.adjust(check$pvalue, method = "BH")
check <- cbind(check, FDR)
colnames(check)[1] <- "Protein"
check <- check[, c("group1", "group2", "Protein", "gene__name", "logFC", "method", "statistic", "pvalue", "FDR")]

## Exporting summary table ---------------------------------------------------
write.csv(
  x = check,
  file = file.path(path_data, "Summary_Cefotaxime.csv"),
  row.names = FALSE)
return(check)
}
  
## Control vs Ciprofloxacin --------------------------------------------------

### logFC and pvalue
if (Treatment == "Ciprofloxacin") { check <- data_summarized |>
  tidyr::pivot_wider(names_from = group,
                     values_from = value) |>
  group_by(protein__id, gene__name) |>
  summarise(
    group1="Control",
    group2=Treatment,
    logFC = mean(Ciprofloxacin - Control),
    method="t.test",
    statistic = t.test(Ciprofloxacin, Control, var.equal = TRUE)$statistic,
    pvalue = t.test(Ciprofloxacin, Control, var.equal = TRUE)$p.value,
    .groups = "drop"
  )
FDR <- p.adjust(check$pvalue, method = "BH")
check <- cbind(check, FDR)
colnames(check)[1] <- "Protein"
check <- check[, c("group1", "group2", "Protein", "gene__name", "logFC", "method", "statistic", "pvalue", "FDR")]

## Exporting summary table -----------------------------------------------
write.csv(
  x = check,
  file = file.path(path_data, "Summary_Ciprofloxacin.csv"),
  row.names = FALSE)
return(check)
}}

tableSummary("Ampicillin")
tableSummary("Impipenem")
tableSummary("Cefotaxime")
tableSummary("Ciprofloxacin")
-----------------------------------------------------------------------------

### now we can perform different analyses
  # uploading summary table
Ampicillin <- read.csv("Summary_Ampicillin.csv")
Imipenem <- read.csv("Summary_Imipenem.csv")
Cefotaxime <- read.csv("Summary_Cefotaxime.csv")
Ciprofloxacin <- read.csv("Summary_Ciprofloxacin.csv")

# ----USAGE: to perform the selected analysis -> corresponding function(treatment)
#   in which treatment= Ampicillin|Impipenem|Cefotaxime|Ciprofloxacin------------

### perform traditional analysis with |logFC| > 0.5 & pvalue < 0.05
P <- function(treatment) {
df1 <-  treatment %>%
  filter(((logFC) >= (0.05)) & (pvalue) <= (0.05)) %>%
  arrange(desc(logFC)) %>%
  head(10)

df2 <-  treatment %>%
  filter(((logFC) <= -(0.05)) & (pvalue) <= (0.05)) %>%
  arrange(logFC) %>%
  head(10)

df <- rbind(df1, df2[ ,])
return(df)}

### Perform analysis pvalue based, without a cutoff on LogFC
Pbased <- function(treatment) {
  df <-  treatment %>%
    filter((pvalue) <= (0.05)) %>%
    arrange(pvalue) %>%
    head(20)
  return(df)}

### perform more restrictive analysis with |logFC| > 0.5 & pvalue < 0.005
newP <- function(treatment) {
df1 <- treatment  %>%
  filter((logFC) >= (0.5) & (pvalue) <= (0.005)) %>%
  arrange(desc(logFC)) %>%
  head(10)

df2 <-  treatment %>%
  filter(((logFC) <= -(0.5)) & (pvalue) <= (0.005)) %>%
  arrange(logFC) %>%
  head(10)

df <- rbind(df1, df2[ ,])
return(df)}

### Perform analysis based on more restrictive Pvalue, without a cutoff on LogFC
newPbased <- function(treatment) {
  df <-  treatment %>%
    filter((pvalue) <= (0.005)) %>%
    arrange(pvalue) %>%
    head(20)
  return(df)}

### perform analysis with |logFC| > 0.5 & p value adjustment (FDR < 0.1)
FDR <- function(treatment) {
df1 <-  treatment %>%
  filter(((logFC) >= (0.5)) & (FDR) <= (0.1)) %>%
  arrange(desc(logFC)) %>%
  head(10)

df2 <-  treatment %>%
  filter(((logFC) <= -(0.5)) & (FDR) <= (0.1)) %>%
  arrange(logFC) %>%
  head(10)

df <- rbind(df1, df2[ ,])
return(df)}

### Perform analysis FDR based, without a cutoff on LogFC
FDRbased <- function(treatment) {
df <-  treatment %>%
  filter((FDR) <= (0.1)) %>%
  arrange(FDR) %>%
  head(20)
return(df)}

## to see same genes
subset(newP(Ciprofloxacin)$gene__name, newP(Ciprofloxacin)$gene__name == FDR(Ciprofloxacin)$gene__name)

library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(tidyverse)
library(openxlsx)

# Reading Data ------------------------------------------------------------
path_data <- "./"
data <- read_xlsx(file.path(path_data, "all_results.xlsx"))


# ----USAGE: to perform the selected analysis -> corresponding function(treatment)
#   in which treatment= "Ampicillin"|"Impipenem"|"Cefotaxime"|"Ciprofloxacin"------------

### perform traditional analysis with |logFC| > 0.5 & pvalue < 0.05
P <- function(treatment, name="no", method="Two Sample T Test") {
  data_filtered <- data |>
    filter(Method==method, Treatment %in% c(treatment)) |>
    select(c(Protein, Gene, logFC, P_value))

df1 <- data_filtered %>%
  filter(((logFC) >= (0.5)) & ((P_value) <= (0.05))) %>%
  arrange(desc(logFC))
  if (name=="yes") {df1 <- df1 %>% head(10)}

df2 <-  data_filtered %>%
  filter(((logFC) <= -(0.5)) & ((P_value) <= (0.05))) %>%
  arrange(logFC)
  if (name=="yes") {df2 <- df2 %>% head(10)}

df <- rbind(df1, df2[ ,])
return(df)}

### Perform analysis pvalue based, without a cutoff on LogFC
Pbased <- function(treatment, name="no", method="Two Sample T Test") {
  data_filtered <- data |>
    filter(Method==method, Treatment %in% c(treatment)) |>
    select(c(Protein, Gene, logFC, P_value))
  
  df <-  data_filtered %>%
    filter((P_value) <= (0.05)) %>%
    arrange(P_value)
  if (name=="yes") {df <- df %>% head(20)}
  return(df)}

### perform more restrictive analysis with |logFC| > 0.5 & pvalue < 0.005
newP <- function(treatment, name="no", method="Two Sample T Test") {
  data_filtered <- data |>
    filter(Method==method, Treatment %in% c(treatment)) |>
    select(c(Protein, Gene, logFC, P_value))
  
df1 <- data_filtered  %>%
  filter((logFC) >= (0.5) & (P_value) <= (0.005)) %>%
  arrange(desc(logFC))
if (name=="yes") {df1 <- df1 %>% head(10)}

df2 <-  data_filtered %>%
  filter(((logFC) <= -(0.5)) & (P_value) <= (0.005)) %>%
  arrange(logFC) 
if (name=="yes") {df2 <- df2 %>% head(10)}

df <- rbind(df1, df2[ ,])
return(df)}

### Perform analysis based on more restrictive Pvalue, without a cutoff on LogFC
newPbased <- function(treatment, name="no", method="Two Sample T Test") {
  data_filtered <- data |>
    filter(Method==method, Treatment %in% c(treatment)) |>
    select(c(Protein, Gene, Method, logFC, P_value))
  
  df <-  data_filtered %>%
    filter((P_value) <= (0.005)) %>%
    arrange(P_value) 
  if (name=="yes") {df <- df %>% head(20)}
  return(df)}

### perform analysis with |logFC| > 0.5 & p value adjustment (FDR < 0.1)
FDR <- function(treatment, name="no", method="Two Sample T Test") {
  data_filtered <- data |>
    filter(Method==method, Treatment %in% c(treatment)) |>
    select(c(Protein, Gene, logFC, FDR))
  
df1 <-  data_filtered %>%
  filter((logFC) >= (0.5) & (FDR) <= (0.1)) %>%
  arrange(desc(logFC))
if (name=="yes") {df1 <- df1 %>% head(10)}

df2 <-  data_filtered %>%
  filter(((logFC) <= -(0.5)) & (FDR) <= (0.1)) %>%
  arrange(logFC)
if (name=="yes") {df2 <- df2 %>% head(10)}

df <- rbind(df1, df2[ ,])
return(df)}

### Perform analysis FDR based, without a cutoff on LogFC
FDRbased <- function(treatment, name="no", method="Two Sample T Test") {
  data_filtered <- data |>
    filter(Method==method, Treatment %in% c(treatment)) |>
    select(c(Protein, Gene, logFC, FDR))
df <-  data_filtered %>%
  filter((FDR) <= (0.1)) %>%
  arrange(FDR) 
if (name=="yes") {df <- df %>% head(20)}
return(df)}

## Protein Differential Expressed in each Treatment based on Different Approaches - Bar Chart

# Ampicillin -----------------------------------------------------------------------------------

Ampicillin_proteinDE <- tibble(nrow(P("Ampicillin")), nrow(Pbased("Ampicillin")), nrow(newP("Ampicillin")), nrow(newPbased("Ampicillin")), nrow(FDR("Ampicillin")), nrow(FDRbased("Ampicillin")))
Ampicillin <- as.numeric(Ampicillin_proteinDE)
Approaches <- c("P", "Pbased", "newP", "newPbased", "FDR", "FDRbased")
Ampicillin_proteinDE <- cbind(Approaches, Ampicillin)
Ampicillin_proteinDE <- as.data.frame(Ampicillin_proteinDE)

# summary table
Pdrug <- P("Ampicillin")
P_based_drug <- Pbased("Ampicillin")
newP_drug <- newP("Ampicillin")
new_Pbased_drug <- newPbased("Ampicillin")
FDR_drug<- FDR("Ampicillin")
FDR_based_drug <- FDRbased("Ampicillin")
# export xlsx
list_amp <- list("P" = Pdrug, "Pbased" = P_based_drug, "new_P" = newP_drug, "newP_based" = new_Pbased_drug, "FDR" = FDR_drug, "FDR_based" = FDR_based_drug)
                 write.xlsx(list_amp, file = "Approaches_Ampiccilin.xlsx")
# plotting
p <- ggplot(Ampicillin_proteinDE, aes(x=Approaches, y=as.numeric(Ampicillin))) + 
  geom_bar(position = position_dodge(), stat="identity", width=.5, fill="tomato3") + 
  ylab(expression("# Protein DE")) +
  geom_text(aes(label = Ampicillin), position=position_dodge(width=0.9), size=4) +
  labs(title="Ampicillin Bar Chart", 
       subtitle="Approaches comparison") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
p
## export plot
save_plot(file.path("./", "Ampicillin_BarChart.png"), p)

# Imipenem -----------------------------------------------------------------------------------

Impipenem_proteinDE <- tibble(nrow(P("Impipenem")), nrow(Pbased("Impipenem")), nrow(newP("Impipenem")), nrow(newPbased("Impipenem")), nrow(FDR("Impipenem")), nrow(FDRbased("Impipenem")))
Imipenem <- as.numeric(Impipenem_proteinDE)
Aproaches <- c("P", "Pbased", "newP", "newPbased", "FDR", "FDRbased")
Impipenem_proteinDE <- cbind(Approaches, Imipenem)
Impipenem_proteinDE <- as.data.frame(Impipenem_proteinDE)

# summary table
Pdrug <- P("Impipenem")
P_based_drug <- Pbased("Impipenem")
newP_drug <- newP("Impipenem")
new_Pbased_drug <- newPbased("Impipenem")
FDR_drug<- FDR("Impipenem")
FDR_based_drug <- FDRbased("Impipenem")
# export xlsx
list_imi <- list("P" = Pdrug, "Pbased" = P_based_drug, "new_P" = newP_drug, "newP_based" = new_Pbased_drug, "FDR" = FDR_drug, "FDR_based" = FDR_based_drug)
write.xlsx(list_imi, file = "Approaches_Imipenem.xlsx")

# plotting

p <- ggplot(Impipenem_proteinDE, aes(x=Approaches, y=as.numeric(Imipenem))) + 
  geom_bar(position = position_dodge(), stat="identity", width=.5, fill="yellow") + 
  ylab(expression("# Protein DE")) +
  geom_text(aes(label = Imipenem), position=position_dodge(width=0.9), size=4) +
  labs(title="Imipenem Bar Chart", 
       subtitle="Approaches comparison") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
p
## export plot
save_plot(file.path("./", "Imipenem_BarChart.png"), p)

# Cefotaxime -----------------------------------------------------------------------------------

Cefotaxime_proteinDE <- tibble(nrow(P("Cefotaxime")), nrow(Pbased("Cefotaxime")), nrow(newP("Cefotaxime")), nrow(newPbased("Cefotaxime")), nrow(FDR("Cefotaxime")), nrow(FDRbased("Cefotaxime")))
Cefotaxime <- as.numeric(Cefotaxime_proteinDE)
Approaches <- c("P", "Pbased", "newP", "newPbased", "FDR", "FDRbased")
Cefotaxime_proteinDE <- cbind(Approaches, Cefotaxime)
Cefotaxime_proteinDE <- as.data.frame(Cefotaxime_proteinDE)

# summary table
Pdrug <- P("Cefotaxime")
P_based_drug <- Pbased("Cefotaxime")
newP_drug <- newP("Cefotaxime")
new_Pbased_drug <- newPbased("Cefotaxime")
FDR_drug<- FDR("Cefotaxime")
FDR_based_drug <- FDRbased("Cefotaxime")
# export xlsx
list_cef <- list("P" = Pdrug, "Pbased" = P_based_drug, "new_P" = newP_drug, "newP_based" = new_Pbased_drug, "FDR" = FDR_drug, "FDR_based" = FDR_based_drug)
write.xlsx(list_cef, file = "Approaches_Cefotaxime.xlsx")

# plotting 

p <- ggplot(Cefotaxime_proteinDE, aes(x=Approaches, y=as.numeric(Cefotaxime))) + 
  geom_bar(position = position_dodge(), stat="identity", width=.5, fill="green") + 
  ylab(expression("# Protein DE")) +
  geom_text(aes(label = Cefotaxime), position=position_dodge(width=0.9), size=4) +
  labs(title="Cefotaxime Bar Chart", 
       subtitle="Approaches comparison") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
p
## export plot
save_plot(file.path("./", "Cefotaxime_BarChart.png"), p)

# Ciprofloxacin -----------------------------------------------------------------------------------

Ciprofloxacin_proteinDE <- tibble(nrow(P("Ciprofloxacin")), nrow(Pbased("Ciprofloxacin")), nrow(newP("Ciprofloxacin")), nrow(newPbased("Ciprofloxacin")), nrow(FDR("Ciprofloxacin")), nrow(FDRbased("Ciprofloxacin")))
Ciprofloxacin <- as.numeric(Ciprofloxacin_proteinDE)
Approaches <- c("P", "Pbased", "newP", "newPbased", "FDR", "FDRbased")
Ciprofloxacin_proteinDE <- cbind(Approaches, Ciprofloxacin)
Ciprofloxacin_proteinDE <- as.data.frame(Ciprofloxacin_proteinDE)

# summary table
Pdrug <- P("Ciprofloxacin")
P_based_drug <- Pbased("Ciprofloxacin")
newP_drug <- newP("Ciprofloxacin")
new_Pbased_drug <- newPbased("Ciprofloxacin")
FDR_drug<- FDR("Ciprofloxacin")
FDR_based_drug <- FDRbased("Ciprofloxacin")
# export xlsx
list_cip <- list("P" = Pdrug, "Pbased" = P_based_drug, "new_P" = newP_drug, "newP_based" = new_Pbased_drug, "FDR" = FDR_drug, "FDR_based" = FDR_based_drug)
write.xlsx(list_cip, file = "Approaches_Ciprofloxacin.xlsx")

# plotting

p <- ggplot(Ciprofloxacin_proteinDE, aes(x=Approaches, y=as.numeric(Ciprofloxacin))) + 
  geom_bar(position = position_dodge(), stat="identity", width=.5, fill="blue") + 
  ylab(expression("# Protein DE")) +
  geom_text(aes(label = Ciprofloxacin), position=position_dodge(width=0.9), size=4)
  labs(title="Ciprofloxacin Bar Chart", 
       subtitle="Approaches comparison") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))
p
## export plot
save_plot(file.path("./", "Ciprofloxacin_BarChart.png"), p)

# Summary Bar Chart of how many differential expressed protein each approaches identify for each treatment

## preparing table of data
sumupCC <- right_join(Cefotaxime_proteinDE, Ciprofloxacin_proteinDE, by="Approaches")
sumupAI <- right_join(Ampicillin_proteinDE, Impipenem_proteinDE, by="Approaches")
sumup <- right_join(sumupAI, sumupCC, by="Approaches")
graph <- sumup %>% pivot_longer(!Approaches, names_to = 'Treatment')
head(graph)

## plotting
tot <- graph %>%
  ggplot(aes(x = Approaches, y = as.numeric(value), 
             fill = Treatment)) +
  ylab(expression("# Protein DE")) +
  geom_bar(position = position_dodge(), stat="identity") +
  geom_text(aes(label = value), position=position_dodge(width=0.9), size=4)
tot
save_plot(file.path("./", "ProteinDE_BarChart.png"), tot)

# see topGene | choose approaches and treatment
FDR("Ciprofloxacin", "yes")$Gene
P("Ampicillin", "yes")$Gene

# see topGene | choose approaches and treatment
P("Ampicillin", "yes")$Protein

# Comparison between two selected approaches -- Common Genes
intersect((P("Ciprofloxacin", "yes"))$Gene, (FDR("Ciprofloxacin", "yes")$Gene))
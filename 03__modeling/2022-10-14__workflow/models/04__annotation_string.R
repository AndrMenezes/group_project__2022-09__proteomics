library(STRINGdb)
library(openxlsx)
library(readr)
library(ggplot2)
library(dplyr)

data_de <- read.csv("de_proteins.csv")
head(data_de)
dim(data_de)
data_de
data_filtered <- data_de |> 
  dplyr::filter(p_value <= 0.005 & abs(log_fc) >= 0.5) |> 
  dplyr::as_tibble()

data_filtered

data_de |> 
  dplyr::filter(p_value <= 0.005 & abs(log_fc) >= 0.5) |> 
  dplyr::group_by(method_label, group) |> 
  dplyr::count() |> 
  ggplot(aes(x = group, y = n, fill = method_label)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.6,
           col = "black") +
  geom_text(aes(label = n), vjust = -0.05,
            position = position_dodge(width = 0.9), size = 6) +
  labs(x = "", y = "", fill = "") +
  scale_y_continuous(limits = c(0, 55)) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank())

# string_db <- STRINGdb$new( version="11.5", species=9606,
#                            score_threshold=200, network_type="functional", input_directory="./")
# 
# example1_mapped <- string_db$map( Ampicillin, "protein", removeUnmappedRows = TRUE )
# string_proteins <- string_db$get_proteins()
# tp53 = string_db$get_proteins()


## Ampicillin

Proteins_Ampicillin <- subset(data_filtered, data_filtered$group == "Ampicillin")

# upload this file on https://string-db.org/ and download mapping
write.table(
  Proteins_Ampicillin$protein,
  file = "Proteins_Ampicillin.csv",
  row.names = FALSE,
  col.names = FALSE)

uniprotAMP <- read_tsv("string_mappingAMP.tsv")
uniprot_dataAMP <- uniprotAMP |> dplyr::select(c(preferredName, annotation)) |>
  dplyr::rename(gene = preferredName)

amp <- cbind(Proteins_Ampicillin[!grepl("ACQ", Proteins_Ampicillin$protein , ignore.case = FALSE), ], uniprot_dataAMP) # STRING doesn't recognize proteins called "ACQ..."
ampACQ <- cbind(Proteins_Ampicillin[grepl("ACQ", Proteins_Ampicillin$protein , ignore.case = FALSE), ], gene = "", annotation = "") # don't lose any proteins
ampALL <- rbind(ampACQ, amp)
# export table
write.table(
  ampALL,
  file = "Proteins_AmpicillinSTRING.tsv",
  row.names = FALSE)

-------------------------------------------------------------------------------------------------------------
## Cefotaxime

Proteins_Cefotaxime <- subset(data_filtered, data_filtered$group == "Cefotaxime")

# upload this file on https://string-db.org/ and download mapping
write.table(
  Proteins_Cefotaxime$protein,
  file = "Proteins_Cefotaxime.csv",
  row.names = FALSE,
  col.names = FALSE)

uniprotCEF <- read_tsv("string_mappingCEF.tsv")
uniprot_dataCEF <- uniprotCEF |> dplyr::select(c(preferredName, annotation)) |>
  dplyr::rename(gene = preferredName)

cef <- cbind(Proteins_Cefotaxime[!grepl("ACQ", Proteins_Cefotaxime$protein , ignore.case = FALSE), ], uniprot_dataCEF) # STRING doesn't recognize proteins called "ACQ..."
write.table(
  cef,
  file = "Proteins_CefotaximeSTRING.tsv",
  row.names = FALSE)
# cefACQ <- cbind(Proteins_Cefotaxime[grepl("ACQ", Proteins_Cefotaxime$protein , ignore.case = FALSE), ], gene = "", annotation = "") # don't lose any proteins
# cefALL <- rbind(cefACQ, cef)
# 
# # export table
# write.table(
#   cefALL,
#   file = "Proteins_CefotaximeSTRING.tsv",
#   row.names = FALSE)

-------------------------------------------------------------------------------------------------------------
## Ciprofloxacin

Proteins_Ciprofloxacin <- subset(data_filtered, data_filtered$group == "Ciprofloxacin")

# upload this file on https://string-db.org/ and download mapping
write.table(
  Proteins_Ciprofloxacin$protein,
  file = "Proteins_Ciprofloxacin.csv",
  row.names = FALSE,
  col.names = FALSE)

uniprotCIP <- read_tsv("string_mappingCIP.tsv")
uniprot_dataCIP <- uniprotCIP |> dplyr::select(c(preferredName, annotation)) |>
  dplyr::rename(gene = preferredName)

cip <- cbind(Proteins_Ciprofloxacin[!grepl("ACQ", Proteins_Ciprofloxacin$protein , ignore.case = FALSE), ], uniprot_dataCIP) # STRING doesn't recognize proteins called "ACQ..."
cipACQ <- cbind(Proteins_Ciprofloxacin[grepl("ACQ", Proteins_Ciprofloxacin$protein , ignore.case = FALSE), ], gene = "", annotation = "") # don't lose any proteins
cipALL <- rbind(cipACQ, cip)
# export table
write.table(
  cipALL,
  file = "Proteins_CiprofloxacinSTRING.tsv",
  row.names = FALSE)

------------------------------------------------------------------------------------------------------------------------------------------------
## Impipenem

Proteins_Impipenem <- subset(data_filtered, data_filtered$group == "Impipenem")

# upload this file on https://string-db.org/ and download mapping
write.table(
  Proteins_Impipenem$protein,
  file = "Proteins_Impipenem.csv",
  row.names = FALSE,
  col.names = FALSE)

uniprotIMI <- read_tsv("string_mappingIMI.tsv")
uniprot_dataIMI <- uniprotIMI |> dplyr::select(c(preferredName, annotation)) |>
  dplyr::rename(gene = preferredName)

imi <- cbind(Proteins_Impipenem[!grepl("ACQ", Proteins_Impipenem$protein , ignore.case = FALSE), ], uniprot_dataIMI) # STRING doesn't recognize proteins called "ACQ..."
imiACQ <- cbind(Proteins_Impipenem[grepl("ACQ", Proteins_Impipenem$protein , ignore.case = FALSE), ], gene = "", annotation = "") # don't lose any proteins
imiALL <- rbind(imiACQ, imi)
# export table
write.table(
  imiALL,
  file = "Proteins_ImpipenemSTRING.tsv",
  row.names = FALSE)


# one file
allProteins <- rbind(ampALL, cef, cipALL, imiALL)
write.table(
  allProteins,
  file = "Proteins_STRING.tsv",
  row.names = FALSE)

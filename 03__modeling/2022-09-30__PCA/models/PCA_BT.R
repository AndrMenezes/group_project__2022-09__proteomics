library(dplyr)
library(readxl)
library(ggplot2)
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
    lqf_Ampicillin_1 = `LFQ intensity Ampicillin Rep 1`,
    lqf_Ampicillin_2 = `LFQ intensity Ampicillin Rep 2`,
    lqf_Ampicillin_3 = `LFQ intensity Ampicillin Rep 3`,
    lqf_Cefotaxime_1 = `LFQ intensity Cefotaxime Rep 1`,
    lqf_Cefotaxime_2 = `LFQ intensity Cefotaxime Rep 2`,
    lqf_Cefotaxime_3 = `LFQ intensity Cefotaxime Rep 3`,
    lqf_Impipenem_1 = `LFQ intensity Impipenem Rep 1`,
    lqf_Impipenem_2 = `LFQ intensity Impipenem Rep 2`,
    lqf_Impipenem_3 = `LFQ intensity Impipenem Rep 3`,
    lqf_Ciprofloxacin_1 = `LFQ intensity Ciprofloxacin Rep 1`,
    lqf_Ciprofloxacin_2 = `LFQ intensity Ciprofloxacin Rep 2`,
    lqf_Ciprofloxacin_3 = `LFQ intensity Ciprofloxacin Rep 3`,
    lqf_Control_1 = `LFQ intensity Control Rep 1`,
    lqf_Control_2 = `LFQ intensity Control Rep 2`,
    lqf_Control_3 = `LFQ intensity Control Rep 3`,
  )

data_processed <- data_raw |>
  select(
    protein__name,
    gene__name,
    lqf_Ampicillin_1,
    lqf_Ampicillin_2,
    lqf_Ampicillin_3,
    lqf_Cefotaxime_1,
    lqf_Cefotaxime_2,
    lqf_Cefotaxime_3,
    lqf_Impipenem_1,
    lqf_Impipenem_2,
    lqf_Impipenem_3,
    lqf_Ciprofloxacin_1,
    lqf_Ciprofloxacin_2,
    lqf_Ciprofloxacin_3,
    lqf_Control_1,
    lqf_Control_2,
    lqf_Control_3
  )
# Exporting -------------------------------------------------------------------
write.csv(
  x = data_processed,
  file = file.path(path_data, "data_processed.csv"),
  row.names = FALSE)

# PC analysis -----------------------------------------------------------------
data_processed <- read.csv("data_processed.csv")
data <- subset(data_processed, select = -c(protein__name,gene__name))
data <- t(data)
pca <- prcomp(data)

# PC variance -----------------------------------------------------------------
pca_var <- pca$sdev^2
pca_var_percent <- round(pca_var / sum(pca_var) * 100, 1)

# Plotting --------------------------------------------------------------------
pca_plot <- data.frame(pca$x[, 1], pca$x[, 2])

## create a treatment subset to colour and shape the plot
treatment <- (c(row.names(pca_plot)))
treatment[1:15] <- rep(c('Ampicillin', 'Cefotaxime', 'Imipenem', 'Ciprofloxacine', 'Control'), each  = 3)
treatment <- as.factor(treatment)

## create a dataframe t plot
pca_data <- data.frame(treatment, 
                       X = pca_plot$pca.x...1., 
                       Y = pca_plot$pca.x...2.)
## PCA Plot
plot <- ggplot(pca_data, aes(X, Y)) +
                 geom_point(aes(colour=treatment, shape=treatment, fill=treatment), size = 5, colour="black") +
                 scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
                 scale_fill_manual(values = c( "red", "blue", "green", "purple", "yellow")) +
                 xlab(paste("Component 1 (", pca_var_percent[1], "%)", sep = "")) +
                 ylab(paste("Component 2 (", pca_var_percent[2], "%)", sep = "")) +
                 theme_bw() +
                 ggtitle("PCA Analysis")

# Exporting --------------------------------------------------------------------
save_plot(file.path("./", "PCA.png"), plot)

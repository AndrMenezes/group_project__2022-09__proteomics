Rscript -e  'renv::run("./03__modeling/2022-10-14__workflow/processing/00__imputation.R")'
Rscript -e  'renv::run("./03__modeling/2022-10-14__workflow/processing/01__normalization.R")'
Rscript -e  'renv::run("./03__modeling/2022-10-14__workflow/models/01__scran_feature_selection.R")'
Rscript -e  'renv::run("./03__modeling/2022-10-14__workflow/models/02__effect_size_feature_selection.R")'

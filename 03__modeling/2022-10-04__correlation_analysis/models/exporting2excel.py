import pandas as pd
from pandas import ExcelWriter

path_local = "./03__modeling/2022-10-04__correlation_analysis/%s"
path_res = path_local % "models/results/%s"

data = pd.read_csv(path_res % "estimated_correlation.csv")
del data["correlation"]

xlsx_file = path_res % "estimated_correlation.xlsx"
with ExcelWriter(xlsx_file) as writer:
    for group, data_curr in data.groupby("group"):
        print(group)
        data_curr.to_excel(writer, group, index=False)

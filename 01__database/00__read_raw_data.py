import re
import pandas as pd


# Import all peptides data
raw_data = pd.read_table("./01__database/raw_data/allPeptides.txt")

# Selecting the important columns
raw_data = raw_data[["Raw file", "Sequence", "Proteins", "Intensity"]].copy()

# Filtering empty sequences and proteins
raw_data = raw_data[
    (raw_data["Sequence"] != " ") & (raw_data["Proteins"] != " ")].copy()

lt_to_merge = []
cols = ["Sequence", "Proteins", "Intensity"]
for raw_file in raw_data["Raw file"].unique():
    print(raw_file)
    new_name = "Intensity." + re.sub("Ecoli_*", "", raw_file)
    tmp = raw_data[raw_data["Raw file"] == raw_file][cols].copy().rename(
        columns={"Intensity": new_name})
    tmp.set_index(["Sequence", "Proteins"], inplace=True)
    lt_to_merge.append(tmp)

merged_data = pd.concat(lt_to_merge)

merged_data.isna().sum(axis=1)
merged_data.isna().sum(axis=0) / merged_data.shape[0]
raw_data[raw_data["Proteins"] == "P0AFG6"]
merged_data_reseted = merged_data.reset_index()
chupa = merged_data_reseted[merged_data_reseted["Proteins"] == "P0AFG6"]
chupa["Intensity.Amp2"].isna().sum()

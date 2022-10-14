import pandas as pd

data = pd.read_table("./01__database/raw_data/allPeptides.txt")
data.iloc[0, :]
data["Raw file"].value_counts()

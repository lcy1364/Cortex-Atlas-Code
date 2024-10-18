# ls Supercluster*h5ad |xargs -I {} -P 24 python3 2_splitBySubclass.py {}
import scanpy as sc
import numpy as np
import sys

#np.random.seed(32)
file = sys.argv[1]
print(file)

adata = sc.read_h5ad(file)
adata = adata[adata.obs["Region"] != "MTG", :]
adata
# Column to split by
split_column = 'CrossArea_subclass'

# Get unique values in the column
unique_values = adata.obs[split_column].unique()

# Split the AnnData object and save each subset
for val in unique_values:
    subset = adata[adata.obs[split_column] == val, :]
    subset
    filename = f"{val}_subclass.h5ad".replace("/","_")
    print(f"Saved {val} subset to {filename}")
    subset.write_h5ad(filename)
    print("done")

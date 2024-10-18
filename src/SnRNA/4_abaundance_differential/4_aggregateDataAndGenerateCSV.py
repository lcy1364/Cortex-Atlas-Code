# sbatch 9_2_batch2
import sys
import os
import time
import random
import anndata as ad
import pickle
import pandas as pd
from datetime import datetime

pwd = os.getcwd()
# Print all command line arguments
for i, arg in enumerate(sys.argv):
    print(f"Argument {i}: {arg}",flush = True)

# Print all command line arguments
random_sleep_time = random.uniform(1,5)

print("Start : %s" % time.ctime(),flush = True)
time.sleep( random_sleep_time )
print("End : %s" % time.ctime(),flush = True)

region = "cross_areal"
k = "region"
l = sys.argv[2] 
abundances = ad.read_h5ad(sys.argv[1]) 

print(abundances, flush = True)
outfile = sys.argv[1].replace(".h5ad","") + "_results.csv"
print(outfile, flush = True)
print("process " + k + ": " + l, flush = True)
cell_types = abundances.var.index
print(os.path.join(pwd, "output", region, k, l, outfile))

if not os.path.exists(os.path.join(pwd, "output", region, k, l, outfile)):
    
    while i < cell_types.shape[0]:    
    
        ct = cell_types[i]
        print(str(datetime.now()) + " process " + l + " " + ct, flush = True)
    
        pickle_file = os.path.join(pwd, "output", region, k, l, ct.replace("/", " ") + ".pkl")
        with open(pickle_file, 'rb') as pickle_load:
            results = pickle.load(pickle_load)
    
        accepted = results.sample_stats["is_accepted"].to_numpy()
    
        if accepted.sum() / accepted.shape[1] < 0.6:
            print(l + " " + ct + " Did not achieve acceptance threshold")
            continue
    
        else:
            print(str(datetime.now()) + " Converged " + l + " " + ct, flush = True)
            i = i + 1
    
        current_results = results.summary_prepare()[1]
        current_results = current_results.reset_index()
        current_results["Reference Cell Type"] = ct
    
        try:
            results_table = pd.concat([results_table, current_results], axis=0)
        except:
            results_table = current_results.copy()
    
    results_table.to_csv(os.path.join(pwd, "output", region, k, l, outfile)) 
    print(str(datetime.now()) + l + "done", flush = True)
    
else:
    print("done and skip")
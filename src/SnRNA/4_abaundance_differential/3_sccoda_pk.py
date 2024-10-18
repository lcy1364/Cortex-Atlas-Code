import sys
import os
import time
import random
pwd = os.getcwd()
# Print all command line arguments
for i, arg in enumerate(sys.argv):
    print(f"Argument {i}: {arg}",flush = True)

# Print all command line arguments
random_sleep_time = random.uniform(5,10)

print("Start : %s" % time.ctime(),flush = True)
time.sleep( random_sleep_time )
print("End : %s" % time.ctime(),flush = True)

region = "cross_areal"
k = "region"
l = sys.argv[2]
ct = sys.argv[3].replace("="," ")
jobNum = sys.argv[4]
if os.path.exists(os.path.join(pwd, "output", region, k, l)) is False:
    os.makedirs(os.path.join(pwd, "output", region, k, l))
    
pickle_file = os.path.join(pwd, "output", region, k, l, ct.replace("/", " ") + ".pkl")

i = 0
j = 0
acceptedRate_pre = 0
while i < 1:
    if not os.path.exists(pickle_file):
        import anndata as ad
        import pandas as pd
        import numpy as np
        import scipy as sp
        import seaborn as sns
        import copy
        import re
        import glob
        import matplotlib.pyplot as plt
        import warnings
        from datetime import datetime
        from igraph import *
        from sccoda.util import comp_ana as mod
        from sccoda.util import cell_composition_data as dat
        from sccoda.util import data_visualization as viz
        import arviz as az
        abundances = ad.read_h5ad(sys.argv[1])#
        abundances = abundances[abundances.obs.region != "MTG",] 
        print("Reference cell type: " + ct,flush = True)
        formula = "sex + "
        
        random_sleep_time = random.uniform(1,5)
        
        time.sleep(random_sleep_time)
        models = mod.CompositionalAnalysis(abundances, formula=formula + "C(" + k + ", Treatment('" + l + "'))", reference_cell_type=ct)
        results = models.sample_hmc()
        j = j + 1 
        print("now run" + str(j) + " times!", flush = True)
        accepted = results.sample_stats["is_accepted"].to_numpy()
        acceptedRate = accepted.sum() / accepted.shape[1]
        
        if acceptedRate_pre < acceptedRate:
            acceptedRate_pre = acceptedRate
            results_best = results
        print(str(acceptedRate),flush = True)
    
        if acceptedRate < 0.6:
            print("unacceptable!",flush = True)

        else:
            if not os.path.exists(pickle_file):
                results_best.save(pickle_file)
            i = i + 1
            print("passed!",flush = True)
    else:
            print("file exsit, passed!",flush = True)
            i = i + 1
  


print("Finish at: %s" % time.ctime(),flush = True)



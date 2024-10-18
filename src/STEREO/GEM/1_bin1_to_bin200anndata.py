# find . -name "*gem.gz" | xargs -P 5 -I {} ~/miniconda3/envs/st/bin/python3 ~/cortex/STEREO/GEM/1_bin1_to_bin200anndata.py {}
import warnings
warnings.filterwarnings('ignore')

import stereo as st

import os
import sys
#import glob
import pandas as pd
import datetime
os.chdir("/home/luomeng/cortex/STEREO/GEM")

file = sys.argv[1]

print("work with {} at ".format(file) + str(datetime.datetime.now()), flush = True)
data = st.io.read_gem(file_path=file, bin_type="bins", bin_size=200)
data.tl.raw_checkpoint()
st.io.stereo_to_anndata(data,flavor='seurat',output=file.replace("gem.gz","h5ad")) 
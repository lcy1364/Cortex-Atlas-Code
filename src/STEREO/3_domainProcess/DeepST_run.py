#  ######  run DeepST
# cd ~/DATA/luomeng/STEREO/DeepST
# ls data | grep 200.h5ad | xargs -I {} -P 5 ~/miniconda3/envs/deepst/bin/python DeepST.py {}
import os 
import sys
import datetime
print(f"{sys.argv[1]} start at {datetime.datetime.now()}", flush=True)
fileName = '/media/desk16/lwx/Epilepsy/upload/A01895D1/h5ad/A01895D1tissue.h5ad'
fileName = sys.argv[1]

if os.path.exists(fileName.replace("data2","result")):
    print("result exist, skip")
else:
    os.environ['OPENBLAS_NUM_THREADS'] = '5'
    import matplotlib.pyplot as plt
    from pathlib import Path
    import scanpy as sc
    import torch
    import pickle
    os.chdir("/home/luomeng/DATA/luomeng/STEREO/DeepST/")

    sys.path.append("/home/luomeng/DATA/luomeng/STEREO/DeepST/DeepST-main/deepst/")

    from DeepST import run


    #import glob
    #binFiles = glob.glob("data/*100.h5ad")
    #binFiles.__len__()

    print(f"{fileName} start at {datetime.datetime.now()}", flush=True)
    n_domains = 9 ###### the number of spatial domains.
    deepen = run(save_path = fileName.replace("data","result"),
             platform = "stereoseq",
             pca_n_comps = 200,
             pre_epochs = 500, #### According to your own hardware, choose the number of training
             epochs = 500, #### According to your own hardware, choose the number of training
             Conv_type="GCNConv", #### you can choose GNN types. 
             linear_encoder_hidden = [64, 16],
             conv_hidden=[64, 16])
    # adata = deepen._get_adata(data_path, data_name)


    adata = sc.read(fileName)
    adata.X = adata.X.toarray()
    #    data.tl.raw_checkpoint()

    adata = deepen._get_augment(adata, adjacent_weight = 0.3, neighbour_k = 4, weights="weights_matrix_nomd")
    graph_dict = deepen._get_graph(adata.obsm["spatial"], distType="KDTree",  k=20)
    adata = deepen._fit(adata, graph_dict, pretrain = True, dim_reduction = True)
    #adata = deepen._get_cluster_data(adata, n_domains = n_domains, priori=True)
    # adata = deepen._get_cluster_data(adata, n_domains = n_domains, priori=True) ###### without using prior knowledge, setting priori = False.

    ###### resolution=0.5

    sc.pp.neighbors(adata, use_rep='DeepST_embed')
    sc.tl.leiden(adata, key_added="DeepST_domain_5", resolution=0.5)
    print(f"done 5 job {fileName} at {datetime.datetime.now()}", flush = True)
    ###### resolution = 1
    sc.pp.neighbors(adata, use_rep='DeepST_embed')
    sc.tl.leiden(adata, key_added="DeepST_domain_10", resolution=1)
    print(f"done 10 job {fileName} at {datetime.datetime.now()}", flush = True)
    ######## save data
    adata.write(fileName.replace("data2","result"), compression="gzip")

print(f"done {fileName} at {datetime.datetime.now()}", flush = True)


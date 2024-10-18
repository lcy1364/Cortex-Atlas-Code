import numpy as np
import pandas as pd

from numba import njit
from scipy.sparse import csr_matrix
from anndata import AnnData

from typing import Optional
from os import listdir
from pathlib import Path


@njit
def bin_indices(coords: np.ndarray, coord_min: np.ndarray = np.array([0,0]), binsize: int = 50) -> int:
    """
    Take a DNB coordinate, the mimimum coordinate and the binsize, 
    
    calculate the index of bins for the current coordinate. 

    Parameters
    ----------
    coord : `np.ndarray`
        Current x or y coordinate.
    coord_min : Union[list, np.ndarray]
        Minimal value for the current x or y coordinate on the entire tissue slide measured by the spatial
        transcriptomics.
    binsize : binsize
        Size of the bins to aggregate data.

    Returns
    -------
    num: `int`
        The bin index for the current coordinate.
    """
    num = np.floor((coords - coord_min) / binsize)
    return num.astype(np.uint32)


@njit
def centroids(bin_indices: np.ndarray, coord_min: float = 0, binsize: int = 50) -> float:
    """
    Take a bin index, the mimimum coordinate and the binsize, 
    
    calculate the centroid of the current bin.

    Parameters
    ----------
    bin_indices : `np.ndarray`
        The bin index for the current coordinate.
    coord_min : `float`
        Minimal value for the current x or y coordinate on the entire tissue slide measured by the spatial
        transcriptomics.
    binsize : `int`
        Size of the bins to aggregate data.

    Returns
    -------
    `int`
        The bin index for the current coordinate.
    """
    coord_centroids = coord_min + bin_indices * binsize + binsize / 2
    return coord_centroids


def get_points_convex(df: pd.DataFrame, type='Multipoint'):
    """get cell bin points property

    Parameters
    ----------
    df : pd.DataFrame
        grouped data frame

    Returns
    -------
    pd.DataFrame
        set of property
    """
    from shapely.geometry import MultiPoint
    points = df[['x', 'y']].values
    hull = MultiPoint(points).convex_hull
    x_centroid, y_centroid = (np.floor([hull.centroid.x, 
                                        hull.centroid.y])
                                .astype(np.int32))
    return pd.Series([str(hull), hull.area, x_centroid, y_centroid],
                     index=['contour', 'area', 'center_x', 'center_y'])
    

def read_gem(
    gem_file: str,
    bin_size: int = 50,
    sample: Optional[str] = None,
    label_column: Optional[str] = None,
    tag_column: Optional[str] = None,
    n_workers: Optional[int]  = None, 
    verbose: bool = False,
    sep: str ='\t'):  
    """
    read gem file and return anndata object

    Parameters
    ----------
    gem_file : `str`
        gem file
    bin_size : `int`, optional
        bin size, by default 50
    sample : `str`, optional
        sample name, by default None
    label_column : `str`, optional
        label name(eg. cell name), by default None
    n_workers : `int`, optional
        num of workers(only take effect when label column is not None), by default 8

    Returns
    -------
    Anndata
        anndata object
        
    Examples
    --------
    - for bin_type data 
    
    >>> adata = read_gem('data/GF_spatial/FP200000495BL_A4.bin1.Lasso.gem.gz', bin_size=50)
    
    - for cell bin data
    
    >>> adata = read_gem('data/GF_spatial/FP200000495BL_A4.bin1.Lasso.gem.gz', 
    ...                   label_column='label', n_workers=12) # use n_workers to speed up
    """

    # here modified according to spateo
    dtype = {
        "geneID": "category",  # geneID
        "x": np.uint32,  # x
        "y": np.uint32,  # y
        # Multiple different names for MID counts
        "MIDCounts" : np.uint16,
        "MIDCount"  : np.uint16,
        "UMICount"  : np.uint16,
        "UMICounts" : np.uint16,
        "EXONIC"    : np.uint16,  # spliced
        "INTRONIC"  : np.uint16,  # unspliced,
        'label'     : np.uint32
    }
    rename = {
        "MIDCounts" : "UMICounts",
        "MIDCount"  : "UMICounts",
        "UMICount"  : "UMICounts",
        "EXONIC"    : "spliced",
        "INTRONIC"  : "unspliced",
    }
    
    if label_column:
        dtype.update({label_column: np.str0})
        rename.update({label_column: "label"})
        
    if verbose:
        print('-> read gem file...')    
    df = pd.read_csv(
        gem_file, compression='infer',
        header=0, sep=sep, quotechar='"', 
        dtype=dtype, comment="#").rename(columns=rename)
       
    df.drop_duplicates(subset=['geneID', "x", "y"], inplace=True)
    props = None
      
    if verbose:
        print('get bin index...') # TODO it may has other tags like batch
    if label_column is not None:   
        bin_size = 'cellbin'
        # df = df[df['label'] > 0]   
        if verbose:    
            print('\tcaclulate cell bin coordinate...')
            
        if n_workers is not None:
            try:
                from pandarallel import pandarallel
            except ImportError:
                raise ImportError(
                    'Please install the pandarallel, `pip3 install pandarallel')
                
            pandarallel.initialize(progress_bar=True, nb_workers=n_workers)    
            props = pd.DataFrame(df
                                 .groupby('label')
                                 .parallel_apply(get_points_convex))
        else:
            props = (df
                     .groupby('label')
                     .apply(get_points_convex))  
    
    else:    
        if bin_size > 1:
            df[['x', 'y']] = bin_indices(
                df[['x', 'y']].values, np.array([0, 0]), bin_size)
            
            df['label'] = (
                df['x'].astype(np.str0) + '-' + df['y'].astype(np.str0)
                )
            
        else:
            raise ValueError('bin_size must be greater than 1')

    cell_lst = sorted(df['label'].unique())
    gene_lst = sorted(df['geneID'].unique())
    
    cell_dict = dict(zip(cell_lst, range(len(cell_lst))))
    gene_dict = dict(zip(gene_lst, range(len(gene_lst))))

    rows = df['label'].map(cell_dict).astype(np.int32).values
    cols = df['geneID'].map(gene_dict).astype(np.int32).values

    # 存成稀疏矩阵
    if verbose:
        print('save sparse matrix...')
    shape = (len(cell_lst), len(gene_lst))
    exp_mtx = csr_matrix(
        (df['UMICounts'].values, (rows, cols)), shape=shape
        )

    layers = {}
    layers['counts'] = exp_mtx.copy()
    
    if "spliced" in df.columns:
        layers['spliced'] = csr_matrix(
            (df["spliced"].values, (rows, rows)), shape=shape)
        
    if "unspliced" in df.columns:
        layers['unspliced'] = csr_matrix(
            (df["unspliced"].values, (rows, rows)), shape=shape)
        
    obs = pd.DataFrame(index=cell_lst)
    var = pd.DataFrame(index=gene_lst)
    
    if verbose:
        print('initialize AnnData object...', end='\n\n')
    adata = AnnData(X=exp_mtx, obs=obs, var=var, layers=layers)
    #! If we centroid the points per bins here 
    #! the distance between two bins is bin_size not 1
    if props is not None:
        adata.obs = props.loc[adata.obs_names, :]
        pos = adata.obs[['center_x', 'center_y']].values
        
    else:
        pos = np.array([x.split('-') for x in cell_lst], dtype=np.uint32)
        
    adata.obsm['spatial'] = pos
    adata.obs[['x', 'y']] = pos
    adata.uns["sample_info"] = dict()
    adata.uns["sample_info"]['bin_size'] = bin_size if label_column is None else 'cellbin'
    adata.uns["sample_info"]['sample_name'] = sample if sample else None

    print(adata) # rich print
    return adata


def read_C4_mtx(file_path: str,
                prefix: Optional[str] = None,
                dtype = np.uint32, 
                transpose: bool=True) -> AnnData:
    """
    Read C4 data.

    Parameters
    ----------
    file path
        The filepath.
    dtype
        Numpy data type.
        
    Returns
    -------
    Anndata
        anndata object
        
    Examples
    --------
    >>> adata = read_C4_mtx('path/with/barcode_mtx_feature')
    """
    
    path = Path(file_path)
    prefix = '' if prefix is None else prefix.lower()
    # could be rewritten accounting for dtype to be more performant
    files = listdir(path)
    mtx_file =  [i for i in files if f'{prefix}matrix' in i.lower()]
    barcode_file =  [i for i in files if f'{prefix}barcodes' in i.lower()]
    gene_file = [i for i in files 
                 if (f'{prefix}features' in i.lower() 
                     or f'{prefix}genes' in i.lower())]
    if any([len(f) > 1 for f in [mtx_file, barcode_file, gene_file]]):
        raise ValueError('Multiple files found')
    
    # for .mtx format
    from scipy.io import mmread
    X = mmread(path.joinpath(mtx_file[0])).astype(dtype)
    
    # initialize anndata
    X = csr_matrix(X.T) if transpose else csr_matrix(X)
    var = pd.read_csv(path.joinpath(gene_file[0]), header=None, names=['geneID'], index_col=0)
    obs = pd.read_csv(path.joinpath(barcode_file[0]), header=None, names=['cell'], index_col=0)
    
    adata = AnnData(X, obs=obs, var=var, dtype=dtype)
    adata.layers['counts'] = adata.X.copy()
    return adata



def output_DGE_table(
    adata: AnnData, 
    key:str = 'clusters', 
    top_n: int = 20, 
    save_path: Optional[str] = None
    ):
    # TODO: scanpy do not caculate pts with default param, catch this error

    columns_1 = ['names', 'scores', 'pvals',
                 'pvals_adj', 'logfoldchanges']
    columns_2 = ['pts', 'pts_rest']

    with pd.ExcelWriter(r'DGE_table.xlsx') as writer:
        for i in np.unique(adata.obs[key]):
            data_1 = [[a[i] for a in adata.uns['rank_genes_groups'][column][:top_n]]
                      for column in columns_1]
            gene_names = data_1[0]
            data_2 = [adata.uns['rank_genes_groups'][column].loc[gene_names, i].values
                      for column in columns_2]

            columns = columns_1 + columns_2
            data = data_1 + data_2

            df = pd.DataFrame(dict(zip(columns, data)))
            df.to_excel(writer, sheet_name=i, index=False)


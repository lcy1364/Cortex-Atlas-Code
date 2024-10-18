import read_io as gem
import sys
import os
gemfile = sys.argv[1]

#gemfile = "./B02111D1.gem"
import os
os.environ['JOBLIB_TEMP_FOLDER'] = '/tmp' 
adata = gem.read_gem(gemfile,label_column="CellID")
ha5ad_file=gemfile.replace(".gem",".cellbin.h5ad")
adata.write_h5ad(ha5ad_file)
import loompy as lp
import numpy as np
import scanpy as sc
import csv
import pandas as pd
x = sc.read_csv("./data1.csv").T
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create("./data1.loom",x.X.transpose(),row_attrs,col_attrs)
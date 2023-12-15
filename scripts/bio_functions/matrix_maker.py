import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import muon
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import csv

adataTH1 = muon.read_10x_h5("./data/updated/raw_data/TH1/filtered_feature_bc_matrix.h5")
adataTH1.var_names_make_unique()

sparse_array = adataTH1.mod['rna'].X.toarray()
mtx_file = './data/TH1_test/matrix.mtx'

print("Starting...")

with open(mtx_file, 'w') as file:
    file.write("%%MatrixMarket matrix coordinate real general\n")

    # Write dimensions
    rows, cols = len(sparse_array), len(sparse_array[0])
    file.write(f"{rows} {cols} {np.count_nonzero(sparse_array)}\n")

    for i, row in enumerate(sparse_array):
        print(f"Reached row {i}")
        for j, value in enumerate(row):
            if value != 0:
                file.write(f"{i + 1} {j + 1} {value}\n")

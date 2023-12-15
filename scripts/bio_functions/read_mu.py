import mudata as md
# import muon as mu

# mdata = MuData({'rna': adata_rna, 'atac': adata_atac})

data_file_path = "../data/Project_Brand_Perkins_copy.h5mu"

mudata = md.read(data_file_path)
print(mudata)
# available_datasets = mudata.list_datasets()

# Print the names of the datasets
# for dataset_name in available_datasets:
#   print(dataset_name)

# Close the MuData file
# mudata.close()

# -----

# mu.read_10x_h5(data_file_path)

# mudata = mu.MuData(data_file_path, mode='r')

# MuData object with n_obs × n_vars = 10000 × 80000 
# 2 modalities
#   rna:	10000 x 30000
#     var:	'gene_ids', 'feature_types', 'genome', 'interval'
#   atac:	10000 x 50000
#     var:	'gene_ids', 'feature_types', 'genome', 'interval'
#     uns:	'atac', 'files'

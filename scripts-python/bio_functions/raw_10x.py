# system imports
import os
# genetic imports
import muon
import mudata as md
import anndata as ad
# ML imports
import numpy as np
import pandas as pd
from scipy.sparse import vstack
# general imports
import gzip
from pathlib import Path
from typing import List, Optional


def generate_zip(path: str):
    """
    Use gzip to generate a compressed version of the matrix and associated files
    :return:
    """

    print("generating zips")

    for filename in os.listdir(path):
        file_path = os.path.join(path, filename)

        # Check if the path is a file (not a subdirectory)
        if os.path.isfile(file_path):
            print(file_path)
            with open(file_path, "rb") as f_in:
                with gzip.open(f"{file_path}.gz", "wb") as f_out:
                    f_out.writelines(f_in)


def initialize_output_folder(path: str):
    """
    Remove all files from output folder
    :param path:
    :return:
    """

    print("deleting all files from output folder")

    if not os.path.exists(path):
        os.makedirs(path)

    for filename in os.listdir(path):
        file_path = os.path.join(path, filename)

        if os.path.isfile(file_path):
            os.remove(file_path)
            print(f"Removed: {file_path}")


def cleanup_folder(path: str):
    """
    Remove all non-gz files from a folder
    :param path:
    :return:
    """

    print("cleaning up folder")

    for filename in os.listdir(path):
        file_path = os.path.join(path, filename)

        # Check if the path is a file and does not end with ".gz"
        if os.path.isfile(file_path) and not filename.endswith(".gz"):
            os.remove(file_path)
            print(f"Removed: {file_path}")


def write_tsv(obj: pd.DataFrame, name: str, index: bool = True) -> None:
    """
    Utility function to write to tsv files
    :param obj:
    :param name:
    :param index:
    :return:
    """

    print(f"writing {name}")

    directory = os.path.dirname(name)
    if not os.path.exists(directory):
        os.makedirs(directory)

    obj.to_csv(name, sep="\t", header=False, index=index)


def write_matrix(sparse_array: np.ndarray, name: str) -> None:
    """
    Utility function to write to matrix files

    :param sparse_array:
    :param name:
    :return:
    """

    print(f"writing {name}")

    sparse_nonzero = np.nonzero(sparse_array)
    nonzero_count = np.count_nonzero(sparse_array)

    directory = os.path.dirname(name)
    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(name, "w") as file:
        file.write("%%MatrixMarket matrix coordinate real general\n")

        # Write dimensions
        rows, cols = len(sparse_array), len(sparse_array[0])
        file.write(f"{cols} {rows} {np.count_nonzero(sparse_array)}\n")

        for i in range(nonzero_count):
            col = sparse_nonzero[1][i]
            row = sparse_nonzero[0][i]

            file.write(f"{col + 1} {row + 1} {sparse_array[row][col]}\n")


def raw_10x(
    subsets: List[str],
    output_folder_path: str,
    input_folder_path: str,
    cell_sample_size: Optional[int] = None,
    gene_sample_size: Optional[int] = None,
    balanced_dist: bool = False
) -> None:
    """
    Main function
    :return:
    """

    # NOTE - Equivalent:
    # adataTH1.var[adataTH1.var["feature_types"] == "Gene Expression"]  # Antibody Capture
    # adataTH1.var[adataTH1.var["genome"] != ""]
    # adataTH1.mod["rna"].var

    if len(subsets) == 0:
        print("Must have at least one subset specified!")
        return

    print("Starting extraction process...")

    adata_obs_l = []
    adata_rna_l = []
    adata_protein_l = []

    rna_var = None
    prot_var = None

    seed_value = 42  # Seed for reproducibility

    for subset in subsets:
        print(f"Importing {subset}...")

        adata_sub = muon.read_10x_h5(
            Path(f"{input_folder_path}/raw_data/{subset}/filtered_feature_bc_matrix.h5")
        )
        adata_sub.var_names_make_unique()
        adata_sub.obs.index = [name + subset for name in adata_sub.obs_names]

        # balanced selection
        if balanced_dist:
            if cell_sample_size:
                selected_rows = adata_sub.obs.sample(n=cell_sample_size, random_state=seed_value)
                adata_sub = adata_sub[selected_rows.index]

            if gene_sample_size:
                selected_cols = adata_sub.var.sample(n=gene_sample_size, random_state=seed_value)
                adata_sub = adata_sub[:, selected_cols.index]

        adata_obs_l.append(adata_sub.obs)
        adata_rna_l.append(adata_sub.mod["rna"].X)
        adata_protein_l.append(adata_sub.mod["prot"].X)

        rna_var = adata_sub.mod["rna"].var
        prot_var = adata_sub.mod["prot"].var

    obs_merged = pd.concat(adata_obs_l)
    rna_merged = vstack(adata_rna_l)
    prot_merged = vstack(adata_protein_l)

    adata_merged = md.MuData({
        "rna": ad.AnnData(X=rna_merged, var=rna_var),
        "prot": ad.AnnData(X=prot_merged, var=prot_var)}
    )
    adata_merged.obs = obs_merged

    if not balanced_dist:  # balanced selection already handled above
        if cell_sample_size:
            selected_rows = adata_merged.obs.sample(n=cell_sample_size, random_state=seed_value)
            adata_merged = adata_merged[selected_rows.index]

        if gene_sample_size:
            selected_cols = adata_merged.var.sample(n=gene_sample_size, random_state=seed_value)
            adata_merged = adata_merged[:, selected_cols.index]

    sparse_array = adata_merged.mod["rna"].X.toarray()
    # rows, cols = len(sparse_array), len(sparse_array[0])

    # Prepare vars dataframe
    sub_vars = adata_merged.mod["rna"].var.reset_index()[["gene_ids", "index", "feature_types"]]

    # Store in a local folder
    write_matrix(sparse_array, f"{output_folder_path}/matrix.mtx")  # Count matrix
    write_tsv(adata_merged.obs, f"{output_folder_path}/barcodes.tsv")  # Cell barcodes
    write_tsv(sub_vars, f"{output_folder_path}/features.tsv", index=False)  # Gene labels

    generate_zip(output_folder_path)
    cleanup_folder(output_folder_path)

import os
import time
import muon
from pathlib import Path
import gzip
from typing import List, Optional
import numpy as np
from pandas import DataFrame


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
            with open(file_path, 'rb') as f_in:
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

        # Check if the path is a file and does not end with '.gz'
        if os.path.isfile(file_path) and not filename.endswith('.gz'):
            os.remove(file_path)
            print(f"Removed: {file_path}")


def write_tsv(obj: DataFrame, name: str, index: bool = True) -> None:
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

    with open(name, 'w') as file:
        file.write("%%MatrixMarket matrix coordinate real general\n")

        # Write dimensions
        rows, cols = len(sparse_array), len(sparse_array[0])
        file.write(f"{cols} {rows} {np.count_nonzero(sparse_array)}\n")

        for i in range(nonzero_count):
            col = sparse_nonzero[1][i]
            row = sparse_nonzero[0][i]

            file.write(f"{col + 1} {row + 1} {sparse_array[row][col]}\n")


def main(
    subsets: List[str],
    output_folder_path: str,
    input_folder_path: str,
    cell_sample_size: Optional[int] = None,
    gene_sample_size: Optional[int] = None
) -> None:
    """
    Main function
    :return:
    """

    # NOTE - Equivalent:
    # adataTH1.var[adataTH1.var['feature_types'] == "Gene Expression"]  # Antibody Capture
    # adataTH1.var[adataTH1.var['genome'] != ""]
    # adataTH1.mod['rna'].var

    print("Starting extraction process...")

    adata_sub = None

    # TODO: add merging code
    for subset in subsets:
        adata_sub = muon.read_10x_h5(
            Path(f"{input_folder_path}/raw_data/{subset}/filtered_feature_bc_matrix.h5")
        )
        adata_sub.var_names_make_unique()

    if cell_sample_size:
        seed_value = 42  # Seed for reproducibility
        selected_rows = adata_sub.obs.sample(n=cell_sample_size, random_state=seed_value)
        adata_sub = adata_sub[selected_rows.index]

    if gene_sample_size:
        seed_value = 42  # Seed for reproducibility
        selected_cols = adata_sub.var.sample(n=gene_sample_size, random_state=seed_value)
        adata_sub = adata_sub[:, selected_cols.index]

    sparse_array = adata_sub.mod['rna'].X.toarray()
    # rows, cols = len(sparse_array), len(sparse_array[0])

    # Prepare vars dataframe
    sub_vars = adata_sub.mod['rna'].var.reset_index()[["gene_ids", "index", "feature_types"]]

    write_matrix(sparse_array, f"{output_folder_path}/matrix.mtx")  # Count matrix
    write_tsv(adata_sub.obs, f"{output_folder_path}/barcodes.tsv")  # Cell barcodes
    write_tsv(sub_vars, f"{output_folder_path}/features.tsv", index=False)  # Gene labels

    generate_zip(output_folder_path)
    cleanup_folder(output_folder_path)


# store in a local folder

if __name__ == "__main__":
    start_time = time.process_time()

    output_folder = "./output_10x"
    input_folder = "D:/genomics/data/updated"

    initialize_output_folder(output_folder)

    try:
        main(
            subsets=["TH1"],
            input_folder_path=input_folder,
            output_folder_path=output_folder,
            gene_sample_size=1000
        )
    except Exception as e:
        print(e)

    total_time = time.process_time() - start_time
    print(f"Total CPU time: {total_time:.3f} seconds")

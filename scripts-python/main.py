import time
from bio_functions.raw_10x import raw_10x, initialize_output_folder

if __name__ == "__main__":
    start_time = time.process_time()

    output_folder = "./output_10x"
    input_folder = "E:/genomics/data/updated"

    initialize_output_folder(output_folder)

    try:
        raw_10x(
            subsets=["BM4", "PB2", "TH1", "TH2"],
            input_folder_path=input_folder,
            output_folder_path=output_folder,
            # gene_sample_size=10000,
            cell_sample_size=500,
            balanced_dist=True
        )

    except Exception as e:
        print(e)

    total_time = time.process_time() - start_time
    print(f"Total CPU time: {total_time:.3f} seconds")

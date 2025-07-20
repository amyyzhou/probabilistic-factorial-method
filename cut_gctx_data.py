# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python [conda env:base] *
#     language: python
#     name: conda-base-py
# ---

# %%
import pandas as pd
import numpy as np
import os # For checking file existence
import h5py # For reading HDF5 files directly

# --- Configuration ---
# IMPORTANT: Ensure these file paths are correct for your local environment.
file_paths = {
    'sig_info': 'GSE70138_Broad_LINCS_sig_info_2017-03-06.txt',
    'pert_info': 'GSE70138_Broad_LINCS_pert_info.txt',
    'gene_info': 'GSE70138_Broad_LINCS_gene_info_2017-03-06.txt',
    'gctx_data': 'GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx'
}

# --- Parameters for cutting down the GCTX file ---
# Maximum number of A375 active drug samples to include in the subset.
# Set to None to include all available A375 active drug samples.
# You can adjust this based on your system's memory.
MAX_SAMPLES_TO_EXTRACT = 2000

# Number of genes to randomly select and extract.
# Set to None to extract all genes that are in the filtered samples.
# You can adjust this based on your system's memory and analysis needs.
NUM_GENES_TO_EXTRACT = 1000 

# Output file path for the cut-down data.
# We'll save it as an HDF5 file, which is efficient for numerical data.
OUTPUT_FILE_NAME = 'GSE70138_A375_subset_expression.h5' 

# Set random seed for reproducibility of random sampling
np.random.seed(42)

print("--- Starting GCTX File Reduction Script (using h5py) ---")

# --- 1. Load All Necessary Metadata Files ---
try:
    sig_info_df = pd.read_csv(file_paths['sig_info'], sep='\t')
    pert_info_df = pd.read_csv(file_paths['pert_info'], sep='\t')
    gene_info_df = pd.read_csv(file_paths['gene_info'], sep='\t')
    print("All metadata files loaded successfully.")
except FileNotFoundError as e:
    print(f"Error loading a metadata file: {e}")
    print("Please ensure all .txt files are in the correct directory.")
    exit("Exiting due to missing metadata files.")


# --- 2. Identify Target Samples (A375 active drug treatments) ---
print("\n--- Identifying Target Samples (A375 active drug treatments) ---")

# Filter sig_info_df to get only A375 samples
a375_signatures = sig_info_df[sig_info_df['cell_id'] == 'A375'].copy()
print(f"Initial A375 signatures found: {len(a375_signatures)}")

# Filter these A375 signatures to include only 'trt_cp' (compound treatments)
# and explicitly exclude 'DMSO' (control vehicle).
active_drug_signatures = a375_signatures[
    (a375_signatures['pert_type'] == 'trt_cp') &
    (a375_signatures['pert_iname'] != 'DMSO')
].copy()

print(f"A375 signatures with active compound treatments (excluding DMSO): {len(active_drug_signatures)}")

# Get the list of sig_ids for these relevant samples
# IMPORTANT: Convert to string and strip whitespace from sig_ids
target_sig_ids = [str(sid).strip() for sid in active_drug_signatures['sig_id'].tolist()]

# Apply sample limitation if configured
if MAX_SAMPLES_TO_EXTRACT is not None and len(target_sig_ids) > MAX_SAMPLES_TO_EXTRACT:
    print(f"Limiting target samples to {MAX_SAMPLES_TO_EXTRACT} randomly selected ones.")
    target_sig_ids = np.random.choice(target_sig_ids, MAX_SAMPLES_TO_EXTRACT, replace=False).tolist()
else:
    print(f"Using all {len(target_sig_ids)} identified A375 active drug samples.")


# --- 3. Identify Target Genes ---
print("\n--- Identifying Target Genes ---")

# For simplicity, let's target 'landmark' genes (pr_is_lm = 1) as they are well-measured.
# You could also choose all genes, or a specific list of genes if desired.
# IMPORTANT: Convert to string and strip whitespace from gene_ids
landmark_gene_ids = [str(gid).strip() for gid in gene_info_df[gene_info_df['pr_is_lm'] == 1]['pr_gene_id'].tolist()]
print(f"Total landmark genes identified: {len(landmark_gene_ids)}")

# Apply gene limitation if configured
genes_to_extract_ids = []
if NUM_GENES_TO_EXTRACT is not None and len(landmark_gene_ids) > NUM_GENES_TO_EXTRACT:
    print(f"Randomly selecting {NUM_GENES_TO_EXTRACT} genes to extract.")
    genes_to_extract_ids = np.random.choice(landmark_gene_ids, NUM_GENES_TO_EXTRACT, replace=False).tolist()
else:
    print(f"Extracting all {len(landmark_gene_ids)} landmark genes.")
    genes_to_extract_ids = landmark_gene_ids

# --- 4. Load and Filter GCTX Data using h5py ---
print("\n--- Loading and Filtering GCTX Data using h5py ---")

expression_subset_df = None
try:
    # Open the GCTX file in read mode
    with h5py.File(file_paths['gctx_data'], 'r') as f:
        # --- UPDATED PATHS BASED ON INSPECTION ---
        # Load the ENTIRE data matrix into a NumPy array first to allow 2D fancy indexing
        print("Loading full data matrix into memory... (This might take a moment)")
        full_data_matrix_np = f['/0/DATA/0/matrix'][()] # [()] loads the entire dataset into a NumPy array
        
        # IMPORTANT: Convert to string and strip whitespace from IDs read from GCTX
        all_gctx_sample_ids = [cid.decode('utf-8').strip() for cid in f['/0/META/COL/id'][()]] # These are the sample IDs (sig_ids)
        all_gctx_gene_ids = [rid.decode('utf-8').strip() for rid in f['/0/META/ROW/id'][()]] # These are the gene IDs (pr_gene_ids)

        # --- DEBUGGING: Print example IDs ---
        print("\n--- Debugging ID Mismatch ---")
        print("Example target_sig_ids (from metadata):", target_sig_ids[:5])
        print("Example all_gctx_sample_ids (from GCTX):", all_gctx_sample_ids[:5])
        print("Example genes_to_extract_ids (from metadata):", genes_to_extract_ids[:5])
        print("Example all_gctx_gene_ids (from GCTX):", all_gctx_gene_ids[:5])
        print("-----------------------------\n")

        # Find the indices for the target samples and genes that actually exist in the GCTX file
        # Use a set for faster lookup
        gctx_sample_id_set = set(all_gctx_sample_ids)
        gctx_gene_id_set = set(all_gctx_gene_ids)

        # Filter target_sig_ids and genes_to_extract_ids to only include those present in the GCTX file
        filtered_target_sig_ids = [sid for sid in target_sig_ids if sid in gctx_sample_id_set]
        filtered_genes_to_extract_ids = [gid for gid in genes_to_extract_ids if gid in gctx_gene_id_set]

        if not filtered_target_sig_ids or not filtered_genes_to_extract_ids:
            # Added more specific error message to help diagnose
            missing_samples = len(target_sig_ids) - len(filtered_target_sig_ids)
            missing_genes = len(genes_to_extract_ids) - len(filtered_genes_to_extract_ids)
            error_msg = f"No matching genes or samples found in GCTX file after filtering. " \
                        f"({missing_samples} samples, {missing_genes} genes missing from GCTX IDs)."
            raise ValueError(error_msg)

        # Get the numerical indices for the filtered IDs
        # Important: The GCTX matrix is (samples, genes), so sample_indices are rows, gene_indices are columns
        # We get the indices in the original order first.
        sample_indices_original_order = [all_gctx_sample_ids.index(sig_id) for sig_id in filtered_target_sig_ids]
        gene_indices_original_order = [all_gctx_gene_ids.index(gene_id) for gene_id in filtered_genes_to_extract_ids]

        # --- IMPORTANT: Convert indices to NumPy arrays for direct fancy indexing ---
        # Use sorted indices with np.ix_ to ensure compatibility with HDF5/NumPy indexing.
        sorted_sample_indices_np = np.array(sorted(sample_indices_original_order))
        sorted_gene_indices_np = np.array(sorted(gene_indices_original_order))

        # --- DEBUGGING: Print index array info ---
        print(f"sample_indices_np shape: {sorted_sample_indices_np.shape}, dtype: {sorted_sample_indices_np.dtype}")
        print(f"gene_indices_np shape: {sorted_gene_indices_np.shape}, dtype: {sorted_gene_indices_np.dtype}")
        print(f"full_data_matrix_np shape: {full_data_matrix_np.shape}")
        print("-----------------------------\n")

        # Load the subset of the data matrix using np.ix_ on the full NumPy array.
        # This will extract a rectangular block.
        subset_data_transposed_raw = full_data_matrix_np[np.ix_(sorted_sample_indices_np, sorted_gene_indices_np)]

        # Create a temporary DataFrame with the data loaded using sorted indices.
        # The index (rows) will be the sample IDs corresponding to sorted sample indices.
        # The columns will be the gene IDs corresponding to sorted gene indices.
        temp_df_transposed = pd.DataFrame(
            subset_data_transposed_raw,
            index=[all_gctx_sample_ids[i] for i in sorted_sample_indices_np],
            columns=[all_gctx_gene_ids[i] for i in sorted_gene_indices_np]
        )

        # --- Reindex to get the desired original order for samples and genes ---
        # This step reorders the data to match the order of `filtered_target_sig_ids`
        # and `filtered_genes_to_extract_ids` as derived from your metadata.
        # Then, transpose to get (genes x samples)
        expression_subset_df = temp_df_transposed.reindex(
            index=filtered_target_sig_ids, # Reorder rows (samples)
            columns=filtered_genes_to_extract_ids # Reorder columns (genes)
        ).T # Final transpose to get genes as rows, samples as columns
    
    print(f"Successfully loaded and filtered subset using h5py: {expression_subset_df.shape} "
          f"(rows: {expression_subset_df.shape[0]} genes, columns: {expression_subset_df.shape[1]} samples)")

except FileNotFoundError:
    print(f"GCTX file not found at {file_paths['gctx_data']}.")
    print("Cannot proceed with gene expression analysis without the GCTX file.")
except KeyError as e:
    print(f"KeyError: Could not find expected HDF5 dataset in GCTX file: {e}")
    print("Double-check the HDF5 paths ('/0/DATA/0/matrix', '/0/META/COL/id', '/0/META/ROW/id').")
except ValueError as e:
    print(f"ValueError during data extraction: {e}")
except Exception as e:
    print(f"An unexpected error occurred during GCTX loading/filtering with h5py: {e}")
    print("Consider checking file integrity or memory availability. You might need to reduce NUM_GENES_TO_EXTRACT or MAX_SAMPLES_TO_EXTRACT further.")

if expression_subset_df is None:
    exit("Exiting because GCTX data subset could not be loaded using h5py.")


# --- 5. Save the Cut-Down Data ---
print(f"\n--- Saving Cut-Down Data to {OUTPUT_FILE_NAME} ---")

try:
    if OUTPUT_FILE_NAME.endswith('.csv'):
        expression_subset_df.to_csv(OUTPUT_FILE_NAME, index=True)
    elif OUTPUT_FILE_NAME.endswith('.h5'):
        # For HDF5, you might need to specify a key (e.g., 'expression_data')
        expression_subset_df.to_hdf(OUTPUT_FILE_NAME, key='expression_data', mode='w')
    else:
        print(f"Unsupported output file format: {OUTPUT_FILE_NAME}. Please use .csv or .h5.")
        exit("Exiting: Unsupported output format.")
    
    print(f"Cut-down data saved successfully to {OUTPUT_FILE_NAME}.")
    print(f"New file size: {os.path.getsize(OUTPUT_FILE_NAME) / (1024*1024):.2f} MB")

except Exception as e:
    print(f"Error saving the cut-down data: {e}")
    exit("Exiting due to save error.")

print("\n--- GCTX File Reduction Complete ---")
print(f"You can now use '{OUTPUT_FILE_NAME}' in your main analysis script.")
print("Remember to update the 'file_paths' dictionary in your main script to point to this new file.")

# %%

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
import h5py
import os

gctx_filepath = 'GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx'

print(f"--- Inspecting GCTX file: {gctx_filepath} ---")

if not os.path.exists(gctx_filepath):
    print(f"Error: File not found at {gctx_filepath}")
    print("Please ensure the GCTX file is in the same directory as this script.")
else:
    try:
        with h5py.File(gctx_filepath, 'r') as f:
            print("\nInternal structure of the GCTX file:")
            
            # Function to recursively print groups and datasets
            def print_hdf5_item(name, obj):
                if isinstance(obj, h5py.Group):
                    print(f"  Group: /{name}")
                elif isinstance(obj, h5py.Dataset):
                    print(f"    Dataset: /{name} (Shape: {obj.shape}, Dtype: {obj.dtype})")
            
            f.visititems(print_hdf5_item)

            print("\nCommon GCTX paths to look for:")
            print("  - Data Matrix: Usually '/data/matrix' or similar.")
            print("  - Row IDs (Genes): Usually '/meta/row_ids', '/row_ids', or similar.")
            print("  - Column IDs (Samples): Usually '/meta/col_ids', '/col_ids', or similar.")

    except Exception as e:
        print(f"An error occurred while inspecting the GCTX file: {e}")
        print("This might happen if the file is corrupted or not a valid HDF5 file.")

print("\n--- Inspection Complete ---")

# %%

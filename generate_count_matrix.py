#!/usr/bin/env/python

import pandas as pd
import os
from sys import argv

path = argv[1]
data_frames = []

for file in os.listdir(path):
    # Join the path and file together to create the full path to the file
    full_path = os.path.join(path, file)
    
    if file != ".ipynb_checkpoints":
        try:
            data = pd.read_csv(full_path, sep="\t", header=0)
        
            # Check if the DataFrame is empty
            if data.empty:
                print(f"Skipping empty file: {file}")
                continue
        
            sample_name = file.replace(".bracken", "")
            data["Sample"] = sample_name
            data = data.drop(["taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads", "added_reads", "fraction_total_reads" ], axis=1)
            data_frames.append(data)
        
        except pd.errors.EmptyDataError:
            print(f"Skipping empty file: {file}")
            continue

df = pd.concat(data_frames, ignore_index=True)

result = df.pivot_table(index="name", columns="Sample", values="new_est_reads", fill_value=0)
result = result.reset_index()

print("count_matrix saved to analysis/plots")


result.to_csv("analysis/plots/bracken_count_matrix.csv", index=False)
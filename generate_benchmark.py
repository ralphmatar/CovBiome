import pandas as pd
import os

path = "analysis/benchmarking/speed"
data_frames = []

for file in os.listdir(path):
    # Join the path and file together to create the full path to the file
    full_path = os.path.join(path, file)
    
    if file != ".ipynb_checkpoints":
            data = pd.read_csv(full_path, sep="\t", header=0)
            sample_name = file.replace(".tsv", "")
            data["Sample"] = sample_name
            data = data.drop(["max_rss", "max_vms", "max_uss", "mean_load", "io_in", "io_out", "max_pss", "h:m:s"], axis=1)
            data_frames.append(data)
        

df = pd.concat(data_frames)
df = df.set_index("Sample")


df.to_csv("analysis/benchmarking/speed/performance.csv", index=True)
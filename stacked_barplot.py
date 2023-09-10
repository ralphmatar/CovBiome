import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV file into a DataFrame
df = pd.read_csv("analysis/benchmarking/speed/performance.csv")

# Extract the sample number from the index
df['Name'] = df['Sample'].str.split('_').str[-1]

# Convert the 'Name' column to numeric type
df['Name'] = pd.to_numeric(df['Name'])

# Exclude specific prefixes
excluded_prefixes = ['map2hb_', 'umhyb_', 'umhum_', 'umcov_', 'map2hu_']
df_map2cov = df[~df['Sample'].str.startswith(tuple(excluded_prefixes))]

# Calculate the sum for each sample
df_map2cov = df_map2cov.groupby('Name').sum()
df_map2cov = df_map2cov.drop('Sample', axis=1)
df_map2cov.rename(columns={'s': 's_map2cov', 'cpu_time': 'cpu_map2cov'}, inplace=True)


# Exclude specific prefixes
excluded_prefixes_2 = ['umhyb_', 'umhum_', 'umcov_', 'map2hu_', 'map2cov_']
df_map2hb = df[~df['Sample'].str.startswith(tuple(excluded_prefixes_2))]

# Calculate the sum for each sample
df_map2hb = df_map2hb.groupby('Name').sum()
df_map2hb = df_map2hb.drop('Sample', axis=1)
df_map2hb.rename(columns={'s': 's_map2hb', 'cpu_time': 'cpu_map2hb'}, inplace=True)

# Exclude specific prefixes
excluded_prefixes_3 = ['map2hb_', 'umhyb_', 'umhum_', 'umcov_', 'map2cov_']
df_map2hu = df[~df['Sample'].str.startswith(tuple(excluded_prefixes_3))]

# Calculate the sum for each sample
df_map2hu = df_map2hu.groupby('Name').sum()
df_map2hu = df_map2hu.drop('Sample', axis=1)
df_map2hu.rename(columns={'s': 's_map2hu', 'cpu_time': 'cpu_map2hu'}, inplace=True)

# Exclude specific prefixes
excluded_prefixes_4 = ['map2hb_', 'umhum_', 'umcov_', 'map2cov_', 'map2hu']
df_umhyb = df[~df['Sample'].str.startswith(tuple(excluded_prefixes_4))]

# Calculate the sum for each sample
df_umhyb = df_umhyb.groupby('Name').sum()
df_umhyb = df_umhyb.drop('Sample', axis=1)
df_umhyb.rename(columns={'s': 's_umhyb', 'cpu_time': 'cpu_umhyb'}, inplace=True)

# Exclude specific prefixes
excluded_prefixes_5 = ['map2hb_', 'umhyb_', 'umcov_', 'map2cov_', 'map2hu']
df_umhum = df[~df['Sample'].str.startswith(tuple(excluded_prefixes_5))]

# Calculate the sum for each sample
df_umhum = df_umhum.groupby('Name').sum()
df_umhum = df_umhum.drop('Sample', axis=1)
df_umhum.rename(columns={'s': 's_umhum', 'cpu_time': 'cpu_umhum'}, inplace=True)

# Exclude specific prefixes
excluded_prefixes_6 = ['map2hb_', 'umhyb_', 'umhum_', 'map2cov_', 'map2hu']
df_umcov = df[~df['Sample'].str.startswith(tuple(excluded_prefixes_6))]

# Calculate the sum for each sample
df_umcov = df_umcov.groupby('Name').sum()
df_umcov = df_umcov.drop('Sample', axis=1)
df_umcov.rename(columns={'s': 's_umcov', 'cpu_time': 'cpu_umcov'}, inplace=True)

# adjust missing values
final_df = pd.concat([df_umcov, df_umhum, df_umhyb, df_map2hb, df_map2hu, df_map2cov], axis=1)
final_df.at[95, "s_umhyb"] = 8.4858
final_df.at[95, "cpu_umhyb"] = 2.22
final_df.at[599, "s_umhyb"] = 1.8005
final_df.at[599, "cpu_umhyb"] = 0.45

final_df.to_csv('analysis/benchmarking/speed/final_performance.csv')
final_df

# barplot
fig, ax = plt.subplots()
x = np.arange(len(final_df.index))
width = 0.35
gap = 0.05
ax.bar(x - width/2 - gap/2, final_df['s_umcov'], width, label='get unmapped to cov')
ax.bar(x - width/2 - gap/2, final_df['s_umhum'], width, bottom=final_df['s_umcov'], label='get unmapped to hg', color='yellow')
ax.bar(x - width/2 - gap/2, final_df['s_map2cov'], width, bottom=final_df['s_umcov'] + final_df['s_umhum'], label='map to cov', color='green')
ax.bar(x - width/2 - gap/2, final_df['s_map2hu'], width, bottom=final_df['s_umcov'] + final_df['s_umhum'] + final_df['s_map2cov'], label='map to hg', color='orange')
ax.bar(x + width/2 + gap/2, final_df['s_map2hb'], width, label='map to hybrid', color="0.3")
ax.bar(x + width/2 + gap/2, final_df['s_umhyb'], width, bottom=final_df['s_map2hb'], label='get unmapped to hybrid', color="0.6")
ax.set_xticks(x, final_df.index, size=8)
ax.set_xlabel('Sample')
ax.set_ylabel('Real time in seconds')
ax.set_title('Performance of concatenated and separate genomes')
ax.legend(prop={'size': 6})

fig.savefig("analysis/benchmarking/speed/performance_bp.pdf", dpi=400)
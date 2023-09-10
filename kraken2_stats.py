import os
import pandas as pd

def extract_percentages(directory):
    data = []
    for filename in os.listdir(directory):
        if filename.endswith('.kraken2'):
            sample_name = filename.split('.')[0]
            with open(os.path.join(directory, filename), 'r') as f:
                lines = f.readlines()
                unclassified_percentage = float(lines[0].split()[0])
                classified_percentage = 100 - unclassified_percentage
                data.append([sample_name, classified_percentage, unclassified_percentage])
    df = pd.DataFrame(data, columns=['Sample', 'Classified', 'Unclassified'])
    df.set_index('Sample', inplace=True)
    return df

df = extract_percentages('analysis/classification/')
print(df)
import os
import re
import pandas as pd

def extract_data_from_txt(filepath):
    with open(filepath, 'r') as file:
        content = file.read()

    number_of_reads_match = re.search(r'number of reads = ([\d,]+)', content)
    number_of_mapped_reads_match = re.search(r'number of mapped reads = ([\d,]+)', content)
    mapped_reads_percentage_match = re.search(r'number of mapped reads = [\d,]+ \((\d+\.\d+)%\)', content)
    mean_coverage_data_match = re.search(r'mean coverageData = ([\d.]+)X', content)

    number_of_reads = int(number_of_reads_match.group(1).replace(',', '')) if number_of_reads_match else None
    number_of_mapped_reads = int(number_of_mapped_reads_match.group(1).replace(',', '')) if number_of_mapped_reads_match else None
    mapped_reads_percentage = float(mapped_reads_percentage_match.group(1)) if mapped_reads_percentage_match else None
    mean_coverage_data = float(mean_coverage_data_match.group(1)) if mean_coverage_data_match else None

    return number_of_reads, number_of_mapped_reads, mapped_reads_percentage, mean_coverage_data

def main():
    root_dir = os.path.abspath('analysis/QC/qualimap')
    data = []

    for folder_name in os.listdir(root_dir):
        folder_path = os.path.join(root_dir, folder_name)
        if os.path.isdir(folder_path):
            txt_file_path = os.path.join(folder_path, 'genome_results.txt')
            if os.path.exists(txt_file_path):
                number_of_reads, number_of_mapped_reads, mapped_reads_percentage, mean_coverage_data = extract_data_from_txt(txt_file_path)
                data.append([folder_name, number_of_reads, number_of_mapped_reads, mapped_reads_percentage, mean_coverage_data])

    df = pd.DataFrame(data, columns=['Sample', 'Number of Reads', 'Number of Mapped Reads', 'Mapped Reads Percentage', 'Mean CoverageData'])
    
    
    p1 = df.head(8)
    p2 = df.tail(8)
    dfs = [p1, p2]
    
    df_thesis = pd.concat(dfs, ignore_index=True)
    df_thesis.to_csv('analysis/QC/qualimap/qualimap_stats.csv', index=False)
    df_thesis.to_latex('analysis/QC/qualimap/qualimap_stats.txt', index=False)
    df.to_csv('analysis/QC/qualimap/complete_stats.csv', index=False)

if __name__ == "__main__":
    main()
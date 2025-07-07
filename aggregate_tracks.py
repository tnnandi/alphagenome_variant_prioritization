# aggregate the tracks for each variant and save the summary to a csv file
# aggregation is done using the quantile score of the tracks
# the summary includes the number of extreme tracks, upregulated tracks, and downregulated tracks
# extreme tracks are tracks with a quantile score greater than 0.95 or less than -0.95
# upregulated tracks are tracks with a quantile score greater than 0.95
# downregulated tracks are tracks with a quantile score less than -0.95


import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from pdb import set_trace

# tracks_df = pd.read_csv("./predictions/chr_17_7571752_T_G_EUR_15820_scores.csv")
# tracks_df.columns: ['variant_id', 'scored_interval', 'gene_id', 'gene_name', 'gene_type',
    #    'gene_strand', 'junction_Start', 'junction_End', 'output_type',
    #    'variant_scorer', 'track_name', 'track_strand', 'Assay title',
    #    'ontology_curie', 'biosample_name', 'biosample_type',
    #    'transcription_factor', 'histone_mark', 'gtex_tissue', 'raw_score',
    #    'quantile_score']


input_dir = "/grand/GeomicVar/tarak/MVP_alphaGenome/predictions"
output_summary = "variant_quantile_summary_parallel.csv"

def process_file(file_path):
    try:
        df = pd.read_csv(file_path, usecols=["variant_id", "quantile_score"])
        if df.empty:
            return None

        variant_id = df["variant_id"].iloc[0]
        n_extreme = (df["quantile_score"].abs() > 0.95).sum() # number of tracks with a quantile score greater than 0.95 or less than -0.95 (extreme tracks)    
        n_up = (df["quantile_score"] > 0.95).sum() # number of tracks with a quantile score greater than 0.95
        n_down = (df["quantile_score"] < -0.95).sum() # number of tracks with a quantile score less than -0.95

        return {
            "filename": os.path.basename(file_path),
            "variant_id": variant_id,
            "n_extreme_tracks": n_extreme,
            "n_upregulated_tracks": n_up,
            "n_downregulated_tracks": n_down,
        }
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None

# Collect all *_scores.csv files
all_files = [
    os.path.join(input_dir, f)
    for f in os.listdir(input_dir)
    if f.endswith("_scores.csv")
]

# Use all available CPU threads
max_threads = os.cpu_count()

# Process files with progress bar
results = []
with ThreadPoolExecutor(max_workers=max_threads) as executor:
    futures = {executor.submit(process_file, f): f for f in all_files}
    for future in tqdm(as_completed(futures), total=len(futures), desc="Processing files"):
        result = future.result()
        if result:
            results.append(result)

# Save results to CSV
summary_df = pd.DataFrame(results)
summary_df.sort_values(by="n_extreme_tracks", ascending=False, inplace=True)
summary_df.to_csv(output_summary, index=False)
print(f"\nSummary saved to: {output_summary}")


set_trace()

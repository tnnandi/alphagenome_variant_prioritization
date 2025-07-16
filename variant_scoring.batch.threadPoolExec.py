import os
import re
import threading
import pandas as pd
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from concurrent.futures import ThreadPoolExecutor, as_completed # for parallel execution

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# Get API key from environment variable
API_KEY = os.getenv('ALPHAGENOME_API_KEY')
if not API_KEY:
    raise ValueError("ALPHAGENOME_API_KEY environment variable is not set. Please set it with your API key.")

# Configuration
INPUT_FILE = "Data_S1.xlsx" # path to the MVP study file downloaded from https://datadryad.org/dataset/doi:10.5061/dryad.zgmsbcck4
PIP_THRESHOLD = 0.95 # threshold for filtering variants by PIP from the MVP study
OUTPUT_DIR = "predictions" # directory to save the predictions to
MAX_WORKERS = 25  # Number of parallel threads
SEQUENCE_SIZE = '1MB'  # Options: "2KB", "16KB", "100KB", "500KB", "1MB" (length of the sequence around the variant to use for scoring)

LOAD_ALL_VARIANTS = False # set to True to load all variants from the MVP study, otherwise load the filtered variants from the filtered_variants.xlsx file
DOWNLOAD_PREDICTIONS = True # set to True to save the predictions to csv files


os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load variant data 
if LOAD_ALL_VARIANTS: # only once, to generate the filtered_variants.xlsx file
    df = pd.read_excel(INPUT_FILE, engine="openpyxl", header=1)
    # column names:
    """ Trait - The trait for which the signal was mapped
    Category - The parent category of the phenotype
    Description - Long-form description of the phenotype
    MVP ID - The internal MVP ID of the fine-mapped variant
    RSID - The rsID of the fine-mapped variant
    BP - The base pair position of the fine-mapped variant in hg19
    BP38 - The base pair position of the fine-mapped variant in GRCh38
    VEP Annotation - The most severe variant consequence annotated by the Variant Effect Predictor (VEP)
    Locus - The full locus range in which the signal was mapped
    Merged Signal - The number of the signal within the locus. The numbers are used to distinguish the multiple signals at a locus across the populations
    Population - The population for which the variant was fine-mapped. Variants may be mapped for a single trait in multiple populations and appear on multiple lines
    Population Signal - The number of the signal within the relevant population
    EAF Population - The effect allele frequency of the variant within the mapped population
    Beta Population - The GWAS variant beta within the mapped population
    SE Population - The GWAS variant standard error within the mapped population
    N Population - The size of the mapped population
    Overall PIP - The PIP of the variant within the mapped population summed over all credible sets
    CS-Level PIP - The PIP of the variant within the mapped population for the credible set in which it was mapped
    mu - The variant effect of the mapped variant under the SuSiE fine-mapping model
    mu2 - The second moment of the variant under the SuSiE fine-mapping model
    CS log(Bayes Factor) - The log of the Bayes Factor for the variant under the SuSiE model
    Previously Unidentified if High PIP - If the overall PIP > 0.95, this column describes whether the variant was previously reported """

    df["Overall PIP"] = pd.to_numeric(df["Overall PIP"], errors="coerce")

    # Plot distribution
    plt.figure(figsize=(10, 5))
    plt.hist(df["Overall PIP"], bins=20, color="gray", edgecolor="black")
    plt.xlabel("Overall PIP")
    plt.ylabel("Number of Variants")
    plt.xlim(0, 1)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.title("Distribution of Overall PIPs")
    plt.savefig("overall_pip_distribution.png")

    filtered_df = df[df["Overall PIP"] > PIP_THRESHOLD] # filter variants by PIP threshold
    filtered_df.to_excel("filtered_variants.xlsx", index=False)
    print(f"Filtered variants saved to filtered_variants.xlsx with {len(filtered_df)} entries.")
else:
    filtered_df = pd.read_excel("filtered_variants.xlsx", engine="openpyxl")
    print(f"Loaded {len(filtered_df)} filtered variants.")

# Setup model and scorers 
organism = dna_client.Organism.HOMO_SAPIENS
all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS

# choose which tracks to use for scoring
scorer_selections = {
    'rna_seq': True,
    'cage': True,
    'procap': True,
    'atac': True,
    'dnase': True,
    'chip_histone': True,
    'chip_tf': True,
    'polyadenylation': True,
    'splice_sites': True,
    'splice_site_usage': True,
    'splice_junctions': True,
}

selected_scorers = [
    all_scorers[key]
    for key in all_scorers
    if scorer_selections.get(key.lower(), False)
]
# remove scorers that are not supported for the organism
unsupported_scorers = [
    scorer for scorer in selected_scorers
    if (
        organism.value not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
        or (scorer.requested_output == dna_client.OutputType.PROCAP and organism == dna_client.Organism.MUS_MUSCULUS)
    )
]

# remove unsupported scorers
for scorer in unsupported_scorers:
    selected_scorers.remove(scorer)

# get the sequence length for the scoring
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f'SEQUENCE_LENGTH_{SEQUENCE_SIZE}']

# Worker function to score a variant
lock = threading.Lock() # lock to avoid race conditions when saving the predictions to csv files

def score_variant_worker(idx, row, api_key, selected_scorers, sequence_length, organism):
    try:
        dna_model = dna_client.create(api_key) # Initialize the API client

        mvp_id = row["MVP ID"] # MVP ID of the variant
        bp38 = row["BP38"] # base pair position of the variant in GRCh38
        population = row["Population"] # population of the variant
        chrom, _, ref, alt = mvp_id.split(":") # split the MVP ID into chromosome, position, reference base, and alternate base
        chrom = "chr" + chrom # add "chr" to the chromosome name
        pos = int(bp38) # use position from GRCh38 (used in alphagenome) and not GRCh37 as the MVP IDs use

        # create a variant object
        variant = genome.Variant(
            chromosome=chrom,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt,
        )
        # create a reference interval object
        interval = variant.reference_interval.resize(sequence_length)

        # score the variant
        scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=selected_scorers,
        )

        # create a tidy dataframe of the scores
        df_scores = variant_scorers.tidy_scores(scores)

        if DOWNLOAD_PREDICTIONS:
            safe_id = re.sub(r'[:>]', '_', mvp_id) # replace colons and greater than signs with underscores to avoid issues with filenames
            filename = f"chr_{safe_id}_{population}_{idx}_scores.csv" # create a filename for the scores
            output_path = os.path.join(OUTPUT_DIR, filename) # create a path to the output file
            with lock: # lock to avoid race conditions when saving the predictions to csv files
                tqdm.write(f"Saving scores to {output_path} with {len(df_scores)} tracks")
            df_scores.to_csv(output_path, index=False) # lock isn't needed as we are writing to a separate file for each variant

        return idx, "Success"

    except Exception as e:
        return idx, f"Error: {mvp_id} → {str(e)}"

# Parallel execution 
print(f"Scoring {len(filtered_df)} variants using {MAX_WORKERS} threads...")

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = {
        executor.submit(
            score_variant_worker, idx, row, API_KEY, selected_scorers, sequence_length, organism
        ): idx
        for idx, row in filtered_df.iterrows()
    }

    for future in tqdm(as_completed(futures), total=len(futures), desc="Scoring variants (parallel)"):
        idx, status = future.result()
        if "Error" in status:
            print(status)

print("Scoring complete.")

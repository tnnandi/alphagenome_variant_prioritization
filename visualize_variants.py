# code to visualize the variants
from alphagenome import colab_utils
from alphagenome.data import gene_annotation, genome, track_data, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pdb import set_trace

variant_file = "predictions/chr_14_101293528_A_C_EUR_15361_scores.csv" 

dna_model = dna_client.create(API_KEY)

# Load metadata objects for human.
output_metadata = dna_client.create(API_KEY).output_metadata(
    organism=dna_client.Organism.HOMO_SAPIENS
)

########### only for gene/transcript annotation and not for variant prediction/visualization
# Load gene annotations (from GENCODE).
gtf = pd.read_feather(
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)

# Filter to protein-coding genes and highly supported transcripts.
gtf_transcript = gene_annotation.filter_transcript_support_level(
    gene_annotation.filter_protein_coding(gtf), ['1']
)

# Extractor for identifying transcripts in a region.
transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)

# Also define an extractor that fetches only the longest transcript per gene.
gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(
    gtf_transcript
)
longest_transcript_extractor = transcript.TranscriptExtractor(
    gtf_longest_transcript
)

############

# Read the variant scores CSV file
scores_df = pd.read_csv(variant_file)

# Calculate absolute quantile scores and get top tracks
scores_df['abs_quantile_score'] = scores_df['quantile_score'].abs()
top_scores = scores_df.nlargest(10, 'abs_quantile_score')

print("Top 10 tracks by absolute quantile score:")
for idx, row in top_scores.iterrows():
    print(f"Score: {row['quantile_score']:.6f}, Output: {row['output_type']}, "
          f"Ontology: {row['ontology_curie']}, Biosample: {row['biosample_name']}")

# Extract unique ontology terms and output types from top scores
top_ontology_terms = top_scores['ontology_curie'].unique().tolist()
top_output_types = top_scores['output_type'].unique().tolist()

print(f"\nUnique ontology terms: {top_ontology_terms}")
print(f"Unique output types: {top_output_types}")

# Define interval to make predictions for (using the variant coordinates from CSV)
# Parse the scored_interval column which has format "chr14:100302903-101351479:."
scored_interval = scores_df['scored_interval'].iloc[0]
parts = scored_interval.split(':')
chromosome = parts[0]
start_end = parts[1].split('-')
interval_start = int(start_end[0])
interval_end = int(start_end[1].split('.')[0])  # Remove the trailing ':.'
interval_center = (interval_start + interval_end) // 2
print(f"Chromosome: {chromosome}, Interval start: {interval_start}, Interval end: {interval_end}, Interval center: {interval_center}")

# Create interval around the variant position
interval = genome.Interval(chromosome, interval_center - 50000, interval_center + 50000).resize(
    dna_client.SEQUENCE_LENGTH_1MB
)

print(f"Using interval: {interval}")

# Convert output types to DNA client output types
output_type_mapping = {
    'RNA_SEQ': dna_client.OutputType.RNA_SEQ,
    'CAGE': dna_client.OutputType.CAGE,
    'ATAC': dna_client.OutputType.ATAC,
    # 'CHIP': dna_client.OutputType.CHIP
}

requested_outputs = set()
for output_type in top_output_types:
    if output_type in output_type_mapping:
        requested_outputs.add(output_type_mapping[output_type])

print(f"Requested outputs: {requested_outputs}")

# Make predictions for the top tracks
output = dna_model.predict_interval(
    interval=interval,
    requested_outputs=requested_outputs,
    ontology_terms=top_ontology_terms,
)

# Extract the longest transcripts per gene for this interval.
longest_transcripts = longest_transcript_extractor.extract(interval)

# Build plot components
plot_components_list = [
    plot_components.TranscriptAnnotation(longest_transcripts),
]

# Add tracks for each output type
if hasattr(output, 'rna_seq') and output.rna_seq is not None:
    plot_components_list.append(
        plot_components.Tracks(
            tdata=output.rna_seq,
            ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
        )
    )

if hasattr(output, 'cage') and output.cage is not None:
    plot_components_list.append(
        plot_components.Tracks(
            tdata=output.cage,
            ylabel_template='CAGE: {biosample_name} ({strand})\n{name}',
        )
    )

if hasattr(output, 'atac') and output.atac is not None:
    plot_components_list.append(
        plot_components.Tracks(
            tdata=output.atac,
            ylabel_template='ATAC: {biosample_name} ({strand})\n{name}',
        )
    )


# Build plot
plot = plot_components.plot(
    plot_components_list,
    interval=interval,
    title=f'Predicted tracks for variant chr14:100827191:A>C\nTop absolute quantile scores from {len(top_scores)} tracks',
)

plot.savefig("plot_top_quantile_scores.png", dpi=300, bbox_inches='tight')
print("Plot saved as plot_top_quantile_scores.png")

# Print gene names in the region
gene_names = [t.info['gene_name'] for t in longest_transcripts if t.strand == '+']
print(f"Genes in the region: {gene_names}")

# Define the variant of interest
# get variant string from the variant_id column in the scores_df
variant_string = scores_df['variant_id'].iloc[0]    
variant = genome.Variant.from_str(variant_string)
print(f"Variant: {variant}")

set_trace()
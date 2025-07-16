# code to visualize the variants
import os
from alphagenome import colab_utils
from alphagenome.data import gene_annotation, genome, track_data, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pdb import set_trace
from datetime import datetime

date_time = datetime.now().strftime("%Y%m%d_%H%M%S")

# Get API key from environment variable
API_KEY = os.getenv('ALPHAGENOME_API_KEY')
if not API_KEY:
    raise ValueError("ALPHAGENOME_API_KEY environment variable is not set. Please set it with your API key.")

# choose the variant file to visualize
variant_file = "predictions/chr_14_101293528_A_C_EUR_15361_scores.csv" # the top variant in the variant_quantile_summary_parallel.csv file
variant_file = "predictions/chr_17_7571752_T_G_EUR_15820_scores.csv" 

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

#################################################################################################

# Read the scores for different tracks from the CSV file
scores_df = pd.read_csv(variant_file)

# Calculate absolute quantile scores and get top tracks
scores_df['abs_quantile_score'] = scores_df['quantile_score'].abs()
top_scores = scores_df.nlargest(10, 'abs_quantile_score')

print("Top 10 tracks by absolute quantile score:")
print("=" * 80)
for i, (idx, row) in enumerate(top_scores.iterrows(), 1):
    print(f"{i:2d}. Score: {row['quantile_score']:.6f} (abs: {row['abs_quantile_score']:.6f})")
    print(f"    Output: {row['output_type']}")
    print(f"    Ontology: {row['ontology_curie']}")
    print(f"    Biosample: {row['biosample_name']}")
    print(f"    Track: {row['track_name']} ({row['track_strand']})")
    print()

# Extract unique ontology terms and output types from top scores
top_ontology_terms = top_scores['ontology_curie'].unique().tolist()
top_output_types = top_scores['output_type'].unique().tolist()

print(f"\nUnique ontology terms: {top_ontology_terms}")
print(f"Unique output types: {top_output_types}")

# Define the variant of interest
# get variant string from the variant_id column in the scores_df
# Note that the file names follow the MVP ID (based on hg19), whereas the variant string is based on hg38
variant_string = scores_df['variant_id'].iloc[0]    
variant = genome.Variant.from_str(variant_string)
print(f"Variant: {variant}")

# Create interval around the variant position for transcript annotation
variant_center = variant.position
interval = genome.Interval(variant.chromosome, variant_center - 50000, variant_center + 50000).resize(
    dna_client.SEQUENCE_LENGTH_1MB
) # resize extends the interval to 1MB around the variant

print(f"Using interval: {interval}")

# Convert output types to DNA client output types (to be used as "requested_outputs" in predict_variant)
output_type_mapping = {
    'RNA_SEQ': dna_client.OutputType.RNA_SEQ,
    'CAGE': dna_client.OutputType.CAGE,
    'DNASE': dna_client.OutputType.DNASE,
    'ATAC': dna_client.OutputType.ATAC,
    'CHIP_HISTONE': dna_client.OutputType.CHIP_HISTONE,
    'CHIP_TF': dna_client.OutputType.CHIP_TF,
    'PROCAP': dna_client.OutputType.PROCAP,
    'CONTACT_MAPS': dna_client.OutputType.CONTACT_MAPS,
    'SPLICE_SITES': dna_client.OutputType.SPLICE_SITES,
    'SPLICE_SITE_USAGE': dna_client.OutputType.SPLICE_SITE_USAGE,
    'SPLICE_JUNCTIONS': dna_client.OutputType.SPLICE_JUNCTIONS,
}

requested_outputs = set()
for output_type in top_output_types:
    if output_type in output_type_mapping:
        requested_outputs.add(output_type_mapping[output_type])

print(f"Requested outputs: {requested_outputs}")

# Make variant predictions for the top tracks
output = dna_model.predict_variant(
    variant=variant,
    interval=interval,
    requested_outputs=requested_outputs,
    ontology_terms=top_ontology_terms,
)

# Extract the longest transcripts per gene for this interval.
longest_transcripts = longest_transcript_extractor.extract(interval)

# --- Filter output to only top 10 tracks ---

# Get the set of (track_name, track_strand) for the top tracks
top_track_keys = set(zip(top_scores['track_name'], top_scores['track_strand']))
print(f"Top track keys to filter for: {top_track_keys}")

# Print detailed info about top scores
print("\nDetailed top scores info:")
for i, (idx, row) in enumerate(top_scores.iterrows(), 1):
    print(f"{i}. Track: '{row['track_name']}' | Strand: '{row['track_strand']}' | Score: {row['quantile_score']:.6f}")

def filter_trackdata_to_top(trackdata):
    if trackdata is None:
        print("No trackdata available")
        return None
    
    print(f"Available tracks in metadata:")
    for _, row in trackdata.metadata.iterrows():
        print(f"  - '{row['name']}' ({row['strand']})")
    
    # Filter the metadata DataFrame to only rows in top_track_keys
    mask = [
        (row['name'], row['strand']) in top_track_keys
        for _, row in trackdata.metadata.iterrows()
    ]
    
    print(f"Found {sum(mask)} matching tracks out of {len(mask)} total")
    
    # If none match, return None
    if not any(mask):
        print("No matching tracks found")
        return None
    
    # Filter values and metadata
    filtered_metadata = trackdata.metadata[mask].reset_index(drop=True)
    filtered_values = trackdata.values[:, mask] if trackdata.values.ndim == 2 else trackdata.values[mask] # to handle chromatin tracks too
    
    print(f"Successfully filtered to {len(filtered_metadata)} tracks")
    
    # Rebuild TrackData object
    from alphagenome.data.track_data import TrackData
    return TrackData(
        values=filtered_values,
        metadata=filtered_metadata,
        resolution=trackdata.resolution,
        interval=trackdata.interval,
        uns=trackdata.uns
    )

# Define the colors for REF and ALT predictions.
ref_alt_colors = {'REF': 'dimgrey', 'ALT': 'red'}

# --- Build plot components with REF/ALT tracks for top 10 ---
plot_components_list = [
    plot_components.TranscriptAnnotation(longest_transcripts),
]
# set_trace()
# Filter and add RNA-seq tracks if present
if hasattr(output.reference, 'rna_seq') and output.reference.rna_seq is not None:
    print("RNA-seq tracks found")
    ref_rna_seq = filter_trackdata_to_top(output.reference.rna_seq)
    alt_rna_seq = filter_trackdata_to_top(output.alternate.rna_seq)
    
    if ref_rna_seq is not None and alt_rna_seq is not None:
        plot_components_list.append(
            plot_components.OverlaidTracks(
                tdata={
                    'REF': ref_rna_seq,
                    'ALT': alt_rna_seq,
                },
                colors=ref_alt_colors,
                ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
            )
        )
    else:
        print("No RNA-seq tracks found")
# set_trace()
# Filter and add CAGE tracks if present
if hasattr(output.reference, 'cage') and output.reference.cage is not None:
    ref_cage = filter_trackdata_to_top(output.reference.cage)
    alt_cage = filter_trackdata_to_top(output.alternate.cage)
    
    if ref_cage is not None and alt_cage is not None:
        plot_components_list.append(
            plot_components.OverlaidTracks(
                tdata={
                    'REF': ref_cage,
                    'ALT': alt_cage,
                },
                colors=ref_alt_colors,
                ylabel_template='CAGE: {biosample_name} ({strand})\n{name}',
            )
        )
    else:
        print("No CAGE tracks found")
# Filter and add ATAC tracks if present
if hasattr(output.reference, 'atac') and output.reference.atac is not None:
    ref_atac = filter_trackdata_to_top(output.reference.atac)
    alt_atac = filter_trackdata_to_top(output.alternate.atac)
    
    if ref_atac is not None and alt_atac is not None:
        plot_components_list.append(
            plot_components.OverlaidTracks(
                tdata={
                    'REF': ref_atac,
                    'ALT': alt_atac,
                },
                colors=ref_alt_colors,
                ylabel_template='ATAC: {biosample_name} ({strand})\n{name}',
            )
        )
    else:
        print("No ATAC tracks found")
# Build plot with variant annotation
plot = plot_components.plot(
    plot_components_list,
    annotations=[plot_components.VariantAnnotation([variant])],
    interval=interval,
    title=f'Effect of variant {variant} on predicted tracks\nTop 10 absolute quantile scores',
)

plot.savefig(f"plot_variant_effect_top10_{date_time}.png", dpi=300, bbox_inches='tight')
print(f"Plot saved as plot_variant_effect_top10_{date_time}.png")


# # Plot the difference between REF and ALT for the top track
# delta_values = alt_rna_seq.values - ref_rna_seq.values
# delta_track = track_data.TrackData(
#     values=delta_values,
#     metadata=ref_rna_seq.metadata,
#     resolution=ref_rna_seq.resolution,
#     interval=ref_rna_seq.interval
# )



# Plot in a smaller region around the variant
zoom_interval = genome.Interval(
    variant.chromosome,
    variant_center - 5_000,
    variant_center + 5_000
)

# Extract the longest transcripts in the zoomed region
zoom_longest = longest_transcript_extractor.extract(zoom_interval)

# Build the plot with the zoomed interval
plot = plot_components.plot(
    # [
    #     plot_components.TranscriptAnnotation(zoom_longest),
    #     plot_components_list
    # ],
    plot_components_list,
    annotations=[plot_components.VariantAnnotation([variant])],
    interval=zoom_interval,      # <-- only this range will be drawn
    title=f'Zoom ±5 kb around {variant}'
)

plot.savefig("zoom_5kb.png", dpi=300)

# Print gene names in the region
gene_names = [t.info['gene_name'] for t in longest_transcripts if t.strand == '+']
print(f"Genes in the region: {gene_names}")

# Print summary of top 10 tracks
print("\n" + "=" * 80)
print("SUMMARY OF TOP 10 TRACKS")
print("=" * 80)
print(f"Variant: {variant}")
print(f"Total tracks analyzed: {len(scores_df)}")
print(f"Top 10 tracks by absolute quantile score:")
print()

# Create a summary table
summary_data = []
for i, (idx, row) in enumerate(top_scores.iterrows(), 1):
    summary_data.append({
        'Rank': i,
        'Score': f"{row['quantile_score']:.6f}",
        'Abs_Score': f"{row['abs_quantile_score']:.6f}",
        'Output_Type': row['output_type'],
        'Biosample': row['biosample_name'],
        'Track': f"{row['track_name']} ({row['track_strand']})"
    })

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))

print("\n" + "=" * 80)
print("Analysis complete!")

set_trace()
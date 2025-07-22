import re
import time
from collections import defaultdict
from tools import generate_bed

from intervaltree import IntervalTree, Interval


def get_walks(gfa_filepath, seg_id_lookup, reference_sample):
    '''
    Maps Walk paths for a GFA file

    :param gfa_filepath: filepath of our target GFA file
    :param seg_id_lookup: A dictionary mapping each segment IDs to segment information.
    :param reference_sample: The sample ID of the reference genome to exclude from walks.
    :return walks: A defaultdict where keys are sample IDs (str), and values are another defaultdict.
                   The inner defaultdict's keys are sequence IDs (str, e.g., 'chr1#hap1'), and values
                   are IntervalTrees. Each IntervalTree contains Intervals representing segments
                   along that specific walk, with their coordinates (0-based, half-open) and info.
    '''
    print('Parsing target Walk paths...')
    time.sleep(3)

    walks = defaultdict(lambda: defaultdict(IntervalTree))

    with (open(gfa_filepath, 'r') as f):
        for line in f:
            if not line.startswith('W'):
                continue

            fields = line.strip().split('\t')

            if len(fields) < 7:
                print(f"Malformed W-line (expected at least 7 fields): {line.strip()}. Skipping.")
                continue

            sample_id = fields[1]

            # Do not map reference genome
            if sample_id == reference_sample:
                continue

            sequence_id = fields[3]  # This is the sequence ID of our sample

            walk_start_field = fields[4]
            walk_end_field = fields[5]

            current_position = 0  # 0-based current position on the walk

            walk_path_string = fields[6]  # The actual path string like '>123<456'

            # If walk_start is specified, initialize current_position with it (assuming 0-based)
            if walk_start_field != "*":
                try:
                    current_position = int(walk_start_field)
                except ValueError:
                    print(
                        f"Malformed field (expecting an integer for walk_start): {walk_start_field}. Skipping walk for {sample_id}/{sequence_id}.")
                    continue  # Skip this walk if start is malformed

            segments = re.findall(r'([<>])(\d+)', walk_path_string)

            if not segments:
                print(
                    f"Warning: W-line for {sample_id}/{sequence_id} has no valid segments in path string '{walk_path_string}'. Skipping this walk.")
                continue

            for orientation_char, seg_id in segments:
                orientation = '+' if orientation_char == '>' else '-'

                if seg_id not in seg_id_lookup:
                    print(
                        f"Error: Segment '{seg_id}' in walk for {sample_id}/{sequence_id} not found in S-line data (seg_id_lookup). Skipping this entire walk.")
                    break  # Break out of inner loop to skip the walk

                seg_info_from_lookup = seg_id_lookup[seg_id]
                segment_length = seg_info_from_lookup['segment_length']

                # Store information about this segment instance as it appears in the walk
                segment_info_in_walk = {
                    'segment_id': seg_id,
                    'segment_length': segment_length,
                    'orientation': orientation,
                    'current_start_position_on_walk': current_position,  # 0-based inclusive
                    'current_end_position_on_walk': current_position + segment_length  # 0-based exclusive
                }

                # Add this segment instance to the IntervalTree for the current sample's walk
                walks[sample_id][sequence_id].add(
                    Interval(current_position, current_position + segment_length, segment_info_in_walk)
                )
                current_position += segment_length

    print('Successfully parsed Walk paths.\n')
    time.sleep(3)
    return walks


def map_genes_to_target_walks(gene_ref_mappings, all_walks_data, seg_lookup):
    """
    Maps gene features from reference segments onto specific target genome walks.
    Identifies basic variants like novel insertions.

    :param gene_ref_mappings: List of dictionaries, output from map_genes_to_segments (0-based, half-open coordinates).
    :param all_walks_data: Dictionary, output from get_walks (sample_id -> contig_id -> IntervalTree).
    :param seg_id_lookup: Dictionary, mapping segment_id to its full info (from segment_mapping).
    :param reference_sample_id_to_exclude: The sample ID representing the reference genome to skip.
    :return: A defaultdict (sample_id -> contig_id -> list of mapped_gene_info dicts).
    """
    mapped_genes_on_target_walks = defaultdict(lambda: defaultdict(list))

    print("Starting gene mapping onto target genome walks...")
    time.sleep(3)

    # Iterate through each sample's walks
    for sample_id, sequence_ids in all_walks_data.items():
        print(f"  Processing walks for sample: {sample_id}")

        # Iterate through each sequence (haplotype path) within the sample
        for sequence_id, segments_in_walk in sequence_ids.items():  # segments_in_walk is an IntervalTree
            print(f"    Processing sequence: {sequence_id}")

            segments_by_id_in_current_walk = defaultdict(list)
            for interval_in_walk in segments_in_walk:
                segments_by_id_in_current_walk[interval_in_walk.data['segment_id']].append(interval_in_walk.data)

            # For each gene annotation mapped on a reference segment:
            for gene_feature_mapping in gene_ref_mappings:
                ref_segment_id = gene_feature_mapping['segment_id']
                original_ref_chrom = gene_feature_mapping['sequence_id']
                original_gene_start = gene_feature_mapping.get('original_feat_start')
                original_gene_end = gene_feature_mapping.get('original_feat_end')

                gene_start_relative_to_ref_seg = \
                    gene_feature_mapping['overlap_start'] - gene_feature_mapping['seg_start_on_ref']
                gene_end_relative_to_ref_seg = \
                    gene_feature_mapping['overlap_end'] - gene_feature_mapping[
                        'seg_start_on_ref']

                ref_seg_original_info = seg_lookup.get(ref_segment_id)
                if not ref_seg_original_info:
                    print(
                        f"Warning: Reference segment '{ref_segment_id}' not found in seg_id_lookup. Skipping gene {gene_feature_mapping['gene']}.")
                    continue

                # Retrieve the full segment information for the reference segment
                found_segments_in_walk = segments_by_id_in_current_walk.get(ref_segment_id, [])

                if not found_segments_in_walk:
                    continue

                for walk_segment_info in found_segments_in_walk:
                    walk_seg_start = walk_segment_info['current_start_position_on_walk']
                    walk_seg_end = walk_segment_info['current_end_position_on_walk']
                    walk_seg_orientation = walk_segment_info['orientation']

                    # If there is no sequence_id, this means the segment is a variant
                    has_tags = seg_lookup[ref_segment_id]['sequence_id'] is not None

                    # --- Coordinate Translation based on Orientation ---
                    if walk_seg_orientation == '+':
                        # Gene's relative position maps directly (example, seg starts at 150, gene on reference starts
                        # at 50. Gene starts at 200 on target)
                        gene_start_on_target = walk_seg_start + gene_start_relative_to_ref_seg
                        gene_end_on_target = walk_seg_start + gene_end_relative_to_ref_seg
                    else:
                        gene_start_on_target_unord = walk_seg_end - gene_end_relative_to_ref_seg
                        gene_end_on_target_unord = walk_seg_end - gene_start_relative_to_ref_seg

                        # Ensure start < end for 0-based, half-open format
                        gene_start_on_target = min(gene_start_on_target_unord, gene_end_on_target_unord)
                        gene_end_on_target = max(gene_start_on_target_unord, gene_end_on_target_unord)

                    # --- Variant Labeling ---
                    variant_label = "No Variants"
                    if not has_tags:
                        variant_label = "Variant"

                    overlap_len_on_target = gene_end_on_target - gene_start_on_target
                    if overlap_len_on_target < 0:
                        print(
                            f"Warning: Calculated negative overlap length for gene {gene_feature_mapping['gene']} on target {sample_id}/{sequence_id}. Skipping.")
                        continue

                    # Store the mapped gene information
                    mapped_gene_info = {
                        'sample_id': sample_id,
                        'target_contig_id': sequence_id,
                        'gene_name': gene_feature_mapping['gene'],
                        'feature_type': gene_feature_mapping['feature_type'],
                        'start_on_target': gene_start_on_target,
                        'end_on_target': gene_end_on_target,
                        'original_ref_chrom': original_ref_chrom,
                        'original_gene_start': original_gene_start,
                        'original_gene_end': original_gene_end,
                        'mapped_segment_id': ref_segment_id,
                        'segment_orientation_in_walk': walk_seg_orientation,
                        'overlap_len_on_target': overlap_len_on_target,
                        'variant_label': variant_label
                    }
                    mapped_genes_on_target_walks[sample_id][sequence_id].append(mapped_gene_info)

    print("Finished gene mapping onto target genome walks.\n")

    for sample_id, sequence_ids in mapped_genes_on_target_walks.items():
        all_genes_for_sample = []
        for sequence_id, gene_list in sequence_ids.items():
            all_genes_for_sample.extend(gene_list)

        if all_genes_for_sample:
            output_filename = f"test_files/{sample_id}.mapped_genes_on_target_walks.bed"
            generate_bed(all_genes_for_sample, output_filename)
        else:
            print(f"No genes mapped for sample {sample_id}. Skipping file generation.")

    return mapped_genes_on_target_walks

import re
import time
from collections import defaultdict
from tools import generate_bed_like_file

from intervaltree import IntervalTree, Interval


def get_walks(gfa_filepath, seg_lookup, reference_sample):
    '''
    Maps Walk paths for a GFA file

    :param gfa_filepath: filepath of our target GFA file
    :param seg_lookup: A dictionary mapping each segment IDs to segment information.
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

                if seg_id not in seg_lookup:
                    print(
                        f"Error: Segment '{seg_id}' in walk for {sample_id}/{sequence_id} not found in S-line data (seg_id_lookup). Skipping this entire walk.")
                    break  # Break out of inner loop to skip the walk

                seg_info_from_lookup = seg_lookup[seg_id]
                segment_length = seg_info_from_lookup['segment_length']

                # Store information about this segment instance as it appears in the walk
                segment_info_in_walk = {
                    'segment_id': seg_id,
                    'segment_length': segment_length,
                    'orientation': orientation,
                    'current_start_position_on_walk': current_position,
                    'current_end_position_on_walk': current_position + segment_length
                }

                # Add this segment instance to the IntervalTree for the current sample's walk
                walks[sample_id][sequence_id].add(
                    Interval(current_position, current_position + segment_length, segment_info_in_walk)
                )
                current_position += segment_length
            break

    print('Successfully parsed Walk paths.\n')
    time.sleep(3)
    return walks




def map_genes_to_target_segments(gene_ref_mappings, all_walks_data, seg_lookup):
    gene_start_on_target_walk = gene_end_on_target_walk = None
    mapped_target_segments = defaultdict(lambda: defaultdict(list))

    for sample_id, walks_for_sample in all_walks_data.items():
        target_output = f"test_files/{sample_id}.gene_mappings.bed"

        for target_sequence_id, segments_in_this_walk in walks_for_sample.items():

            segments_by_id_in_current_walk = defaultdict(list)
            for interval_in_walk in segments_in_this_walk:
                segments_by_id_in_current_walk[interval_in_walk.data['segment_id']].append(interval_in_walk.data)

            for gene_annotation_on_reference in gene_ref_mappings:
                ref_segment_id = gene_annotation_on_reference['segment_id']
                ref_gene = gene_annotation_on_reference['gene']
                ref_feature_type = gene_annotation_on_reference['feature_type']
                ref_feature_id = gene_annotation_on_reference['feature_id']
                gene_start_relative_to_ref_seg = gene_annotation_on_reference['overlap_start'] - \
                                                 gene_annotation_on_reference['seg_start_on_ref']
                gene_end_relative_to_ref_seg = gene_annotation_on_reference['overlap_end'] - \
                                               gene_annotation_on_reference['seg_start_on_ref']

                ref_segments_on_walk = segments_by_id_in_current_walk.get(ref_segment_id, [])
                if not ref_segments_on_walk:
                    continue

                for walk_segment_info in ref_segments_on_walk:
                    walk_seg_start_pos_on_walk = walk_segment_info['current_start_position_on_walk']
                    walk_seg_end_pos_on_walk = walk_segment_info['current_end_position_on_walk']
                    walk_seg_orientation = walk_segment_info['orientation']

                    if walk_seg_orientation == '+':
                        gene_start_on_target_walk = walk_seg_start_pos_on_walk + gene_start_relative_to_ref_seg
                        gene_end_on_target_walk = walk_seg_start_pos_on_walk + gene_end_relative_to_ref_seg
                    elif walk_seg_orientation == '-':
                        gene_start_on_target_walk_unord = walk_seg_end_pos_on_walk - gene_end_relative_to_ref_seg
                        gene_end_on_target_walk_unord = walk_seg_end_pos_on_walk - gene_start_relative_to_ref_seg
                        gene_start_on_target_walk = min(gene_start_on_target_walk_unord, gene_end_on_target_walk_unord)
                        gene_end_on_target_walk = max(gene_start_on_target_walk_unord, gene_end_on_target_walk_unord)

                    if gene_start_on_target_walk >= gene_end_on_target_walk:
                        continue

                    # --- Proper variant check ---
                    variant = False
                    for ovl in segments_in_this_walk.overlap(gene_start_on_target_walk, gene_end_on_target_walk):
                        ovl_id = ovl.data['segment_id']

                        # Skip reference segment itself
                        if ovl_id == ref_segment_id:
                            continue

                        # Variant = any non-reference segment overlapping the gene
                        variant = True
                        break

                    variant_label = 'Variant' if variant else 'No Variant'

                    mapped_target_segments[sample_id][target_sequence_id].append({
                        'sample_id': sample_id,
                        'target_sequence_id': target_sequence_id,
                        'feature_id': ref_feature_id,
                        'gene_name': ref_gene,
                        'feature_type': ref_feature_type,
                        'start_on_target': gene_start_on_target_walk,
                        'end_on_target': gene_end_on_target_walk,
                        'original_ref_chrom': gene_annotation_on_reference['sequence_id'],
                        'original_gene_start': gene_annotation_on_reference.get('original_feat_start'),
                        'original_gene_end': gene_annotation_on_reference.get('original_feat_end'),
                        'mapped_segment_id': ref_segment_id,
                        'segment_orientation_in_walk': walk_seg_orientation,
                        'overlap_len_on_target': gene_end_on_target_walk - gene_start_on_target_walk,
                        'variant_label': variant_label
                    })

        generate_bed_like_file(target_output, mapped_target_segments[sample_id])

    return mapped_target_segments
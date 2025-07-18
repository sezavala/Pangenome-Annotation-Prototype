import re
import time
from collections import defaultdict

from intervaltree import IntervalTree, Interval


def segment_mapping(gfa_filepath):
    '''
    Parses a GFA file to create a mapping of segment IDs to their genomic information.

    :param gfa_filepath: Path to the GFA file (str).
    :return: A dictionary where keys are chromosomes (str) and values are
             IntervalTree object holding genomic information.
    '''
    print('Parsing Segments from .gfa file...')
    time.sleep(3)
    with open(gfa_filepath, 'r') as f:
        segment_tree = defaultdict(IntervalTree)
        segment_id_lookup = {}
        reference_sample = None

        for line in f:
            # Only search for segments
            if not line.startswith('S'):
                continue

            fields = line.split('\t')
            if len(fields) < 3:  # GFA1 'S' line requires at least 3 fields
                print(f"Malformed 'S' line (too few fields): {line}")
                continue

            segment_id = fields[1]

            if fields[2] == '*':
                length = None
            else:
                length = len(fields[2])

            sequence_id = None
            current_start = None

            # Fields 0 - 2 are required, so this will be searching for optional tags
            for field in fields[3:]:
                tags = field.split(':')
                if len(tags) < 3:
                    print(f"Skipping malformed tag in segment {segment_id}: {field}")
                    continue
                tagName = tags[0]
                type_code = tags[1]
                value = tags[2]

                try:
                    if tagName == 'SO' and type_code == 'i':
                        current_start = int(value)
                    elif tagName == 'SN' and type_code == 'Z':
                        # Assumes SN format like 'CHM13#0#chr1'
                        try:
                            SN_fileds = value.split('#')
                            sequence_id = SN_fileds[2]
                            if reference_sample is None:
                                reference_sample = SN_fileds[0]
                        except IndexError:
                            print(f"Unexpected SN tag format for segment {segment_id}: '{value}'. Expected 'X#Y#chrZ'.")
                    elif tagName == 'LN' and type_code == 'i':
                        # If sequence was '*', use this length
                        if length is None:
                            length = int(value)
                except ValueError:
                    print(f"Invalid value type for tag '{tagName}' in segment {segment_id}: '{value}'")

            if length is None:
                print(
                    f"Error: Could not determine length for segment {segment_id} in S-line. Skipping this segment. Line: {line.strip()}")
                continue

            segment_info = {
                'segment_id': segment_id,
                'segment_length': length,
                'sequence_id': sequence_id,
                'seg_start_on_sequence_id': current_start
            }

            if current_start is not None and sequence_id is not None:
                segment_tree[sequence_id].add(
                    Interval(current_start, current_start + length, segment_info)
                )
            else:
                print(f"Missing 'SO' or 'SN' tag for segment '{segment_id}'. Skipping. Line: {line}")

            segment_id_lookup[segment_id] = segment_info

        print('Successfully parsed Segments from .gfa file...\n')
        time.sleep(3)

        return segment_tree, segment_id_lookup, reference_sample


def annotation_mapping(gff_filepath):
    '''
    Parses a GFF file to create a mapping of feature IDs to their genomic information.

    :param gff_filepath: Path to the GFF file (str).:
    :return A dictionary where keys are chromosomes (str), along with feature names (str), and values are
             dictionaries with 'feat_start', 'feat_end'.
             Returns an empty dictionary if no valid annotations are found.:
    '''
    print('Parsing Annotations from .gff file...')
    time.sleep(3)

    annotation_map = defaultdict(lambda: defaultdict(list))

    with open(gff_filepath, 'r') as f:
        for line in f:
            # Skip header
            if line.startswith('#'):
                continue

            # Fields are all required, so we can just parse them from the GFF file
            fields = line.strip().split('\t')
            if len(fields) != 9:
                print(f"Malformed GFF file: {line}")

            sequence_id = fields[0]
            feature_type = fields[2]
            attributes = fields[8].split(';')
            gene = None

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                print(f'Invalid start or end coordinate format in line: {line}')
                continue

            for attribute in attributes:
                if '=' in attribute:
                    tag, value = attribute.split('=', 1)
                    if tag == 'gene':
                        gene = value
                        break

            # Check for gene tag in 9th field
            if gene is None:
                print(f"Missing gene attribute in annotation: {line}")
                continue

            annotation_info = {
                'feat_start': start,
                'feat_end': end,
                'feat_type': feature_type,
            }

            annotation_map[sequence_id][gene].append(annotation_info)

        print('Successfully parsed Annotations from .gff file...\n')
        time.sleep(3)

        return annotation_map


def generate_bed(mappings, output_path):
    '''
    Generates a BED file from a list of mapped gene coordinates.

    :param mappings: A list of dictionaries holding genomic information from annotation mapping.
    :param output_path: The path to the output BED file.
    :return:
    '''

    print(f'Generating BED file...')
    time.sleep(3)

    with open(output_path, 'w') as f:
        for feature in mappings:
            for field in feature.values():
                f.write(f"{field}\t")
            f.write('\n')

    print(f'Generated BED file to {output_path}\n')
    time.sleep(3)


def walk_paths(gfa_filepath, seg_id_lookup, reference_sample):
    '''
    Maps Walk paths for a GFA file

    :param gfa_filepath: filepath of our target GFA file
    :param seg_id_lookup: A dictionary mapping each segment IDs to segment information.
    :return:
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

            sequence_id = fields[3]

            walk_start = fields[4]
            walk_end = fields[5]

            current_position = 0

            walk = fields[6]

            if walk_start != "*":
                try:
                    walk_start = int(walk_start)
                    current_position = walk_start
                except ValueError:
                    print(f"Malformed field (expecting an integer): {walk_start}. Skipping.")

            if walk_end != "*":
                try:
                    walk_end = int(walk_end)
                except ValueError:
                    print(f"Malformed field (expecting an integer): {walk_end}. Skipping.")

            segments = re.findall(r'([<>])(\d+)', walk)

            if not segments:
                print(
                    f"Warning: W-line for {sample_id}/{sequence_id} has no valid segments in path string '{walk}'. Skipping this walk.")
                continue

            for orientation, seg_id in segments:
                orientation_char = '+' if orientation == '>' else '-'

                if seg_id not in seg_id_lookup:
                    print(
                        f"Error: Segment '{seg_id}' in walk for {sample_id}/{sequence_id} not found in S-line data (seg_id_lookup). Skipping this entire walk.")
                    break

                seg_info = seg_id_lookup[seg_id]
                segment_length = seg_info['segment_length']

                segment_info = {
                    'segment_length': segment_length,
                    'orientation': orientation_char,
                    'segment_id': seg_id,
                    'current_start_position_on_walk': current_position,
                    'current_end_position_on_walk': current_position + segment_length
                }

                walks[sample_id][sequence_id].add(Interval(current_position, current_position + segment_length, segment_info))
                current_position += segment_length

    print('Successfully parsed Walk paths.\n')
    time.sleep(3)
    return walks


def map_genes_to_segments(gff_filepath, gfa_filepath):
    '''
    Maps gene coordinates from our reference annotation onto reference segments

    :param gff_filepath: annotated genome reference (str)
    :param gfa_filepath: genome assembly (str)
    :return mappings: A list of each segment ID with its chromosome, feature type, gene start/end positions, and gene
    '''
    mappings = []
    segment_tree, seg_id_lookup, reference_sample = segment_mapping(gfa_filepath)
    annotation_map = annotation_mapping(gff_filepath)

    print("Mapping gene coordinates onto reference segments...")
    time.sleep(3)

    for sequence_id, genes in annotation_map.items():
        if sequence_id not in segment_tree:
            continue
        for gene, features in genes.items():
            for feature in features:
                feat_start = feature['feat_start']
                feat_end = feature['feat_end']
                feat_type = feature['feat_type']
                for overlaps in segment_tree[sequence_id][feat_start:feat_end]:
                    seg_start = overlaps.begin
                    seg_end = overlaps.end
                    seg_info = overlaps.data

                    overlap_start = max(seg_start, feat_start)
                    overlap_end = min(seg_end, feat_end)

                    mappings.append({
                        'sequence_id': sequence_id,
                        'gene': gene,
                        'feature_type': feat_type,
                        'original_feat_start': feat_start,
                        'original_feat_end': feat_end,
                        'segment_id': seg_info['segment_id'],
                        'seg_start': seg_start,
                        'seg_end': seg_end,
                        'overlap_start': overlap_start,
                        'overlap_end': overlap_end,
                        'overlap_len': overlap_end - overlap_start
                    })

    print("Successfully mapped gene coordinates to reference segments...\n")
    time.sleep(3)

    return mappings, seg_id_lookup, reference_sample


def map_genes_to_targets(walk_paths, segment_mappings):
    '''
    Maps gene coordinates from our reference annotation onto target segments

    :param walk_paths: sequence path of each target sample and sequence
    :param segment_mappings: annotated reference segment mappings
    :return:
    '''
    print("Mapping gene coordinates onto target segments...")
    time.sleep(3)

    print("Successfully mapped gene coordinates to target segments...")
    time.sleep(3)


if __name__ == '__main__':
    ref_loc = "test_files/chm13.MANE.gff3"
    target_loc = "test_files/chm13-kolf2.1j.gfa"
    bed_output_loc = ref_loc + ".bed"

    segment_mapping, seg_id_lookup, reference_sample = map_genes_to_segments(ref_loc, target_loc)
    generate_bed(segment_mapping, bed_output_loc)

    walk_paths = walk_paths(target_loc, seg_id_lookup, reference_sample)
    for sample_id in walk_paths:
        for sequence_id in walk_paths[sample_id]:
            print(f"{sample_id}/{sequence_id}")
    # target_mapping = map_genes_to_targets(walk_paths, segment_mapping)

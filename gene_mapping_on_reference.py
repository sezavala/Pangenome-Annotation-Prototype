import time
from collections import defaultdict
from tools import generate_bed

from intervaltree import Interval, IntervalTree


def segment_mapping(gfa_filepath):
    '''
    Parses a GFA file to create a mapping of segment IDs to their genomic information.

    :param gfa_filepath: Path to the GFA file (str).
    :return segment_tree: A dictionary mapping of segment IDS (key) to an IntervalTree (value) holding segment
    start and end coordinates as Intervals and their genomic information.
    :return segment_lookup: A dictionary mapping of all segment IDs to their genomic information.
    :return reference_sample: The sample ID of the reference genome in our pangenome.
    '''
    print('Parsing Segments from .gfa file...')
    time.sleep(3)

    with open(gfa_filepath, 'r') as f:
        segment_tree = defaultdict(IntervalTree)
        segment_lookup = {}
        reference_sample = None

        for line in f:
            # Only search for segments
            if not line.startswith('S'):
                continue

            # Ensure we have at least 3 fields (GFA1 standard)
            fields = line.split('\t')
            if len(fields) < 3:
                print(f"Malformed 'S' line (too few fields): {line}")
                continue

            # Segment ID for sequence
            segment_id = fields[1]

            # If sequence is included, extract length.
            # Else, search for LN tag
            if fields[2] == '*':
                segment_length = None
            else:
                segment_length = len(fields[2])

            sequence_id = None
            current_start = None

            # Fields 0 - 2 are required, so this will be searching for optional tags
            for field in fields[3:]:
                # Verify tag format
                tags = field.split(':')
                if len(tags) < 3:
                    print(f"Skipping malformed tag in segment {segment_id}: {field}")
                    continue

                tag_name = tags[0]
                type_code = tags[1]
                value = tags[2]

                try:
                    # Extract starting index (0-based index)
                    if tag_name == 'SO' and type_code == 'i':
                        current_start = int(value)
                    # Extract sequence name (Assumes SN format like 'CHM13#0#chr1)
                    elif tag_name == 'SN' and type_code == 'Z':
                        try:
                            sn_fields = value.split('#')
                            sequence_id = sn_fields[2]
                            if reference_sample is None:
                                reference_sample = sn_fields[0]
                        except IndexError:
                            print(f"Unexpected SN tag format for segment {segment_id}: '{value}'. Expected 'X#Y#chrZ'.")
                    elif tag_name == 'LN' and type_code == 'i':
                        # If sequence was '*', use this length
                        if segment_length is None:
                            segment_length = int(value)
                except ValueError:
                    print(f"Invalid value type for tag '{tag_name}' in segment {segment_id}: '{value}'")

            if segment_length is None:
                print(
                    f"Error: Could not determine length for segment {segment_id} in S-line. Skipping this segment. Line: {line.strip()}")
                continue

            segment_info = {
                'segment_id': segment_id,
                'segment_length': segment_length,
                'sequence_id': sequence_id,
                'segment_start_on_sequence_id': current_start
            }

            if current_start is not None and sequence_id is not None:
                segment_tree[sequence_id].add(
                    Interval(current_start, current_start + segment_length, segment_info)
                )
            else:
                print(f"Missing 'SO' or 'SN' tag for segment '{segment_id}'. Skipping. Line: {line}")

            segment_lookup[segment_id] = segment_info

        print('Successfully parsed Segments from .gfa file...\n')
        time.sleep(3)

        return segment_tree, segment_lookup, reference_sample


def annotation_mapping(gff_filepath):
    '''
    Parses a GFF file to create a mapping of feature IDs to their genomic information.

    :param gff_filepath: Path to the GFF file (str).:
    :return A dictionary where keys are chromosomes (str), along with gene name (str), and values are
             dictionaries with 'feature_start', 'feature_end', and 'feature_type'.
    '''
    print('Parsing Annotations from .gff file...')
    time.sleep(3)

    annotation_map = defaultdict(lambda: defaultdict(list))

    with open(gff_filepath, 'r') as f:
        for line in f:
            # Skip header
            if line.startswith('#'):
                continue

            # 9 fields are required (GFF standard)
            fields = line.strip().split('\t')
            if len(fields) != 9:
                print(f"Malformed GFF file: {line}")

            sequence_id = fields[0]
            feature_type = fields[2]
            attributes = fields[8].split(';')
            gene = None
            gene_id = None

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                print(f'Invalid start or end coordinate format in line: {line}')
                continue

            # Search for gene tags in 9th field
            for attribute in attributes:
                if '=' in attribute:
                    tag, value = attribute.split('=', 1)
                    if tag == 'gene':
                        gene = value
                    elif tag == 'ID':
                        gene_id = value
                    elif tag == 'Parent' and gene_id is None:
                        gene_id = value

            # Skip lines with no gene attributes
            if gene is None:
                print(f"Missing gene attribute in annotation: {line}")
                continue

            annotation_info = {
                'feature_start': start,
                'feature_end': end,
                'feature_type': feature_type,
                'feature_id': gene_id,
            }

            annotation_map[sequence_id][gene].append(annotation_info)

        print('Successfully parsed Annotations from .gff file...\n')
        time.sleep(3)

        return annotation_map


def map_genes_to_segments(segment_tree, annotation_map, reference_sample):
    '''
    Maps gene coordinates from our reference annotation onto reference segments

    :param segment_tree: A dictionary mapping of segment IDS (key) to an IntervalTree (value) holding segment
    start and end coordinates as Intervals and their genomic information.
    :param annotation_map: A dictionary where keys are chromosomes (str), along with gene name (str), and values are
    dictionaries with 'feature_start', 'feature_end', and 'feature_type'.
    :param reference_sample: The sample ID of the reference genome in our pangenome.
    :return mappings: A list of dictionaries holding the translated gene coordinates on our reference's segments
    and other genomic information.
    '''
    mappings = []
    target_output = "test_files/" + reference_sample + ".reference_gene_coordinate_mappings.bed"

    print("Mapping gene coordinates onto reference segments..\n")
    time.sleep(3)

    for sequence_id, genes in annotation_map.items():
        if sequence_id not in segment_tree:
            continue
        for gene, features in genes.items():
            for feature in features:
                feat_start = feature['feature_start']
                feat_end = feature['feature_end']
                feat_type = feature['feature_type']
                feature_id = feature['feature_id']

                feat_start_for_calc = feat_start - 1
                feat_end_for_calc = feat_end

                for overlaps in segment_tree[sequence_id][feat_start_for_calc:feat_end_for_calc]:
                    seg_start = overlaps.begin
                    seg_end = overlaps.end
                    seg_info = overlaps.data

                    overlap_start = max(seg_start, feat_start)
                    overlap_end = min(seg_end, feat_end)

                    mappings.append({
                        'sequence_id': sequence_id,
                        'gene': gene,
                        'feature_id': feature_id,
                        'feature_type': feat_type,
                        'original_feat_start': feat_start,
                        'original_feat_end': feat_end,
                        'segment_id': seg_info['segment_id'],
                        'seg_start_on_ref': seg_start,
                        'seg_end_on_ref': seg_end,
                        'overlap_start': overlap_start,
                        'overlap_end': overlap_end,
                        'overlap_len': overlap_end - overlap_start
                    })

    generate_bed(mappings, target_output)

    print("Successfully mapped gene coordinates to reference segments...\n")
    time.sleep(3)

    return mappings

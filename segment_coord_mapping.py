from collections import defaultdict

from intervaltree import IntervalTree, Interval


def segment_mapping(gfa_filepath):
    '''
    Parses a GFA file to create a mapping of segment IDs to their genomic information.

    :param gfa_filepath: Path to the GFA file (str).
    :return: A dictionary where keys are chromosomes (str) and values are
             IntervalTree object holding genomic information.
    '''
    with open(gfa_filepath, 'r') as f:
        segment_tree = defaultdict(IntervalTree)

        for line in f:
           # Skip header
            if line.startswith('H'):
                continue

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
            found_length_tag = False

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
                            sequence_id = value.split('#')[2]
                        except IndexError:
                            print(f"Unexpected SN tag format for segment {segment_id}: '{value}'. Expected 'X#Y#chrZ'.")
                    elif tagName == 'LN' and type_code == 'i':
                        # If sequence was '*', use this length
                        if length is None:
                            length = int(value)
                            found_length_tag = True
                except ValueError:
                    print(f"Invalid value type for tag '{tagName}' in segment {segment_id}: '{value}'")

            if length is None and not found_length_tag:
                print(f"Could not determine length for segment {segment_id}. Skipping.")
                continue

            if current_start is not None and sequence_id is not None:
                segment_info = {
                    'seg_start': current_start,
                    'seg_end': current_start + length,
                    'seg_length': length,
                    'seg_id': segment_id,
                }

                segment_tree[sequence_id].add(
                    Interval(current_start, current_start + length, segment_info)
                )
            else:
                print(f"Missing 'SO' or 'SN' tag for segment '{segment_id}'. Skipping. Line: {line}")

        return segment_tree

def annotation_mapping(gff_filepath):
    '''
    Parses a GFF file to create a mapping of feature IDs to their genomic information.

    :param gff_filepath: Path to the GFF file (str).:
    :return A dictionary where keys are chromosomes (str), along with feature names (str), and values are
             dictionaries with 'feat_start', 'feat_end'.
             Returns an empty dictionary if no valid annotations are found.:
    '''
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

        return annotation_map

def generate_bed(mappings, output_path):
    '''
    Generates a BED file from a list of mapped gene coordinates.

    :param mappings: A list of dictionaries holding genomic information from annotation mapping.
    :param output_path: The path to the output BED file.
    :return:
    '''

    print(f'Generating BED file...\n')

    with open(output_path, 'w') as f:
        for feature in mappings:
            for field in feature.values():
                f.write(f"{field}\t")
            f.write('\n')

    print(f'Generated BED file to {output_path}')

def map_genes_to_segments(gff_filepath, gfa_filepath):
    '''
    Maps gene coordinates from our reference annotation onto reference segments

    :param gff_filepath: annotated genome reference (str)
    :param gfa_filepath: genome assembly (str)
    :return mappings: A list of each segment ID with its chromosome, feature type, gene start/end positions, and gene
    '''
    mappings = []
    segment_tree = segment_mapping(gfa_filepath)
    annotation_map = annotation_mapping(gff_filepath)

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
                        'segment_id': seg_info['seg_id'],
                        'seg_start': seg_start,
                        'seg_end': seg_end,
                        'overlap_start': overlap_start,
                        'overlap_end': overlap_end,
                        'overlap_len': overlap_end - overlap_start
                    })

    return mappings

if __name__ == '__main__':
    ref_loc = "test_files/chm13.MANE.gff3"
    target_loc = "test_files/chm13-kolf2.1j.gfa"
    bed_output_loc = "test_files/chm13-kolf2.1j.bed"

    mappings = map_genes_to_segments(ref_loc, target_loc)
    generate_bed(mappings, bed_output_loc)
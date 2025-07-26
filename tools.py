import time


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

def generate_bed_like_file(output_filepath, mapped_genes_for_sample_data):
    """
    Generates a tab-separated BED-like file for mapped genes of a single sample, WITHOUT A HEADER.

    :param output_filepath: The full path to the output file (e.g., "test_files/sampleX.gene_mappings.bed").
    :param mapped_genes_for_sample_data: A defaultdict(list) where keys are target_sequence_id
                                         and values are lists of mapped gene info dictionaries for that sample.
    """
    with open(output_filepath, 'w') as outfile:
        for target_sequence_id, gene_list in mapped_genes_for_sample_data.items():
            for gene_info in gene_list:
                sample_id = gene_info['sample_id']
                feature_id = gene_info['feature_id']
                chrom = gene_info['target_sequence_id']
                start = str(gene_info['start_on_target'])
                end = str(gene_info['end_on_target'])
                name = gene_info['gene_name']
                orientation = gene_info['segment_orientation_in_walk']

                original_ref_chrom = gene_info['original_ref_chrom']
                original_gene_start = str(gene_info['original_gene_start'])
                original_gene_end = str(gene_info['original_gene_end'])
                mapped_segment_id = gene_info['mapped_segment_id']
                overlap_len_on_target = str(gene_info['overlap_len_on_target'])
                variant_label = gene_info['variant_label']
                feature_type = gene_info['feature_type']

                # Create the tab-separated string for the current gene entry
                line_data = [
                    sample_id, chrom, feature_id, start, end, name, feature_type, orientation,
                    original_ref_chrom, original_gene_start, original_gene_end,
                    mapped_segment_id, overlap_len_on_target, variant_label
                ]
                outfile.write("\t".join(line_data) + "\n")

    print(f"Generated BED-like file at: {output_filepath}")
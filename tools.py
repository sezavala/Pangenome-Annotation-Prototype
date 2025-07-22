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
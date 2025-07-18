import re
from collections import defaultdict

def walk_paths(gfa_filepath, seg_id_lookup):
    '''
    Maps Walk paths for a GFA file

    :param gfa_filepath: filepath of our target GFA file
    :return:
    '''
    walks = defaultdict(lambda: defaultdict(list))

    with open(gfa_filepath, 'r') as f:
        for line in f:
            if not line.startswith('W'):
                continue

            fields = line.strip().split('\t')

            if len(fields) < 7:
                print(f"Malformed W-line (expected at least 7 fields): {line.strip()}. Skipping.")
                continue

            sample_id = fields[1]
            sequence_id = fields[3]

            walk_start = fields[4]
            walk_end = fields[5]

            walk = fields[6]

            if walk_start != "*":
                try:
                    walk_start = int(walk_start)
                except ValueError:
                    print(f"Malformed field (expecting an integer): {walk_start}. Skipping.")

            if walk_end != "*":
                try:
                    walk_end = int(walk_end)
                except ValueError:
                    print(f"Malformed field (expecting an integer): {walk_end}. Skipping.")

            segments = re.findall(r'([<>])(\d+)', walk)
            segment_length = 0

            for segment in segments:
                segment_length += 1
                orientation = segment[0]
                segment_id = segment[1]

                segment_info = {
                    'seg_length': segment_length,
                    'orientation': orientation,
                    'segment_id': segment_id
                }

                walks[sample_id][sequence_id].append(segment_info)

    return walks


if __name__ == '__main__':
    target_loc = "test_files/chm13-kolf2.1j.gfa"

    walk_paths = walk_paths(target_loc)
    for sample_id in walk_paths:
        for sequence_id in walk_paths[sample_id]:
            print(f"{sample_id}\t{sequence_id}")
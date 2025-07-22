from gene_mapping_on_reference import segment_mapping, annotation_mapping, map_genes_to_segments
from gene_mapping_on_target import get_walks, map_genes_to_target_walks

if __name__ == '__main__':
    reference_annotation_location = "test_files/chm13.MANE.gff3"
    pangenome_graph_location = "test_files/chm13-kolf2.1j.gfa"

    # Step 1: Map gene coordinates onto reference genome segments and save annotations in a bed file

    # Save segments from reference genome into an Interval Tree for fast overlap search
    # Store all segments in a dictionary for ID lookups
    # Note reference sample ID to avoid when mapping targets
    segment_tree, seg_lookup, reference_sample = segment_mapping(pangenome_graph_location)

    # Save gene annotation information along with their start and end locations
    annotation_map = annotation_mapping(reference_annotation_location)

    # Map gene locations onto reference segments and generate BED-like file
    gene_coordinates_on_reference_mapping = map_genes_to_segments(segment_tree ,annotation_map, reference_sample)

    # Step 2: Map gene coordinates onto target genome(s) segments

    # Save all walks from our target samples into a map
    walks = get_walks(pangenome_graph_location, seg_lookup, reference_sample)
    gene_coordinates_on_targets_mapping = map_genes_to_target_walks(gene_coordinates_on_reference_mapping, walks, seg_lookup)
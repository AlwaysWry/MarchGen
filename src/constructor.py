# a module of build the final March test, using the sub-sequences generated in subseq_creator.py
from subseq_creator import *


def get_linked_sequence_intersection(init_0_pool, init_1_pool):
	seq_intersection = set()
	for seq in init_0_pool:
		if find_identical_objs(seq, init_1_pool, {'ass_init'}):
			seq_intersection.add(seq)

	return seq_intersection


def build_coverage_graph():
	pass


if __name__ == '__main__':
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'])
	flatten_SF_pool = flatten_sf_pool(filtered_SF_pool)
	filtered_unlinked_pool = (filter_redundant_2cF(filtered_SF_pool, classified_pool['2cF_nonCFds_included'],
												   classified_pool['2cF_CFds']))
	seq_pool = create_sequence_pool(flatten_SF_pool, filtered_unlinked_pool, classified_pool['2cF_CFds']['linked'])
	for seq in get_linked_sequence_intersection(seq_pool[0]['Init_0'], seq_pool[0]['Init_1']):
		print(seq.seq_text)

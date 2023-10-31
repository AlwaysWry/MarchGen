# a module of build the final March test, using the sub-sequences generated in subseq_creator.py
from subseq_creator import *

HIT = True
NOT_HIT = False


class CoverageVertex:
	"""data structure of vertices in coverage graph"""

	def __init__(self, prop_dict):
		"""
		:param prop_dict:
		index: the index of current vertex in coverage graph
		coverage: a list contains all sequences that this vertex covers
		"""
		self.__dict__.update(prop_dict)

	def get_march_sequence(self):
		if len(self.__dict__['coverage']) > 1:
			init_seq = self.__dict__['coverage'][0].seq_text[0]
			for seq in self.__dict__['coverage']:
				init_seq += seq.seq_text[1:]
		else:
			init_seq = self.__dict__['coverage'][0].seq_text

		if self.__dict__['coverage'][-1].detect_tag:
			init_seq += 'r' + init_seq[-1]

		return init_seq

# end of vertex definition


def get_linked_sequence_intersection(init_0_pool, init_1_pool):
	seq_intersection = set()
	redundant_0_seq = set()
	redundant_1_seq = set()

	for seq in init_0_pool:
		find_result = find_identical_objs(seq, init_1_pool, {'ass_init'})
		if not isinstance(find_result, int):
			seq_intersection.add(seq)
			redundant_0_seq.add(seq)
			redundant_1_seq.add(find_result)

	init_0_pool -= redundant_0_seq
	init_1_pool -= redundant_1_seq

	return seq_intersection


def find_terminal_seq(seq_intersection, inits):
	# find terminal sequence in the sequence linked_intersection set.
	# terminal sequences will be excluded from the shortest path search.
	# terminal sequence has to meet following conditions:
	# 1) start with a read operation
	# 2) have different start and end states;
	terminal_seq_pool = {}
	for init_key in inits:
		terminal_seq_pool[init_key] = set()

	for seq in seq_intersection:
		if seq.seq_text[1] != 'r':
			continue
		elif seq.seq_text[0] == seq.seq_text[-1]:
			continue
		else:
			key = 'Init_' + str(seq.seq_text[0])
			if key in terminal_seq_pool.keys():
				terminal_seq_pool[key].add(seq)

	return terminal_seq_pool


def filter_vertices(vertex_pool):
	nest_vertices = set(filter(lambda v: len(v.coverage) > 1, vertex_pool))
	unnest_vertices = vertex_pool - nest_vertices
	nest_texts = list(map(lambda v: v.get_march_sequence(), nest_vertices))
	redundant_vertices = set()
	for v_obj in unnest_vertices:
		for n_text in nest_texts:
			if v_obj.get_march_sequence() in n_text:
				redundant_vertices.add(v_obj)
	filtered_pool = vertex_pool - redundant_vertices

	return filtered_pool


def define_vertices(sequence_pool):
	vertex_pool = set()
	for seq in sequence_pool:
		if seq.detect_tag and seq.nest_tag != 'invalid':
			# check if the corresponding nest sequence exists, merge them into 1 vertex
			nest_init = lambda i: '0' if i == '1' else '1'
			nest_text = nest_init(seq.seq_text[0]) + seq.seq_text[1:]
			nest_vertex = list(filter(lambda v: nest_text in map(lambda s: s.seq_text, v.coverage), vertex_pool))
			if len(nest_vertex) > 0:
				if seq.nest_tag == 'donor':
					nest_vertex[0].coverage.insert(0, copy.deepcopy(seq))
				else:
					nest_vertex[0].coverage.append(copy.deepcopy(seq))
				continue

		vertex_dict = {'index': -1, 'coverage': [copy.deepcopy(seq)], 'diff': -1}
		vertex_pool.add(CoverageVertex(vertex_dict))

	return filter_vertices(vertex_pool)


def calculate_diff_value(chain, vertex):
	match_target = vertex.get_march_sequence()
	search_range = min(len(chain[2:]), len(match_target))
	for location in range(-search_range, -1, 2):
		if chain[location:] == match_target[:-location]:
			diff = len(match_target[-location:])
			return diff
	return len(match_target)


def get_vertex_winner(vertex_candidates):
	unnest_pool = set(filter(lambda v: v.coverage[0].nest_tag == 'invalid', vertex_candidates))
	nest_pool = vertex_candidates - unnest_pool

	def check_transition(vertex):
		seq = vertex.get_march_sequence()
		return seq[0] == seq[-1]

	if len(unnest_pool) > 1:
		secondary_pool = set(filter(lambda v: v.coverage[0].dr_tag, unnest_pool))
		if len(secondary_pool) > 1:
			tertiary_pool = set(filter(check_transition, secondary_pool))
		elif len(secondary_pool) == 0:
			tertiary_pool = set(filter(check_transition, unnest_pool - secondary_pool))
		else:
			return next(iter(secondary_pool))
	elif len(unnest_pool) == 0:
		secondary_pool = set(filter(lambda v: v.coverage[0].nest_tag == 'donor', nest_pool))
		if len(secondary_pool) > 1:
			tertiary_pool = set(filter(check_transition, secondary_pool))
		elif len(secondary_pool) == 0:
			tertiary_pool = set(filter(check_transition, nest_pool - secondary_pool))
		else:
			return next(iter(secondary_pool))
	else:
		return next(iter(unnest_pool))

	if len(tertiary_pool) > 0:
		return next(iter(tertiary_pool))
	else:
		return next(iter(secondary_pool - tertiary_pool))


def build_coverage_chain(chain, vertex_pool):
	for v_obj in vertex_pool:
		v_obj.diff = calculate_diff_value(chain, v_obj)

	min_diff = min(map(lambda v: v.diff, vertex_pool))
	vertex_candidates = set(filter(lambda v: v.diff == min_diff, vertex_pool))

	if len(vertex_candidates) > 1:
		vertex_winner = get_vertex_winner(vertex_candidates)
	else:
		vertex_winner = next(iter(vertex_candidates))

	chain_segment = vertex_winner.get_march_sequence()
	if vertex_winner.diff < len(chain_segment):
		chain_appendix = chain_segment[-vertex_winner.diff:]
	else:
		chain_appendix = 'w' + chain_segment

	return vertex_winner, chain_appendix


def terminal_decorator(chain):
	feature_dict = {}
	terminal_feature = chain[0] + chain[-1]
	# TODO: consider the terminal seq and the structure, test the TS-included and TS-excluded case
	pass


def construct_main_element(vertex_pool):
	vertex_candidate_pool = copy.deepcopy(vertex_pool)
	initial_vertex = get_vertex_winner(vertex_candidate_pool)
	vertex_candidate_pool -= initial_vertex
	coverage_chain = initial_vertex.get_march_sequence()

	while len(vertex_candidate_pool) > 0:
		build_result = build_coverage_chain(coverage_chain, vertex_candidate_pool)
		vertex_candidate_pool -= build_result[0]
		coverage_chain += build_result[1]
	pass


if __name__ == '__main__':
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_2cF_pool = (filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']))
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], filtered_2cF_pool)
	flat_SF_pool = flatten_sf_pool(filtered_SF_pool)

	seq_pool = create_sequence_pool(flat_SF_pool, filtered_2cF_pool, classified_pool['2cF_CFds']['linked'])
	linked_intersection = get_linked_sequence_intersection(seq_pool['linked']['Init_0'], seq_pool['linked']['Init_1'])
	# for vertex in define_vertices(linked_intersection):
	# 	print(vertex.get_initial_sequence())

	linked_intersection = linked_intersection
	linked_remainder_0 = seq_pool['linked']['Init_0'] - linked_intersection
	linked_remainder_1 = seq_pool['linked']['Init_1'] - linked_intersection

	terminal_seq = set(find_terminal_seq(linked_intersection, {'Init_0', 'Init_1'}))

	print(terminal_seq)

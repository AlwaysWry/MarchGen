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
		self.index = -1
		self.coverage = []
		self.__dict__.update(prop_dict)

	def get_initial_sequence(self):
		if len(self.coverage) > 1:
			init_seq = self.coverage[-1].seq_text
			if self.coverage[-1].detect_tag:
				init_seq += 'r' + init_seq[-1]
			return init_seq
		else:
			# the first seq must be '0' or '1'
			return self.coverage[0].seq_text


# end of vertex definition


class CoverageEdge:
	"""data structure of edges in coverage graph"""

	def __init__(self, prop_dict):
		"""
		:param prop_dict:
		index: the index of current edge in coverage graph
		weight: the weight of current edge
		operation_map: the operations that need to be added into ME for visiting current edge
		"""
		self.terminal = []
		self.weight = 0
		self.operation_map = []
		self.__dict__.update(prop_dict)


# end of edge definition


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


def nest_available_sequences(linked_seq_pool):
	return linked_seq_pool
	pass


def arbit_sequence_nest_condition(seq_segment, seq_obj):
	# check whether the new sequence forms a nest sensitization segment
	max_len = len(seq_obj.seq_text) - 1
	if len(seq_segment) > 2 * max_len:
		if (seq_segment[0] != seq_segment[-1]) and (seq_segment[1:max_len + 1] == seq_segment[max_len + 1:]):
			return NEST
	else:
		return NOT_NEST


def hit_vertex(seq, sequence_pool):
	# check if current sequence covers a new vertex
	coverage = []
	sequence_CFds_pool = set()
	sequence_nonCFds_pool = set()
	coverage_pool = set(map(lambda s_obj: s_obj.seq_text, sequence_pool))

	for seq_obj in sequence_pool:
		if seq_obj.detect_tag:
			sequence_nonCFds_pool.add(seq_obj)
		else:
			sequence_CFds_pool.add(seq_obj)

	if seq[-2] == 'r':
		# nonCFds case
		search_stop = max(len(max(coverage_pool, key=len)), len(seq[:-2]))
		for index in range(search_stop, 1, -2):
			for seq_obj in sequence_nonCFds_pool:
				if (seq[-index - 2:-2] == seq_obj.seq_text) and (seq_obj.seq_text not in coverage):
					coverage.append(seq_obj)
					break
			else:
				# continue if the inner loop is not break
				continue

			if (seq_obj.nest_tag == 'receiver') and arbit_sequence_nest_condition(seq[-index - 2:-2], seq_obj):
				# when a receiver is covered, check the sequence for the corresponding donor
				for donor_obj in sequence_nonCFds_pool:
					if donor_obj.nest_tag == 'donor' and \
							(donor_obj.seq_text[0] != seq_obj.seq_text[0]) and \
							(donor_obj.seq_text[1:] == seq_obj.seq_text[1:]) and \
							(donor_obj.seq_text not in coverage):
						coverage.insert(-2, donor_obj)
						break
			break
	else:
		# CFds case
		search_stop = max(len(max(coverage_pool, key=len)), len(seq))
		for index in range(search_stop, 3, -2):
			for seq_obj in sequence_CFds_pool:
				if (seq[-index:] == seq_obj.seq_text) and (seq_obj.seq_text not in coverage):
					coverage.append(seq_obj)
					break
			else:
				continue
			break

	vertex_dict = {'index': 0, 'coverage': coverage}

	if len(coverage) > 0:
		return CoverageVertex(vertex_dict)
	else:
		return NOT_HIT


def recursive_search(vertex_pool, derived_pool, edge_pool, root, seq, sequence_pool, vertex_count, search_step, max_step):
	seq_branch = [seq + 'r' + seq[-1], seq + 'w1', seq + 'w0']

	for branch in seq_branch:
		vertex = hit_vertex(branch, sequence_pool)
		if isinstance(vertex, CoverageVertex):
			# TODO: fix the problem that the vertex index cannot be updated correctly
			vertex_count += 1
			vertex.index = vertex_count
			vertex.coverage = root.coverage + vertex.coverage

			find_result = find_identical_objs(vertex, vertex_pool, {'index'})
			edge_dict = {'weight': search_step, 'operation_map': [], 'terminal': []}
			if not isinstance(find_result, type(vertex)):
				derived_pool.append(vertex)
				edge_dict['terminal'].extend([str(root.index), str(vertex.index)])
			else:
				edge_dict['terminal'].extend([str(root.index), str(find_result.index)])

			edge_dict['operation_map'].append(branch.replace(root.get_initial_sequence(), '', 1))
			edge = CoverageEdge(edge_dict)
			edge_pool.append(edge)

			return HIT
		elif search_step > max_step:
			return NOT_HIT
		else:
			search_step += 1
			search_result = recursive_search(vertex_pool, derived_pool, edge_pool, root, branch, sequence_pool, vertex_count,
											 search_step, max_step)
			if isinstance(search_result, dict):
				continue

	return


def build_coverage_graph(sequence_pool, init):
	# 递归地向当前顶点的sequence中添加单个操作，存在3个分支。每个分支在遇到一个新的顶点后停止。
	# 设置一个操作计数器，当计数器超过fault list中存在的故障的最大操作数的2倍时停止
	root_seq = Sequence({'seq_text': init, 'ass_init': '', 'detect_tag': '', 'nest_tag': ''})
	root_dict = {'index': 0, 'coverage': [root_seq]}
	root_vertex = CoverageVertex(root_dict)
	vertex_pool = [root_vertex]
	derived_vertices = []
	edge_pool = []
	max_search_step = len(max(map(lambda seq: seq.seq_text, sequence_pool), key=len))
	max_coverage = len(root_vertex.coverage)
	recursive_search(vertex_pool, derived_vertices, edge_pool, root_vertex, root_vertex.get_initial_sequence(),
					 sequence_pool, len(vertex_pool), 0, max_search_step)
	# TODO: need to traverse the derived vertices recursively. note that the initial sequence for each vertex can be
	#  determined. the last sequence of every possible path is the same, and we can use it (plus a detect operation
	#  maybe) as the initial point of new vertices. Since normally, the coverage can be different if the last N
	#  operations of two paths are different, where N is the maximum number of all sequences' operations. While if the
	#  length of same last sequence is smaller than N, the sequence should has been filtered according to the inclusive
	#  principle. As a result, we can only consider the last sequence.
	while max_coverage < len(sequence_pool):
		for derived_vertex in derived_vertices.copy():
			pass
		vertex_pool += derived_vertices
		derived_vertices.clear()
		recursive_search(vertex_pool, derived_vertices, edge_pool, root_vertex, root_vertex.get_initial_sequence(),
						 sequence_pool, len(vertex_pool), 0, max_search_step)
		max_coverage = len(max(map(lambda vertex: vertex.coverage, vertex_pool), key=len))

	pass


if __name__ == '__main__':
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'])
	flat_SF_pool = flatten_sf_pool(filtered_SF_pool)
	filtered_unlinked_pool = (filter_redundant_2cF(filtered_SF_pool, classified_pool['2cF_nonCFds_included'],
												   classified_pool['2cF_CFds']))
	seq_pool = create_sequence_pool(flat_SF_pool, filtered_unlinked_pool, classified_pool['2cF_CFds']['linked'])
	linked_intersection = get_linked_sequence_intersection(seq_pool['linked']['Init_0'], seq_pool['linked']['Init_1'])

	linked_intersection = nest_available_sequences(linked_intersection)
	linked_remainder_0 = nest_available_sequences(seq_pool['linked']['Init_0'] - linked_intersection)
	linked_remainder_1 = nest_available_sequences(seq_pool['linked']['Init_1'] - linked_intersection)

	terminal_seq = find_terminal_seq(linked_intersection, {'Init_0', 'Init_1'})

	build_coverage_graph(linked_intersection, '0')
	print(terminal_seq)

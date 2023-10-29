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

	def get_initial_sequence(self):
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


def filter_vertices(vertex_pool):
	nest_vertices = set(filter(lambda v: len(v.coverage) > 1, vertex_pool))
	unnest_vertices = vertex_pool - nest_vertices
	nest_texts = list(map(lambda v: v.get_initial_sequence(), nest_vertices))
	redundant_vertices = set()
	vertex_count = 0
	for v_obj in unnest_vertices:
		for n_text in nest_texts:
			if v_obj.get_initial_sequence() in n_text:
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

		vertex_dict = {'index': -1, 'coverage': [copy.deepcopy(seq)]}
		vertex_pool.add(CoverageVertex(vertex_dict))

	return filter_vertices(vertex_pool)


def arbit_sequence_nest_condition(seq_segment, seq_obj):
	# check whether the new sequence forms a nest sensitization segment
	max_len = len(seq_obj.seq_text) - 1
	if len(seq_segment) > 2 * max_len:
		if (seq_segment[0] != seq_segment[-1]) and (seq_segment[1:max_len + 1] == seq_segment[max_len + 1:]):
			return NEST
	else:
		return NOT_NEST


def hit_vertex(seq, vertex_pool, max_step):
	# define edges and their weights
	coverage = []
	search_stop = min((max_step - 1) * 2 + 1, len(seq[2:]))

	for index in range(-search_stop, -3, 2):
		for vertex_obj in vertex_pool:
			if seq[index:] == vertex_obj.get_initial_sequence():
				coverage.append(vertex_obj)
				break
		else:
			continue
		break

	if len(coverage) > 0:
		return coverage[0]
	else:
		return NOT_HIT


def recursive_search(root, vertex_seq, vertex_pool, search_step, max_search_step):
	seq_branch = [vertex_seq + 'r' + vertex_seq[-1], vertex_seq + 'w1', vertex_seq + 'w0']
	edge_candidates = set()
	for branch in seq_branch:
		terminal_vertex = hit_vertex(branch, vertex_pool, max_search_step)
		if isinstance(terminal_vertex, CoverageVertex):
			edge_dict = {'weight': search_step, 'operation_map': [], 'terminal': []}
			edge_dict['terminal'].append(terminal_vertex)
			edge_dict['operation_map'].append(branch.replace(root.get_initial_sequence(), '', 1))
			edge = CoverageEdge(edge_dict)
			edge_candidates.add(edge)
			continue
		elif search_step >= max_search_step:
			continue
		else:
			search_result = recursive_search(root, branch, vertex_pool, search_step + 1, max_search_step)
			edge_candidates.update(search_result)

	return edge_candidates


def get_vertex_mirror_text(vertex_obj):
	middle = ''
	mirror_text = vertex_obj.get_initial_sequence()
	if mirror_text.startswith('0'):
		middle = mirror_text.replace('0', '-')
		middle = middle.replace('1', '0')
		middle = middle.replace('-', '1')
	elif mirror_text.startswith('1'):
		middle = mirror_text.replace('1', '-')
		middle = middle.replace('0', '1')
		middle = middle.replace('-', '0')

	mirror_text = middle
	return mirror_text


def split_vertex_pool(vertex_pool):
	split_0 = set()
	split_1 = set()
	for v_obj in vertex_pool:
		if v_obj.get_initial_sequence().startswith('0'):
			split_0.add(v_obj)
		else:
			split_1.add(v_obj)

	return [split_0, split_1]


def define_edges(vertex_pool, coverage_pool):
	edge_pool = []
	max_operation_num = (len(max(coverage_pool, key=len)) - 1) / 2
	transition_num = 1
	# consider the longest sequence situation: a state-transition operation + nest case
	max_search_step = transition_num + int(max_operation_num) * 2 + 1

	for v_obj in vertex_pool:
		# loops are not permitted in graph
		vertex_temp = vertex_pool - {v_obj}
		edge_temp = recursive_search(v_obj, v_obj.get_initial_sequence(), vertex_temp, 1, max_search_step)
		edge_terminals = set(map(lambda e: e.terminal[0], edge_temp))
		for terminal in edge_terminals:
			# choose the edge with minimum weight in all edge candidates between two vertices
			edge_candidate = set(filter(lambda e: terminal in e.terminal, edge_temp))
			edge_winner = min(edge_candidate, key=lambda e: e.weight)
			edge_winner.terminal.insert(0, v_obj)
			edge_pool.append(edge_winner)

	return edge_pool


def build_coverage_graph(sequence_pool):
	graph_pool = []
	vertex_pool = define_vertices(sequence_pool)
	# v_split_pool = split_vertex_pool(vertex_pool)
	v_split_pool = [vertex_pool]
	coverage_pool = set(map(lambda s: s.seq_text, sequence_pool))

	for v_split in v_split_pool:
		vertex_count = 0
		for v_obj in v_split:
			v_obj.index = vertex_count
			vertex_count += 1
		e_split = define_edges(v_split, coverage_pool)
		graph_pool.append([v_split, e_split])

	return graph_pool


def output_graph_for_hp(graph_info, graph_index):
	vertex_info = sorted(graph_info[0], key=lambda v: v.index)
	edge_info = sorted(graph_info[1], key=lambda e: (e.terminal[0].index, e.terminal[1].index))
	vertex_num = len(vertex_info)
	vertex_name = []
	for v_obj in vertex_info:
		vertex_name.append(v_obj.get_initial_sequence())
	graph_file = "../results/sequence_graph_" + str(graph_index) + ".m2c"
	with open(graph_file, 'w') as g:
		g.write(str(vertex_num) + '\n')
		for v_name in vertex_name:
			g.write(v_name + '\n')
		for edge in edge_info:
			head = str(edge.terminal[0].index)
			tail = str(edge.terminal[1].index)
			weight = str(edge.weight)
			g.write(head + ' ' + tail + ' ' + weight + '\n')

	return


if __name__ == '__main__':
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_2cF_pool = (filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']))
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], filtered_2cF_pool)
	flat_SF_pool = flatten_sf_pool(filtered_SF_pool)

	seq_pool = create_sequence_pool(flat_SF_pool, filtered_2cF_pool, classified_pool['2cF_CFds']['linked'])
	linked_intersection = get_linked_sequence_intersection(seq_pool['linked']['Init_0'], seq_pool['linked']['Init_1'])
	for vertex in define_vertices(linked_intersection):
		print(vertex.get_initial_sequence())

	linked_intersection = nest_available_sequences(linked_intersection)
	linked_remainder_0 = nest_available_sequences(seq_pool['linked']['Init_0'] - linked_intersection)
	linked_remainder_1 = nest_available_sequences(seq_pool['linked']['Init_1'] - linked_intersection)

	terminal_seq = set(find_terminal_seq(linked_intersection, {'Init_0', 'Init_1'}))
	graph = build_coverage_graph(linked_intersection)
	for index, subgraph in enumerate(graph):
		output_graph_for_hp(subgraph, index)

	print(terminal_seq)

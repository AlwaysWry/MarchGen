# a module of build the final March test, using the sub-sequences generated in subseq_creator.py
from subseq_creator import *
import itertools as it

HIT = True


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
		self.index = 0
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
	# find terminal sequence in the sequence intersection set.
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


def hit_vertex(seq, text_pool, max_count):
	vertex = CoverageVertex({})
	# TODO: finish the vertex hit principles
	if seq[-2] == 'r':
		# non-CFds case
		search_stop = min(len(seq) - 2, 2 * max_count + 1)
		for index in range(search_stop, 1, -1):
			if seq[-index - 2:-2] in text_pool:
				return seq[-index - 2:-2]
		pass
	else:
		# CFds case
		pass
	pass


def recursive_search(seq, text_pool, max_count):
	seq_branch = []
	seq_branch[0] = seq + 'r' + seq[-1]
	seq_branch[1] = seq + 'w1'
	seq_branch[2] = seq + 'w0'

	for branch in seq_branch:
		vertex = hit_vertex(branch, text_pool, max_count)
		if isinstance(vertex, CoverageVertex):
			pass
			# TODO: vertex should be added into the graph if hit, also build the edge
			return HIT


def generate_edge_weight(root, max_count):

	pass


def build_coverage_graph(sequence_pool):
	# 递归地向当前顶点的sequence中添加单个操作，存在3个分支。每个分支在遇到一个新的顶点后停止。
	# 设置一个操作计数器，当计数器超过fault list中存在的故障的最大操作数的2倍时停止。
	text_pool = set(map(lambda seq: seq.seq_text, sequence_pool))
	max_count = (len(max(text_pool, key=len)) - 1) / 2
	pass


if __name__ == '__main__':
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'])
	flatten_SF_pool = flatten_sf_pool(filtered_SF_pool)
	filtered_unlinked_pool = (filter_redundant_2cF(filtered_SF_pool, classified_pool['2cF_nonCFds_included'],
												   classified_pool['2cF_CFds']))
	seq_pool = create_sequence_pool(flatten_SF_pool, filtered_unlinked_pool, classified_pool['2cF_CFds']['linked'])
	intersection = get_linked_sequence_intersection(seq_pool[0]['Init_0'], seq_pool[0]['Init_1'])
	terminal_seq = find_terminal_seq(intersection, {'Init_0', 'Init_1'})
	build_coverage_graph(intersection)
	print(terminal_seq)


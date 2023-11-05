# a module of build the final March test, using the sub-sequences generated in subseq_creator.py
from subseq_creator import *

NOT_FOUND = False


class CoverageVertex:
	"""data structure of vertices in coverage graph"""

	def __init__(self, prop_dict):
		"""
		:param prop_dict:
		coverage: a list contains all sequences that this vertex covers
		diff: number of different operations by comparing with target operation sequence
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

class MarchElement:
	content = ''
	initial_state = ''
	final_state = ''
	head_tag = False
	tail_tag = False
	ass_tag = False

	def __init__(self, element_text):
		self.content = element_text
		self.initial_state = element_text[1]
		self.final_state = element_text[-1]

	def update_states(self):
		self.initial_state = self.content[1]
		self.final_state = self.content[-1]

# end of element definition


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


def define_vertices(sequence_pool: set):
	vertex_pool = set()
	for seq in sequence_pool:
		vertex_dict = {'coverage': [copy.deepcopy(seq)], 'diff': -1}
		vertex_pool.add(CoverageVertex(vertex_dict))

	return vertex_pool


def find_nest_match(nest_vertex, vertex_pool):
	# check if the corresponding nest sequence exists
	seq = nest_vertex.coverage[0]
	nest_init = lambda i: '0' if i == '1' else '1'
	nest_text = nest_init(seq.seq_text[0]) + seq.seq_text[1:]
	nest_match = list(filter(lambda v: nest_text in map(lambda s: s.seq_text, v.coverage), vertex_pool))
	if len(nest_match) > 0:
		return nest_match[0]
	else:
		return NOT_FOUND


def filter_vertices_covered_by_nest(nest_vertex, vertex_pool):
	nest_text = nest_vertex.get_march_sequence()
	redundant_vertices = set()
	for v_obj in vertex_pool:
		if v_obj.get_march_sequence() in nest_text:
			redundant_vertices.add(v_obj)

	return redundant_vertices


def calculate_diff_value(chain, vertex):
	match_target = vertex.get_march_sequence()
	# each operation appended to chain may cover a new sequence, search from appending 1 operation, and 2 and 3...
	# find how many operations needed at least to cover the target sequence
	search_range = min(len(chain[2:]), len(match_target[:-2]))
	for location in range(-search_range, 1, 2):
		if chain[location:] == match_target[:-location]:
			diff = len(match_target[-location:])
			return diff
	return len(match_target)


def get_vertex_winner(vertex_candidates: set):
	unnest_pool = set(filter(lambda v: v.coverage[0].nest_tag == 'invalid', vertex_candidates))
	nest_pool = vertex_candidates - unnest_pool

	def check_transition(vertex):
		seq = vertex.get_march_sequence()
		return seq[0] == seq[-1]

	# priority level 1: unnest sequences get higher priority
	if len(unnest_pool) > 1:
		# priority level 2 for unnest_pool is bypassed
		secondary_pool = unnest_pool
		# priority level 3: no-transition sequences get higher priority
		if len(secondary_pool) > 1:
			tertiary_pool = set(filter(check_transition, secondary_pool))
		elif len(secondary_pool) == 0:
			secondary_pool = unnest_pool - secondary_pool
			tertiary_pool = set(filter(check_transition, secondary_pool))
		else:
			return next(iter(secondary_pool))
	# priority level 2-2: donor nest sequences get higher priority
	elif len(unnest_pool) == 0:
		secondary_pool = set(filter(lambda v: v.coverage[0].nest_tag == 'donor', nest_pool))
		# priority level 3 under 2-2: no-transition sequences get higher priority
		if len(secondary_pool) > 1:
			tertiary_pool = set(filter(check_transition, secondary_pool))
		elif len(secondary_pool) == 0:
			secondary_pool = nest_pool - secondary_pool
			tertiary_pool = set(filter(check_transition, secondary_pool))
		else:
			return next(iter(secondary_pool))
	else:
		return next(iter(unnest_pool))

	if len(tertiary_pool) > 0:
		return next(iter(tertiary_pool))
	else:
		return next(iter(secondary_pool - tertiary_pool))


def build_coverage_chain(chain: str, vertex_pool: set):
	# a set records the covered vertices by this build process
	covered_vertices = set()
	for v_obj in vertex_pool:
		v_obj.diff = calculate_diff_value(chain, v_obj)

	min_diff = min(map(lambda v: v.diff, vertex_pool))
	vertex_candidates = set(filter(lambda v: v.diff == min_diff, vertex_pool))

	# if the only closest vertex exists, choose it directly, otherwise use priority function
	if len(vertex_candidates) > 1:
		vertex_winner = get_vertex_winner(vertex_candidates)
	else:
		vertex_winner = next(iter(vertex_candidates))

	# the winner is apparently covered, record it
	covered_vertices.add(vertex_winner)
	vertex_for_chain = copy.deepcopy(vertex_winner)

	if vertex_winner.coverage[0].nest_tag == 'donor':
		receiver = find_nest_match(vertex_winner, vertex_pool)
		if isinstance(receiver, CoverageVertex):
			# if the receiver also exists in the pool, it is covered too
			covered_vertices.add(receiver)
			# merge the donor and the receiver vertices
			vertex_for_chain.coverage.extend(copy.deepcopy(receiver.coverage))
			# the corresponding march sequence of vertex gets longer after merging
			vertex_for_chain.diff += (len(receiver.get_march_sequence()) - 3)
			# check if there are other vertices covered by the current nest chain segment
			covered_vertices.update(filter_vertices_covered_by_nest(vertex_for_chain, vertex_pool))

	chain_segment = vertex_for_chain.get_march_sequence()
	if vertex_for_chain.diff < len(chain_segment):
		chain_appendix = chain_segment[-vertex_for_chain.diff:]
	else:
		chain_appendix = 'w' + chain_segment

	return covered_vertices, chain_appendix


def terminal_decorator(chain: str):
	main_elements = []
	terminal_feature = chain[0] + chain[-1]
	match terminal_feature:
		case '00':
			main_elements.append(MarchElement('r1w0' + chain[1:]))
			main_elements[-1].head_tag = True
			main_elements.append(MarchElement(chain[1:] + 'w1'))
			main_elements[-1].tail_tag = True
		case '01':
			main_elements.append(MarchElement('r1w0' + chain[1:] + 'w0'))
			main_elements[-1].head_tag = True
			main_elements[-1].tail_tag = True
			main_elements.append(MarchElement(chain[1:]))
		case '10':
			main_elements.append(MarchElement('r0w1' + chain[1:] + 'w1'))
			main_elements[-1].head_tag = True
			main_elements[-1].tail_tag = True
			main_elements.append(MarchElement(chain[1:]))
		case '11':
			main_elements.append(MarchElement('r0w1' + chain[1:]))
			main_elements[-1].head_tag = True
			main_elements.append(MarchElement(chain[1:] + 'w0'))
			main_elements[-1].tail_tag = True
		case _:
			pass

	if chain[1] != 'r':
		main_elements[-1].content = 'r' + chain[0] + main_elements[-1].content
		main_elements[-1].update_states()

	return main_elements


def construct_main_elements(vertex_pool: set):
	vertex_candidate_pool = copy.deepcopy(vertex_pool)
	initial_vertex = get_vertex_winner(vertex_candidate_pool)
	vertex_candidate_pool -= {initial_vertex}
	coverage_chain = initial_vertex.get_march_sequence()

	while len(vertex_candidate_pool) > 0:
		build_result = build_coverage_chain(coverage_chain, vertex_candidate_pool)
		vertex_candidate_pool -= build_result[0]
		coverage_chain += build_result[1]

	return terminal_decorator(coverage_chain)


def check_odd_sensitization(elements: list, sequence_pool: set):
	# find all sequences that are sensitized even-number times in ME, which violate the odd-sensitization condition
	element_texts = list(map(lambda e: e.content, elements))
	seq_texts = list(map(lambda s: s.seq_text, sequence_pool))
	violated_pool = set()

	for index, text in enumerate(element_texts):
		count_dict = dict.fromkeys(seq_texts, 0)

		for seq_obj in sequence_pool:
			search_range = len(seq_obj.seq_text)
			search_scope = text[1] + text
			if index > 0 and search_range > 3:
				search_scope = element_texts[index - 1][-search_range + 2:-1] + search_scope

			for location in range(0, len(search_scope[:-search_range + 1]), 2):
				segment = search_scope[location:location + search_range]
				if segment == seq_obj.seq_text:
					count_dict[seq_obj.seq_text] += 1

		even_pool = set(filter(lambda k: count_dict[k] % 2 == 0, count_dict.keys()))

		if index > 0:
			violated_pool = violated_pool.intersection(even_pool)
		else:
			violated_pool.update(even_pool)

	if len(violated_pool) > 0:
		return violated_pool
	else:
		return NOT_FOUND


def check_tail_cover(elements: list, sequence_pool: set):
	# check whether the ME adds tail terminal introduces sequence that violated order condition of 2cF2aa
	tail_cover = set()
	for element in elements:
		element_text = element.content
		search_range = sorted(set(map(lambda s: len(s.seq_text), sequence_pool)), reverse=True)
		for location in search_range:
			order_violation = set(filter(lambda s: s.seq_text == element_text[-location:], sequence_pool))
			if len(order_violation) > 0:
				tail_cover.update(order_violation)
				break

	if len(tail_cover) > 0:
		return tail_cover
	else:
		return NOT_FOUND


def construct_ass_elements():
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

	vertex_intersection = define_vertices(linked_intersection)
	main_elements = construct_main_elements(vertex_intersection)
	print(main_elements)
	print(check_odd_sensitization(main_elements, linked_intersection))
	print(check_tail_cover(list(filter(lambda e: e.tail_tag, main_elements)), linked_intersection))

	terminal_seq = set(find_terminal_seq(linked_intersection, {'Init_0', 'Init_1'}))

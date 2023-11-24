# a module of build the final March test, using the sub-sequences generated in subseq_creator.py
from subseq_creator import *

NOT_FOUND = False
NO_ELEMENT = False


class CoverageVertex:
	"""data structure of vertices in coverage graph"""
	coverage = []
	diff = -1

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
		elif len(self.__dict__['coverage']) == 1:
			init_seq = self.__dict__['coverage'][0].seq_text
		else:
			init_seq = ''
			return init_seq

		if self.__dict__['coverage'][-1].detect_tag:
			init_seq += 'r' + init_seq[-1]

		return init_seq


# end of vertex definition


class MarchElement:
	"""
	content: element text
	initial_state: the state at the beginning of element
	final_state: the state after the element is executed
	head_tag: tag for record if the element is decorated at its head
	tail_tag: tag for record if the element is decorated at its tail
	ass_tag: tag for whether the element is an associate element
	"""
	content = ''
	tied_element = ''
	initial_state = ''
	final_state = ''
	address_order = ''
	head_tag = False
	tail_tag = False
	ass_tag = False

	def __init__(self, element_text):
		self.content = element_text
		if len(element_text) > 0:
			self.initial_state = element_text[1]
			self.final_state = element_text[-1]

	def update_states(self):
		self.initial_state = self.content[1]
		self.final_state = self.content[-1]


# end of element definition


class LinkedElementsBuilder:
	@staticmethod
	def get_vertex_winner(vertex_candidates: set, aux_pool: set, last_winner: CoverageVertex, init: str):
		unnest_pool = set(filter(lambda v: v.coverage[0].nest_tag == 'invalid', vertex_candidates))
		nest_pool = vertex_candidates - unnest_pool

		def check_transition(vertex):
			seq = vertex.get_march_sequence()
			return seq[0] == seq[-1]

		# priority level 1: unnest sequences get higher priority
		if len(unnest_pool) > 1:
			# priority level 2-1 (for unnest_pool) is bypassed
			secondary_pool = unnest_pool
			# priority level 3: no-transition sequences get higher priority
			if len(secondary_pool) > 1:
				tertiary_pool = set(filter(check_transition, secondary_pool))
			elif len(secondary_pool) == 0:
				secondary_pool = unnest_pool - secondary_pool
				tertiary_pool = set(filter(check_transition, secondary_pool))
			else:
				return next(iter(secondary_pool))
		# priority level 2-2 (for nest pool): donor nest sequences get higher priority
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

	@staticmethod
	def terminal_decorator(chain: str):
		main_elements = {'01_me': MarchElement(''), '10_me': MarchElement('')}
		terminal_feature = chain[0] + chain[-1]
		match terminal_feature:
			case '00':
				main_elements['10_me'] = MarchElement('r1w0' + chain[1:])
				main_elements['10_me'].head_tag = True
				main_elements['01_me'] = MarchElement(chain[1:] + 'w1')
				main_elements['01_me'].tail_tag = True
				if chain[1] != 'r':
					main_elements['01_me'].content = 'r' + chain[0] + main_elements['01_me'].content
					main_elements['01_me'].update_states()
			case '01':
				main_elements['10_me'] = MarchElement('r1w0' + chain[1:] + 'w0')
				main_elements['10_me'].head_tag = True
				main_elements['10_me'].tail_tag = True
				main_elements['01_me'] = MarchElement(chain[1:])
				if chain[1] != 'r':
					main_elements['01_me'].content = 'r' + chain[0] + main_elements['01_me'].content
					main_elements['01_me'].update_states()
			case '10':
				main_elements['01_me'] = MarchElement('r0w1' + chain[1:] + 'w1')
				main_elements['01_me'].head_tag = True
				main_elements['01_me'].tail_tag = True
				main_elements['10_me'] = MarchElement(chain[1:])
				if chain[1] != 'r':
					main_elements['10_me'].content = 'r' + chain[0] + main_elements['10_me'].content
					main_elements['10_me'].update_states()
			case '11':
				main_elements['01_me'] = MarchElement('r0w1' + chain[1:])
				main_elements['01_me'].head_tag = True
				main_elements['10_me'] = MarchElement(chain[1:] + 'w0')
				main_elements['10_me'].tail_tag = True
				if chain[1] != 'r':
					main_elements['10_me'].content = 'r' + chain[0] + main_elements['10_me'].content
					main_elements['10_me'].update_states()
			case _:
				pass

		return main_elements


# end of LinkedElementsBuilder definition


def get_linked_CFds_union(init_0_pool, init_1_pool):
	seq_union = set()
	seq_union.update(init_1_pool)

	for seq in init_0_pool:
		find_result = find_identical_objs(seq, seq_union, {'ass_init'})
		if isinstance(find_result, type(DIFFERENT)):
			seq_union.add(seq)

	for seq in seq_union:
		seq.ass_init = -1

	return seq_union


def define_vertices(sequence_pool: set):
	vertex_pool = set()
	for seq in sequence_pool:
		vertex_dict = {'coverage': [copy.deepcopy(seq)], 'diff': -1}
		vertex_pool.add(CoverageVertex(vertex_dict))

	return vertex_pool


def find_nest_match(nest_vertex, vertex_pool):
	# check if the corresponding nest sequence exists
	seq = nest_vertex.coverage[0]

	def nest_init(seq_text):
		if seq_text == '1':
			return '0'
		else:
			return '1'

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


def build_coverage_chain(chain: str, vertex_pool: set, aux_vertex_pool: set, last_winner: CoverageVertex, init: str, builder: type.__name__):
	# a set records the covered vertices by this build process
	covered_vertices = set()
	for v_obj in vertex_pool:
		v_obj.diff = calculate_diff_value(chain, v_obj)

	min_diff = min(map(lambda v: v.diff, vertex_pool))
	vertex_candidates = set(filter(lambda v: v.diff == min_diff, vertex_pool))

	# if the only closest vertex exists, choose it directly, otherwise use priority function
	if len(vertex_candidates) > 1:
		vertex_winner = builder.get_vertex_winner(vertex_candidates, aux_vertex_pool, last_winner, init)
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
			# the corresponding march sequence of vertex gets longer after merging, calculate diff value again
			vertex_for_chain.diff += (len(receiver.get_march_sequence()) - 3)
			# check if there are other vertices covered by the current nest chain segment
			covered_vertices.update(filter_vertices_covered_by_nest(vertex_for_chain, vertex_pool))

	chain_segment = vertex_for_chain.get_march_sequence()
	if vertex_for_chain.diff < len(chain_segment):
		chain_appendix = chain_segment[-vertex_for_chain.diff:]
	else:
		chain_appendix = 'w' + chain_segment

	return covered_vertices, chain_appendix


def construct_main_elements(vertex_pool: set):
	vertex_candidate_pool = copy.deepcopy(vertex_pool)
	initial_vertex = LinkedElementsBuilder.get_vertex_winner(vertex_candidate_pool, set(), CoverageVertex({'coverage': [], 'diff': -1}), 'Init_-1')
	vertex_candidate_pool -= {initial_vertex}
	coverage_chain = initial_vertex.get_march_sequence()

	while len(vertex_candidate_pool) > 0:
		build_result = build_coverage_chain(coverage_chain, vertex_candidate_pool, set(), initial_vertex, 'Init_-1', LinkedElementsBuilder)
		vertex_candidate_pool -= build_result[0]
		coverage_chain += build_result[1]

	return LinkedElementsBuilder.terminal_decorator(coverage_chain)


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


def construct_tail_cover_elements(tail_cover):
	# for each sequence object in tail_cover, it should be regarded as a CFds, since it just violates the
	# no-swap condition for CFds
	tail_cover_pool = set()
	if len(tail_cover) > 0:
		for seq_obj in tail_cover:
			tail_cover_seq = Sequence(seq_obj.__dict__.copy())
			tail_cover_seq.detect_tag = False
			tail_cover_vertex = CoverageVertex({'diff': -1, 'coverage': [tail_cover_seq]})
			tail_cover_pool.add(tail_cover_vertex)

		vertex_candidate_pool = copy.deepcopy(tail_cover_pool)
		initial_vertex = LinkedElementsBuilder.get_vertex_winner(vertex_candidate_pool, set(), CoverageVertex({'coverage': [], 'diff': -1}), 'Init_-1')
		vertex_candidate_pool -= {initial_vertex}
		coverage_chain = initial_vertex.get_march_sequence()

		while len(vertex_candidate_pool) > 0:
			build_result = build_coverage_chain(coverage_chain, vertex_candidate_pool, set(), initial_vertex, 'Init_-1', LinkedElementsBuilder)
			vertex_candidate_pool -= build_result[0]
			coverage_chain += build_result[1]

		tail_cover_me = coverage_chain[1:]
		if not tail_cover_me.startswith('r'):
			tail_cover_me = 'r' + coverage_chain[0] + tail_cover_me

		return tail_cover_me
	else:
		return NO_ELEMENT


def check_odd_sensitization(elements, sequence_pool: set):
	# find all sequences that are sensitized even-number times in ME, which violate the odd-sensitization condition
	element_text_list = list(map(lambda e: e.content, elements.values()))
	seq_texts = list(map(lambda s: s.seq_text, sequence_pool))
	violation_pool = []

	def get_violated_seq_texts(element_list):
		violated_seq_texts = set()
		for index, text in enumerate(element_list):
			count_dict = dict.fromkeys(seq_texts, 0)
			for seq_obj in sequence_pool:
				search_range = len(seq_obj.seq_text)
				search_scope = text[1] + text
				if (index > 0) and (search_range > 3):
					search_scope = element_list[index - 1][-search_range + 2:-1] + search_scope

				for location in range(0, len(search_scope[:-search_range + 1]), 2):
					segment = search_scope[location:location + search_range]
					if segment == seq_obj.seq_text:
						count_dict[seq_obj.seq_text] += 1

			even_pool = set(filter(lambda k: count_dict[k] % 2 == 0, count_dict.keys()))

			if index > 0:
				violated_seq_texts = violated_seq_texts.intersection(even_pool)
			else:
				violated_seq_texts.update(even_pool)

		return set(filter(lambda s: s.seq_text in violated_seq_texts, sequence_pool))

	# save the violations of odd-sensitization for each order case
	violation_pool.append((get_violated_seq_texts(element_text_list), element_text_list[1], element_text_list[0]))
	element_text_list.reverse()
	violation_pool.append((get_violated_seq_texts(element_text_list), element_text_list[1], element_text_list[0]))

	if max(map(lambda t: len(t[0]), violation_pool)) > 0:
		return violation_pool
	else:
		return NOT_FOUND


def construct_odd_sensitization_elements(odd_violation, main_elements):
	# for each sequence object in odd_violation, it should be regarded as a CFds, since it just violates the
	# odd-sensitization condition for CFds
	violation_pool = set()
	violation_me_candidates = set()

	if isinstance(odd_violation, list):
		for odd_case in odd_violation:
			for seq_obj in odd_case[0]:
				violation_seq = Sequence(seq_obj.__dict__.copy())
				violation_seq.detect_tag = False
				violation_vertex = CoverageVertex({'diff': -1, 'coverage': [violation_seq]})
				violation_pool.add(violation_vertex)

			vertex_candidate_pool = copy.deepcopy(violation_pool)
			max_range = max(set(map(lambda s: len(s.seq_text), odd_case[0]))) - 2
			search_range = min(max_range, len(odd_case[1]))
			coverage_chain = odd_case[1][-search_range:] + 'r' + odd_case[1][-1]

			# the initial coverage chain may also cover a sequence
			initial_vertex = set(filter(lambda v: v.coverage[0].seq_text in coverage_chain, vertex_candidate_pool))
			vertex_candidate_pool -= initial_vertex

			while len(vertex_candidate_pool) > 0:
				build_result = build_coverage_chain(coverage_chain, vertex_candidate_pool, set(), CoverageVertex({'coverage': [], 'diff': -1}), 'Init_-1', LinkedElementsBuilder)
				vertex_candidate_pool -= build_result[0]
				coverage_chain += build_result[1]

			# each odd_sensitization candidate has to follow the tied precedent element, or the odd sensitization may be destroyed
			violation_me_candidates.add((coverage_chain[search_range:], odd_case[1], odd_case[2]))

		odd_sensitization = min(violation_me_candidates, key=lambda c: len(c[0]))
		return odd_sensitization
	else:
		return (NO_ELEMENT,)


def construct_ass_elements(main_elements, sequence_pool):
	ass_elements = {'tail_cover_me': MarchElement(''), 'odd_sensitization_me': MarchElement('')}
	# check and build the ME for tail-cover first
	tail_cover = check_tail_cover(list(filter(lambda m: m.tail_tag, main_elements.values())), sequence_pool)
	tail_cover_me_text = construct_tail_cover_elements(tail_cover)

	if isinstance(tail_cover_me_text, str):
		tail_cover_me = MarchElement(tail_cover_me_text)
		tail_cover_me.ass_tag = True
		ass_elements['tail_cover_me'] = tail_cover_me

	# check and build the ME for odd sensitization violation
	odd_violation = check_odd_sensitization(main_elements, sequence_pool)
	construct_result = construct_odd_sensitization_elements(odd_violation, main_elements)
	odd_sensitization_me_text = construct_result[0]
	if isinstance(odd_sensitization_me_text, str):
		odd_sensitization_me = MarchElement(odd_sensitization_me_text)
		odd_sensitization_me.ass_tag = True
		odd_sensitization_me.tied_element = []
		for tied_element_text in construct_result[1:]:
			odd_sensitization_me.tied_element.extend(list(filter(lambda m: m.content == tied_element_text, main_elements.values())))
		ass_elements['odd_sensitization_me'] = odd_sensitization_me

	return ass_elements


def linked_CFds_constructor(linked_pool):
	# For 01/10 type main ME, all linked seq need to be added for constructing; for 00/11 type, init_0 and init_1 pool
	# should be added for constructing separately. however, the different seq set in 00 and 11 ME may introduce the
	# violation to sensitization order. As a result, we only use the 01/10 type ME, which means that all sequences in
	# linked pool should be contained in both 2 ME, and the ass_init of them can be ignored
	union_pool = get_linked_CFds_union(linked_pool['Init_0'], linked_pool['Init_1'])
	vertex_union = define_vertices(union_pool)
	main_mes = construct_main_elements(vertex_union)
	ass_mes = construct_ass_elements(main_mes, union_pool)

	return {'main_me': main_mes, 'ass_me': ass_mes}


if __name__ == '__main__':
	os.chdir("../")
	sys.path.append("src/")
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_2cF_pool = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']['unlinked'])
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], filtered_2cF_pool)
	flat_SF_pool = flatten_sf_pool(filtered_SF_pool)

	seq_pool = create_sequence_pool(flat_SF_pool, filtered_2cF_pool, classified_pool['2cF_CFds']['linked'])

	if len(seq_pool['linked']['Init_0']) + len(seq_pool['linked']['Init_1']) > 0:
		linked_CFds_constructor(seq_pool['linked'])

	pass

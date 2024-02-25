# a module of build the final March test, using the sub-sequences generated in subseq_creator.py
from subseq_creator import *

NOT_FOUND = False
NO_ELEMENT = False
NOT_ALLOWED = False


class CoverageVertex:
	"""data structure of vertices in coverage graph
	coverage: a list contains all sequences that this vertex covers
	diff: number of different operations by comparing with target operation sequence
	"""

	def __init__(self, prop_dict):
		self.coverage = []
		self.diff = -1
		self.__dict__.update(prop_dict)

	def get_march_segment(self):
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

	def __init__(self, element_text):
		self.content = ''
		self.tied_element = ''
		self.initial_state = ''
		self.final_state = ''
		self.address_order = ''
		self.head_tag = False
		self.tail_tag = False
		self.ass_tag = False
		self.transition_tag = False

		self.content = element_text
		if len(element_text) > 0:
			self.initial_state = element_text[1]
			self.final_state = element_text[-1]

	def update_states(self):
		self.initial_state = self.content[1]
		self.final_state = self.content[-1]


# end of element definition


class LinkedMainElementsBuilder:
	@staticmethod
	def get_vertex_winner(vertex_candidates: set, aux_pool=None, last_winner=None, init='Init_-1'):
		donor_pool = set(filter(lambda v: v.coverage[0].nest_tag == 'donor', vertex_candidates))

		def check_transition(vertex):
			seq = vertex.get_march_segment()
			return seq[0] == seq[-1]

		# priority level 1: donor nest sequences get higher priority
		if len(donor_pool) == 1:
			return next(iter(donor_pool))
		elif len(donor_pool) > 1:
			# priority level 2: no-transition sequences get higher priority
			secondary_pool = set(filter(check_transition, donor_pool))
			if len(secondary_pool) > 0:
				return next(iter(secondary_pool))
			else:
				return next(iter(donor_pool - secondary_pool))
		else:
			secondary_pool = set(filter(check_transition, vertex_candidates))
			if len(secondary_pool) > 0:
				return next(iter(secondary_pool))
			else:
				return next(iter(vertex_candidates - secondary_pool))

	@staticmethod
	def terminal_decorator(chain: str):
		main_elements = {'01_me': MarchElement(''), '10_me': MarchElement('')}
		terminal_feature = chain[0] + chain[-1]
		match terminal_feature:
			case '00':
				main_elements['10_me'] = MarchElement('r1w' + chain)
				main_elements['10_me'].head_tag = True
				main_elements['01_me'] = MarchElement(chain[1:] + 'w1')
				main_elements['01_me'].tail_tag = True
				if chain[1] != 'r':
					main_elements['01_me'].content = 'r' + chain[0] + main_elements['01_me'].content
					main_elements['01_me'].update_states()
			case '01':
				main_elements['10_me'] = MarchElement('r1w' + chain + 'w0')
				main_elements['10_me'].head_tag = True
				main_elements['10_me'].tail_tag = True
				main_elements['01_me'] = MarchElement(chain[1:])
				if chain[1] != 'r':
					main_elements['01_me'].content = 'r' + chain[0] + main_elements['01_me'].content
					main_elements['01_me'].update_states()
			case '10':
				main_elements['01_me'] = MarchElement('r0w' + chain + 'w1')
				main_elements['01_me'].head_tag = True
				main_elements['01_me'].tail_tag = True
				main_elements['10_me'] = MarchElement(chain[1:])
				if chain[1] != 'r':
					main_elements['10_me'].content = 'r' + chain[0] + main_elements['10_me'].content
					main_elements['10_me'].update_states()
			case '11':
				main_elements['01_me'] = MarchElement('r0w' + chain)
				main_elements['01_me'].head_tag = True
				main_elements['10_me'] = MarchElement(chain[1:] + 'w0')
				main_elements['10_me'].tail_tag = True
				if chain[1] != 'r':
					main_elements['10_me'].content = 'r' + chain[0] + main_elements['10_me'].content
					main_elements['10_me'].update_states()
			case _:
				pass

		return main_elements


# end of LinkedMainElementsBuilder definition


def get_linked_CFds_union(init_0_pool, init_1_pool):
	seq_union = set()
	seq_union.update(init_1_pool)

	for seq in init_0_pool:
		find_result = find_identical_objs(seq, seq_union, {'detect_tag', 'dr_tag', 'ass_init', 'nest_tag'})
		if isinstance(find_result, type(DIFFERENT)):
			seq_union.add(seq)
		elif seq.detect_tag:
			find_result.detect_tag |= seq.detect_tag
			find_result.dr_tag |= seq.dr_tag
			find_result.nest_tag = seq.nest_tag

	for seq in seq_union:
		seq.ass_init = -1

	return seq_union


def define_vertices(sequence_pool: set):
	vertex_pool = set()
	for seq in sequence_pool:
		vertex_dict = {'coverage': [seq], 'diff': -1}
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


def check_vertices_covered_by_nest(nest_vertex, vertex_pool):
	nest_text = nest_vertex.get_march_segment()
	nest_seq_text = set(map(lambda s: s.seq_text, nest_vertex.coverage))
	redundant_vertices = set()
	for v_obj in vertex_pool:
		# sensitized and also detected, seen as covered
		if v_obj.get_march_segment() in nest_text:
			redundant_vertices.add(v_obj)
		# if just sensitized but not detected, the current nest sequence is not allowed
		elif (v_obj.coverage[0].seq_text in nest_text[:-2]) and (v_obj.coverage[0].seq_text not in nest_seq_text):
			return NOT_ALLOWED

	return redundant_vertices


def calculate_diff_value(chain, vertex):
	match_target = vertex.get_march_segment()
	# each operation appended to chain may cover a new sequence, search from appending 1 operation, and 2 and 3...
	# find how many operations needed at least to cover the target sequence
	search_range = min(len(chain), len(match_target[:-2]))
	for location in range(-search_range, 1, 2):
		if chain[location:] == match_target[:-location]:
			diff = len(match_target[-location:])
			return diff
	return len(match_target)


def build_coverage_chain(chain, seq_check_range, vertex_pool, builder, aux_vertex_pool=None, last_winner=None, init='Init_-1'):
	# a set records the covered vertices by this build process
	covered_vertices = []
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
	covered_vertices.append(vertex_winner)
	vertex_for_chain = copy.deepcopy(vertex_winner)

	if vertex_winner.coverage[0].nest_tag == 'donor':
		chain_check_range = len(chain) - (len(vertex_winner.get_march_segment()) - vertex_winner.diff) + 1

		# if the operations before the current nest sequence is less than the longest sequence, it cannot be known
		# if mis-sensitization sequences exist, since the former ME is not decided yet, as a result, not allowed nest
		# sequence in this case
		if chain_check_range >= seq_check_range:
			receiver = find_nest_match(vertex_winner, vertex_pool)
			if isinstance(receiver, CoverageVertex):
				# merge the donor and the receiver vertices
				vertex_for_chain.coverage.extend(copy.deepcopy(receiver.coverage))
				# check if there are other vertices covered by the current nest chain segment
				check_result = check_vertices_covered_by_nest(vertex_for_chain, vertex_pool)
				if isinstance(check_result, set):
					# if the receiver also exists in the pool, it is covered too
					covered_vertices.append(receiver)
					# the corresponding march sequence of vertex gets longer after merging, calculate diff value again
					vertex_for_chain.diff += (len(receiver.get_march_segment()) - 3)
					covered_vertices.extend(list(check_result))
				else:
					# if there are sequences just sensitized but not detected, the nest sequence is not allowed
					del vertex_for_chain.coverage[-1]

	chain_segment = vertex_for_chain.get_march_segment()
	if vertex_for_chain.diff < len(chain_segment):
		chain_appendix = chain_segment[-vertex_for_chain.diff:]
	else:
		chain_appendix = 'w' + chain_segment

	return covered_vertices, chain_appendix


def construct_main_elements(vertex_pool: set):
	vertex_candidate_pool = copy.deepcopy(vertex_pool)
	covered_vertex_pool = set()
	initial_vertex = LinkedMainElementsBuilder.get_vertex_winner(vertex_candidate_pool)
	vertex_candidate_pool -= {initial_vertex}
	coverage_chain = initial_vertex.get_march_segment()

	nest_check_range = max(map(lambda v: len(v.coverage[0].seq_text), vertex_pool))

	while len(vertex_candidate_pool) > 0:
		build_result = build_coverage_chain(coverage_chain, nest_check_range, vertex_candidate_pool, LinkedMainElementsBuilder)
		vertex_candidate_pool -= set(build_result[0])
		covered_vertex_pool.update(build_result[0])
		coverage_chain += build_result[1]

	return LinkedMainElementsBuilder.terminal_decorator(coverage_chain), coverage_chain

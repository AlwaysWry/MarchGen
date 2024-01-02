from linked_ass_constructor import *


class UnlinkedElementsBuilder:
	@staticmethod
	def get_vertex_winner(vertex_candidates: set, aux_pool: set, last_winner: CoverageVertex, init: str):

		def check_states(vertex):
			seq = vertex.get_march_sequence()
			if len(seq) > 0:
				return seq[0] + seq[-1]
			else:
				return ''

		aux_pool_dict = {'00': set(), '01': set(), '10': set(), '11': set()}
		for feature in aux_pool_dict.keys():
			aux_feature_pool = set(filter(lambda v: check_states(v) == feature, aux_pool))
			aux_pool_dict[feature].update(aux_feature_pool)

		def priority_filter(upper_pool: set, priority_feature: str):
			priority_pool = set(filter(lambda v: check_states(v) == priority_feature, upper_pool))
			aux_priority_pool = aux_pool_dict[priority_feature]

			if len(priority_pool) == 0 and len(aux_priority_pool) == 0:
				priority_pool = upper_pool - priority_pool
				return priority_pool
			elif len(priority_pool) == 0:
				priority_pool.update(aux_priority_pool)
				return next(iter(priority_pool))
			else:
				return next(iter(priority_pool))

		feature_dict = {'Init_0': ['01', '11', '00', '10'], 'Init_1': ['10', '00', '11', '01']}
		feature_list = feature_dict[init]

		state_history = check_states(last_winner)
		if (state_history == feature_list[0]) or (state_history == feature_list[1]):
			primary_result = priority_filter(vertex_candidates, feature_list[1])
			if isinstance(primary_result, CoverageVertex):
				return primary_result

			secondary_result = priority_filter(primary_result, feature_list[3])
			if isinstance(secondary_result, CoverageVertex):
				return secondary_result

			return next(iter(secondary_result))
		else:
			primary_result = priority_filter(vertex_candidates, feature_list[2])
			if isinstance(primary_result, CoverageVertex):
				return primary_result

			secondary_result = priority_filter(primary_result, feature_list[0])
			if isinstance(secondary_result, CoverageVertex):
				return secondary_result

			return next(iter(secondary_result))

	@staticmethod
	def terminal_decorator(chain: str, init: str):
		main_element = MarchElement('')
		terminal_feature = chain[0] + chain[-1]
		match terminal_feature:
			case '00':
				if init == 'Init_0':
					main_element = MarchElement(chain[1:])
				else:
					main_element = MarchElement('r1w0' + chain[1:] + 'w1')
					main_element.head_tag = True
					main_element.tail_tag = True
			case '01':
				if init == 'Init_0':
					main_element = MarchElement(chain[1:] + 'w0')
					main_element.tail_tag = True
				else:
					main_element = MarchElement('r1w0' + chain[1:])
					main_element.head_tag = True
			case '10':
				if init == 'Init_1':
					main_element = MarchElement(chain[1:] + 'w1')
					main_element.tail_tag = True
				else:
					main_element = MarchElement('r0w1' + chain[1:])
					main_element.head_tag = True
			case '11':
				if init == 'Init_1':
					main_element = MarchElement(chain[1:])
				else:
					main_element = MarchElement('r0w1' + chain[1:] + 'w0')
					main_element.head_tag = True
					main_element.tail_tag = True
			case _:
				pass

		if chain[1] != 'r':
			main_element.content = 'r' + chain[0] + main_element.content
			main_element.update_states()

		return main_element


# end of UnlinkedElementsBuilder definition


def construct_unlinked_element(vertex_pool: set, vertex_aux_pool: set, init: str):
	vertex_candidate_pool = copy.deepcopy(vertex_pool)
	initial_vertex = UnlinkedElementsBuilder.get_vertex_winner(vertex_candidate_pool, vertex_aux_pool, CoverageVertex({'coverage': [], 'diff': -1}), init)
	vertex_candidate_pool -= {initial_vertex}
	coverage_chain = initial_vertex.get_march_sequence()

	while len(vertex_candidate_pool) > 0:
		build_result = build_coverage_chain(coverage_chain, -1, vertex_candidate_pool, UnlinkedElementsBuilder, vertex_aux_pool, initial_vertex, init)

		for vertex in build_result[0]:
			if vertex in vertex_candidate_pool:
				vertex_candidate_pool.remove(vertex)
			elif vertex in vertex_aux_pool:
				vertex_aux_pool.remove(vertex)

		coverage_chain += build_result[1]

	return coverage_chain


def unlinked_2cF_constructor(unlinked_pool):
	# since the sequences in unlinked_pool are unnecessary to abide any 2cF3 or 2cF2aa conditions, use 00/11 ME structure
	# to make sure no redundant sequence will be contained. Besides, the instant-detect rule still need to be obeyed for
	# nonCFds*nonCFds-type 2cFs
	vertex_2cF_pool_0 = define_vertices(unlinked_pool['Init_0'])
	vertex_2cF_pool_1 = define_vertices(unlinked_pool['Init_1'])
	vertex_sf_pool = define_vertices(unlinked_pool['Init_-1'])

	# build "00" ME
	unlinked_me_00_text = construct_unlinked_element(vertex_2cF_pool_0, vertex_sf_pool, 'Init_0')
	unlinked_me_00 = UnlinkedElementsBuilder.terminal_decorator(unlinked_me_00_text, 'Init_0')
	# build "11" ME
	unlinked_me_11_text = construct_unlinked_element(vertex_2cF_pool_1, vertex_sf_pool, 'Init_1')
	unlinked_me_11 = UnlinkedElementsBuilder.terminal_decorator(unlinked_me_11_text, 'Init_1')

	return {'00_me': unlinked_me_00, '11_me': unlinked_me_11}


def sf_constructor(unlinked_pool):
	vertex_sf_pool = define_vertices(unlinked_pool['Init_-1'])
	sf_me_candidates = set()

	def get_sf_me_candidates(init, vertex_candidate_pool):
		initial_vertex = UnlinkedElementsBuilder.get_vertex_winner(vertex_candidate_pool, set(), CoverageVertex({'coverage': [], 'diff': -1}), init)
		vertex_candidate_pool -= {initial_vertex}
		coverage_chain = initial_vertex.get_march_sequence()

		while len(vertex_candidate_pool) > 0:
			build_result = build_coverage_chain(coverage_chain, -1, vertex_candidate_pool, UnlinkedElementsBuilder, set(), initial_vertex, init)
			vertex_candidate_pool -= build_result[0]
			coverage_chain += build_result[1]
		# make sure that the SF me also starts with read operation, since terminal decorator is not applied
		if coverage_chain[1] != 'r':
			coverage_chain = coverage_chain[0] + 'r' + coverage_chain

		return MarchElement(coverage_chain[1:])

	sf_me_candidates.add(get_sf_me_candidates('Init_0', copy.deepcopy(vertex_sf_pool)))
	sf_me_candidates.add(get_sf_me_candidates('Init_1', copy.deepcopy(vertex_sf_pool)))

	return sf_me_candidates


if __name__ == '__main__':
	os.chdir("../")
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_2cF_pool = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']['unlinked'])
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], filtered_2cF_pool)
	flat_SF_pool = flatten_sf_pool(filtered_SF_pool)

	seq_pool = create_sequence_pool(flat_SF_pool, filtered_2cF_pool, classified_pool['2cF_CFds']['linked'])

	if len(seq_pool['unlinked']['Init_0']) + len(seq_pool['unlinked']['Init_1']) > 0:
		unlinked_2cF_constructor(seq_pool['unlinked'])
	elif len(seq_pool['unlinked']['Init_-1']) > 0:
		sf_constructor(seq_pool['unlinked'])
	pass

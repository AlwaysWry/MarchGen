from linked_CFds_ass_constructor import *

OVERLAP = True
NO_OVERLAP = False


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

			return next(iter(secondary_result), CoverageVertex({'coverage': [], 'diff': -1}))

	@staticmethod
	def terminal_decorator(chain: str, init: str):
		main_element = MarchElement('')

		if len(chain) > 0:
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
	initial_vertex = UnlinkedElementsBuilder.get_vertex_winner(vertex_candidate_pool, vertex_aux_pool,
															   CoverageVertex({'coverage': [], 'diff': -1}), init)
	vertex_candidate_pool -= {initial_vertex}
	coverage_chain = initial_vertex.get_march_sequence()

	while len(vertex_candidate_pool) > 0:
		build_result = build_coverage_chain(coverage_chain, -1, vertex_candidate_pool, UnlinkedElementsBuilder,
											vertex_aux_pool, initial_vertex, init)

		for vertex in build_result[0]:
			if vertex in vertex_candidate_pool:
				vertex_candidate_pool.remove(vertex)
			elif vertex in vertex_aux_pool:
				vertex_aux_pool.remove(vertex)

		coverage_chain += build_result[1]

	return coverage_chain


# filter the degenerated 2cFs, SFs and the nonCFds*nonCFds faults that are possibly covered by the linked CFds main MEs.
# Besides the segments appended in the linked main MEs, other faults may be covered because of the nest between two segments
def inter_ME_filter(main_elements, main_middle_part, degenerated_seq_pool, sf_seq_pool, undetermined_faults):
	main_element_contents = {'01_me': '', '10_me': ''}
	for element_key in main_elements.keys():
		main_element_contents[element_key] = main_elements[element_key].content[1:]

	def search_seq_redundancy(candidate_seq, init_key):
		# if the sequence comes from a nonCFds FP, the match candidate text should be seq_text and a read-verify operation,
		# and the candidate text should be included in both ME, since the main MEs are 01/10 type.
		if candidate_seq.detect_tag:
			candidate_text = candidate_seq.seq_text + 'r' + candidate_seq.seq_text[-1]
		else:
			candidate_text = candidate_seq.seq_text

		search_result = filter(lambda m: candidate_text in m, main_element_contents)
		if (init_key == 'Init_-1') and (len(search_result) > 0):
			return REDUNDANT
		elif len(search_result) > 1:
			return REDUNDANT

		return NOT_REDUNDANT

	# check the redundancy according to their corresponding test conditions
	# consider the 2cFs in degenerated pool
	degenerated_redundant_pool = set()
	for degenerated_key in degenerated_seq_pool.keys():
		for degenerated_seq in degenerated_seq_pool[degenerated_key]:
			redundancy = search_seq_redundancy(degenerated_seq, degenerated_key)
			if redundancy:
				degenerated_redundant_pool.add(degenerated_seq)

	# consider SFs
	sf_redundancy_pool = set()
	for sf_key in sf_seq_pool.keys():
		for sf_seq in sf_seq_pool[sf_key]:
			redundancy = search_seq_redundancy(sf_seq, sf_key)
			if redundancy:
				sf_redundancy_pool.add(sf_seq)

	# consider undetermined faults, check whether one of the FP is included in the main ME, and no overlapping
	undetermined_redundancy_pool = set()

	def check_overlapping(overlap, undetermined):
		overlap_feature = re.compile(overlap, re.IGNORECASE)
		original_feature = re.compile(undetermined, re.IGNORECASE)
		# check all positions that the overlap sequence and the original sequence appear. Apparently the positions of original
		# sequence is no less than that of overlap sequence, since the original sequence is contained by the overlap sequence.
		# As long as there's a position that original sequence appears independently, the nonCFds*nonCFds 2cF can be determined
		# as covered by the linked CFds main MEs.
		overlap_positions = set(map(lambda r: r.start() + len(sen_text_2) - 1, overlap_feature.finditer(main_middle_part)))
		original_positions = set(map(lambda r: r.start(), original_feature.finditer(main_middle_part)))
		if len(original_positions - overlap_positions) > 0:
			return NO_OVERLAP
		else:
			return OVERLAP

	for undetermined_2cF in undetermined_faults:
		sen_text_1 = undetermined_2cF.comps['comp1'].Sen
		sen_text_2 = undetermined_2cF.comps['comp2'].Sen
		# check the FP2+FP1 type overlapping
		if sen_text_1[0] == sen_text_2[-1]:
			undetermined_text = sen_text_1 + 'r' + sen_text_1[-1]
			overlap_text = sen_text_2[:-1] + undetermined_text
			if not check_overlapping(overlap_text, undetermined_text):
				undetermined_redundancy_pool.add(undetermined_2cF)
		# check the FP1+FP2 type overlapping
		elif sen_text_2[0] == sen_text_1[-1]:
			undetermined_text = sen_text_2 + 'r' + sen_text_2[-1]
			overlap_text = sen_text_1[:-1] + undetermined_text
			if not check_overlapping(overlap_text, undetermined_text):
				undetermined_redundancy_pool.add(undetermined_2cF)
		# if the initial and final states of FP1 and FP2 are not consistent, the overlapping cannot realize.
		# The 2cF is redundant as long as one of the FPs is included in the main MEs' middle part.
		elif (sen_text_1 in main_middle_part) or (sen_text_2 in main_middle_part):
			undetermined_redundancy_pool.add(undetermined_2cF)

	return {'degenerated_seq': degenerated_seq_pool - degenerated_redundant_pool, 'sf_remainder_seq': sf_seq_pool - sf_redundancy_pool,
			'undetermined_fault': undetermined_faults - undetermined_redundancy_pool}


def nonCFds_constructor(degenerated_pool, sf_pool, undetermined_faults):
	# since the sequences in unlinked_pool are unnecessary to abide any 2cF3 or 2cF2aa conditions, use 00/11 ME structure
	# to make sure no redundant sequence will be contained. Besides, the instant-detect rule still need to be obeyed for
	# nonCFds*nonCFds-type 2cFs
	# vertex_2cF_pool_0 = define_vertices(unlinked_pool['Init_0'])
	# vertex_2cF_pool_1 = define_vertices(unlinked_pool['Init_1'])
	# vertex_scf_pool = define_vertices(unlinked_pool['Init_-1'])

	# build "00" ME
	unlinked_me_00_text = construct_unlinked_element(vertex_2cF_pool_0, vertex_scf_pool, 'Init_0')
	unlinked_me_00 = UnlinkedElementsBuilder.terminal_decorator(unlinked_me_00_text, 'Init_0')
	# build "11" ME
	unlinked_me_11_text = construct_unlinked_element(vertex_2cF_pool_1, vertex_scf_pool, 'Init_1')
	unlinked_me_11 = UnlinkedElementsBuilder.terminal_decorator(unlinked_me_11_text, 'Init_1')

	return {'00_me': unlinked_me_00, '11_me': unlinked_me_11}


def scf_constructor(scf_pool):
	vertex_sf_pool = define_vertices(scf_pool)
	scf_me_candidates = set()

	def get_scf_me_candidates(init, vertex_candidate_pool):
		initial_vertex = UnlinkedElementsBuilder.get_vertex_winner(vertex_candidate_pool, set(),
																   CoverageVertex({'coverage': [], 'diff': -1}), init)
		vertex_candidate_pool -= {initial_vertex}
		coverage_chain = initial_vertex.get_march_sequence()

		while len(vertex_candidate_pool) > 0:
			build_result = build_coverage_chain(coverage_chain, -1, vertex_candidate_pool, UnlinkedElementsBuilder,
												set(), initial_vertex, init)
			vertex_candidate_pool -= build_result[0]
			coverage_chain += build_result[1]
		# make sure that the SF me also starts with read operation, since terminal decorator is not applied
		if coverage_chain[1] != 'r':
			coverage_chain = coverage_chain[0] + 'r' + coverage_chain

		return MarchElement(coverage_chain[1:])

	scf_me_candidates.add(get_scf_me_candidates('Init_0', copy.deepcopy(vertex_sf_pool)))
	scf_me_candidates.add(get_scf_me_candidates('Init_1', copy.deepcopy(vertex_sf_pool)))

	return scf_me_candidates


if __name__ == '__main__':
	os.chdir("../")
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	degenerated_2cFs = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'],
											classified_pool['2cF_CFds']['unlinked'])
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], degenerated_2cFs)
	flat_SF_pool = flatten_sf_pool(filtered_SF_pool)

	seq_pool = create_sequence_pool(flat_SF_pool, degenerated_2cFs, classified_pool['2cF_CFds']['linked'],
									classified_pool['2cF_nonCFds_included']['nonCFds_nonCFds'])

	if len(seq_pool['unlinked']['Init_0']) + len(seq_pool['unlinked']['Init_1']) > 0:
		nonCFds_constructor(seq_pool['unlinked'])
	elif len(seq_pool['unlinked']['Init_-1']) > 0:
		scf_constructor(seq_pool['unlinked']['Init_-1'])
	pass

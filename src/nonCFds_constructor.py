from linked_CFds_ass_constructor import *

OVERLAP = True
NO_OVERLAP = False


class UnlinkedElementsBuilder:
	@staticmethod
	def get_vertex_winner(vertex_candidates: set, aux_pool: set, last_winner: CoverageVertex, init: str):

		def check_states(vertex):
			seq = vertex.get_march_segment()
			if len(seq) > 0:
				return seq[0] + seq[-1]
			else:
				return ''

		aux_pool_dict = {'00': set(), '01': set(), '10': set(), '11': set()}
		for feature in aux_pool_dict.keys():
			aux_feature_pool = set(filter(lambda v: check_states(v) == feature, aux_pool))
			aux_pool_dict[feature].update(aux_feature_pool)

		def calculate_coverage(vertex):
			belongings = []
			for seq_obj in vertex.coverage:
				belongings.append(getattr(seq_obj, 'belongings', []))

			return sum(map(lambda b: len(b), belongings))

		def coverage_priority_filter(upper_pool):
			priority_pool = set()
			if len(upper_pool) > 0:
				priority_pool = {max(upper_pool, key=lambda v: calculate_coverage(v))}
			if len(priority_pool) == 0:
				return upper_pool - priority_pool
			else:
				return next(iter(priority_pool))

		def state_priority_filter(upper_pool: set, priority_feature: str):
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

		# check the coverage of the nonCFds*nonCFds faults first, choose the vertex that covers most faults first
		primary_result = coverage_priority_filter(vertex_candidates)
		if isinstance(primary_result, CoverageVertex):
			return primary_result

		# the initial state-based priority
		feature_dict = {'Init_0': ['01', '11', '00', '10'], 'Init_1': ['10', '00', '11', '01']}
		feature_list = feature_dict[init]

		state_history = check_states(last_winner)
		if (state_history == feature_list[0]) or (state_history == feature_list[1]):
			secondary_result = state_priority_filter(primary_result, feature_list[1])
			if isinstance(secondary_result, CoverageVertex):
				return secondary_result

			tertiary_result = state_priority_filter(secondary_result, feature_list[3])
			if isinstance(tertiary_result, CoverageVertex):
				return tertiary_result

			return next(iter(tertiary_result))
		else:
			secondary_result = state_priority_filter(primary_result, feature_list[2])
			if isinstance(secondary_result, CoverageVertex):
				return secondary_result

			tertiary_result = state_priority_filter(secondary_result, feature_list[0])
			if isinstance(tertiary_result, CoverageVertex):
				return tertiary_result

			return next(iter(tertiary_result), CoverageVertex({'coverage': [], 'diff': -1}))

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

			if main_element.content[0] != 'r':
				main_element.content = 'r' + chain[0] + main_element.content
				main_element.update_states()

		return main_element


# end of UnlinkedElementsBuilder definition


def construct_degenerated_segment(original_vertex_pool, degenerated_vertex_pool, degenerated_aux_vertex_pool, undetermined_fault_pool, determined_vertices, init):
	initial_vertex = UnlinkedElementsBuilder.get_vertex_winner(degenerated_vertex_pool, degenerated_aux_vertex_pool, CoverageVertex({'coverage': [], 'diff': -1}), init)
	target_vertex = initial_vertex

	def get_undetermined_adjacency(vertex_obj):
		adjacency = set()
		for seq_obj in vertex_obj.coverage:
			belongings = set(getattr(seq_obj, 'belongings', []))
			adjacency.update(belongings)

		return adjacency

	# the protection is not long enough at the beginning, so the initial vertex cannot be removed here
	max_check_range = max(map(lambda v: len(v.coverage[0].seq_text), original_vertex_pool))
	coverage_chain = initial_vertex.get_march_segment()
	undetermined_finish_flag = False

	while len(degenerated_vertex_pool) > 0:
		build_result = build_coverage_chain(coverage_chain, -1, degenerated_vertex_pool, UnlinkedElementsBuilder, degenerated_aux_vertex_pool, initial_vertex, init)
		target_vertex = build_result[0][0]

		# skip the vertex only comes from undetermined seq pool after all undetermined faults are covered
		if undetermined_finish_flag and target_vertex.coverage[0].undetermined_tag:
			degenerated_vertex_pool.remove(target_vertex)
			continue

		# the protection of the element. To avoid the nonCFds*nonCFds 2cF being disturbed by the mis-sensitization of the 2cF
		# at the tail of last ME and the beginning of current ME, a protection segment needs to be set before the real build
		# process begins. In the protection segment, chosen vertices will not be removed from the pool after operations are appended.
		chain_check_range = len(coverage_chain) - (len(target_vertex.get_march_segment()) - target_vertex.diff) + 1
		# if the protection segment is longer than the longest sequence in the undetermined seq pool, the subsequent operations
		# can be normally built without the disturbance.
		if chain_check_range >= max_check_range:
			if not undetermined_finish_flag:
				undetermined_adjacency = get_undetermined_adjacency(target_vertex)
				undetermined_fault_pool[init] -= undetermined_adjacency
				# if the adjacency is not empty, it means that the chosen vertex comes from undetermined sequence pool,
				# need to be used in sf in-situ filter
				if len(undetermined_adjacency) > 0:
					determined_vertices.add(target_vertex)

			for vertex in build_result[0]:
				if vertex in degenerated_vertex_pool:
					degenerated_vertex_pool.remove(vertex)
				elif vertex in degenerated_aux_vertex_pool:
					degenerated_aux_vertex_pool.remove(vertex)

			if len(undetermined_fault_pool[init]) == 0:
				undetermined_finish_flag = True

		coverage_chain += build_result[1]
	return coverage_chain, target_vertex


def sf_in_situ_filter(sf_init_pool, sf_aux_pool, determined_vertex_pool, init):
	determined_seq_pool = {'Init_0': set(), 'Init_1': set(), 'Init_-1': set()}
	redundant_init_pool = set()
	redundant_aux_pool = set()

	for v_obj in determined_vertex_pool:
		for seq_obj in v_obj.coverage:
			init_key = 'Init_' + seq_obj.ass_init
			determined_seq_pool[init_key].add(seq_obj)

	for v_obj in sf_init_pool:
		for seq_obj in v_obj.coverage:
			find_result = find_inclusive_seq(seq_obj, determined_seq_pool[init], {'ass_init', 'belongings', 'undetermined_tag', 'detect_tag', 'dr_tag', 'nest_tag'})
			if not isinstance(find_result, type(DIFFERENT)):
				if (find_result.seq_text == seq_obj.seq_text) and seq_obj.detect_tag:
					find_result.detect_tag |= seq_obj.detect_tag
					find_result.dr_tag |= seq_obj.dr_tag
					find_result.nest_tag = seq_obj.nest_tag

				redundant_init_pool.add(seq_obj)

	for v_obj in sf_aux_pool:
		for seq_obj in v_obj.coverage:
			find_result = DIFFERENT
			for determined_key in determined_seq_pool.keys():
				find_result = find_inclusive_seq(seq_obj, determined_seq_pool[determined_key], {'ass_init', 'belongings', 'undetermined_tag', 'detect_tag', 'dr_tag', 'nest_tag'})
				if not isinstance(find_result, type(DIFFERENT)):
					break
			if not isinstance(find_result, type(DIFFERENT)):
				if (find_result.seq_text == seq_obj.seq_text) and seq_obj.detect_tag:
					find_result.detect_tag |= seq_obj.detect_tag
					find_result.dr_tag |= seq_obj.dr_tag
					find_result.nest_tag = seq_obj.nest_tag

				redundant_aux_pool.add(seq_obj)

	sf_init_pool -= redundant_init_pool
	sf_aux_pool -= redundant_aux_pool
	return


def construct_nonCFds_element(vertex_pool, vertex_aux_pool, sf_init_vertex_pool, sf_aux_vertex_pool, undetermined_fault_pool, init):
	vertex_candidate_pool = vertex_pool
	# check if the degenerated and undetermined pools are empty. If so, the initial vertex starts from sf pool.
	if len(vertex_pool) > 0:
		# build from degenerated & undetermined pool
		determined_vertices = set()
		construct_result = construct_degenerated_segment(vertex_pool, vertex_candidate_pool, vertex_aux_pool, undetermined_fault_pool, determined_vertices, init)
		coverage_chain = construct_result[0]
		target_vertex = construct_result[1]
		# SFs still need to be filtered by the used nonCFds*nonCFds faults, before the SF sequences are added into coverage chain.
		sf_in_situ_filter(sf_init_vertex_pool, sf_aux_vertex_pool, determined_vertices, init)
	else:
		initial_vertex = UnlinkedElementsBuilder.get_vertex_winner(sf_init_vertex_pool, sf_aux_vertex_pool, CoverageVertex({'coverage': [], 'diff': -1}), init)
		target_vertex = initial_vertex
		coverage_chain = initial_vertex.get_march_segment()
		if initial_vertex in sf_init_vertex_pool:
			sf_init_vertex_pool.remove(initial_vertex)
		elif initial_vertex in sf_aux_vertex_pool:
			sf_aux_vertex_pool.remove(initial_vertex)

	# build from sf pool
	while len(sf_init_vertex_pool) > 0:
		build_result = build_coverage_chain(coverage_chain, -1, sf_init_vertex_pool, UnlinkedElementsBuilder, sf_aux_vertex_pool, target_vertex, init)
		for vertex in build_result[0]:
			if vertex in sf_init_vertex_pool:
				sf_init_vertex_pool.remove(vertex)
			elif vertex in sf_aux_vertex_pool:
				sf_aux_vertex_pool.remove(vertex)

		coverage_chain += build_result[1]

	return coverage_chain


# filter the degenerated 2cFs, SFs and the nonCFds*nonCFds faults that are possibly covered by the linked CFds main MEs.
# Besides the segments appended in the linked main MEs, other faults may be covered because of the nest between two segments
def inter_ME_filter(main_elements, main_middle_part, degenerated_seq_pool, sf_seq_pool, undetermined_faults):
	main_element_contents = {'01_me': '', '10_me': ''}
	for element_key in main_elements.keys():
		main_element_contents[element_key] = main_elements[element_key].content[1:]

	def search_seq_redundancy(candidate_seq, init):
		# if the sequence comes from a nonCFds FP, the match candidate text should be seq_text and a read-verify operation,
		# and the candidate text should be included in both ME, since the main MEs are 01/10 type.
		if candidate_seq.detect_tag:
			candidate_text = candidate_seq.seq_text + 'r' + candidate_seq.seq_text[-1]
		else:
			candidate_text = candidate_seq.seq_text

		search_result = list(filter(lambda m: candidate_text in m, main_element_contents))
		if (init == 'Init_-1') and (len(search_result) > 0):
			return REDUNDANT
		elif len(search_result) > 1:
			return REDUNDANT

		return NOT_REDUNDANT

	# check the redundancy according to their corresponding test conditions
	# consider the 2cFs in degenerated pool
	degenerated_redundant_pool = {'Init_0': set(), 'Init_1': set(), 'Init_-1': set()}
	for degenerated_key in degenerated_seq_pool.keys():
		for degenerated_seq in degenerated_seq_pool[degenerated_key]:
			redundancy = search_seq_redundancy(degenerated_seq, degenerated_key)
			if redundancy:
				degenerated_redundant_pool[degenerated_key].add(degenerated_seq)

	# consider SFs
	sf_redundancy_pool = {'Init_0': set(), 'Init_1': set(), 'Init_-1': set()}
	for sf_key in sf_seq_pool.keys():
		for sf_seq in sf_seq_pool[sf_key]:
			redundancy = search_seq_redundancy(sf_seq, sf_key)
			if redundancy:
				sf_redundancy_pool[sf_key].add(sf_seq)

	# consider undetermined faults, check whether one of the FP is included in the main ME, and no overlapping
	undetermined_redundancy_pool = set()

	def check_overlapping(overlap, undetermined):
		overlap_feature = re.compile(overlap, re.IGNORECASE)
		original_feature = re.compile(undetermined, re.IGNORECASE)
		# check all positions that the overlap sequence and the original sequence appear. Apparently the positions of original
		# sequence is no less than that of overlap sequence, since the original sequence is contained by the overlap sequence.
		# As long as there's a position that original sequence appears independently, the nonCFds*nonCFds 2cF can be determined
		# as covered by the linked CFds main MEs.
		overlap_positions = set(
			map(lambda r: r.start() + len(sen_text_FP2) - 1, overlap_feature.finditer(main_middle_part)))
		original_positions = set(map(lambda r: r.start(), original_feature.finditer(main_middle_part)))
		# if original_positions is 0, it means that the sequence is not included by the main middle part, seen as OVERLAP (not redundant)
		if 0 < len(original_positions) > len(overlap_positions):
			return NO_OVERLAP
		else:
			return OVERLAP

	for undetermined_2cF in undetermined_faults:
		fp1 = undetermined_2cF.comps['comp1']
		fp2 = undetermined_2cF.comps['comp2']
		sen_text_FP1 = fp1.aInit if fp1.CFdsFlag == 1 else fp1.vInit + undetermined_2cF.comps['comp1'].Sen
		sen_text_FP2 = fp2.aInit if fp2.CFdsFlag == 1 else fp2.vInit + undetermined_2cF.comps['comp2'].Sen
		check_result = {'FP1+FP2': sen_text_FP2 not in main_middle_part, 'FP2+FP1': sen_text_FP1 not in main_middle_part}
		if (sen_text_FP1 in main_middle_part) or (sen_text_FP2 in main_middle_part):
			# check the FP1+FP2 type overlapping
			if sen_text_FP2[0] == sen_text_FP1[-1]:
				undetermined_texts = [sen_text_FP2 + 'r' + sen_text_FP2[-1]]
				# consider the nest sensitization cases. If the FP is a donor, the sensitization sequence in main MEs can be
				# the nest style. Check them both.
				if fp2.nestSenFlag == 'donor':
					undetermined_texts.append(sen_text_FP2 + sen_text_FP2[1:] + 'r' + sen_text_FP2[-1])
				for undetermined_text in undetermined_texts:
					overlap_text = sen_text_FP1[:-1] + undetermined_text
					temp_result = check_overlapping(overlap_text, undetermined_text)
					check_result['FP1+FP2'] = temp_result
					# the check finishes as long as any of the check result is NO_OVERLAP
					if not temp_result:
						break
			# check the FP2+FP1 type overlapping
			if sen_text_FP1[0] == sen_text_FP2[-1]:
				undetermined_texts = [sen_text_FP1 + 'r' + sen_text_FP1[-1]]
				if fp1.nestSenFlag == 'donor':
					undetermined_texts.append(sen_text_FP1 + sen_text_FP1[1:] + 'r' + sen_text_FP1[-1])
				for undetermined_text in undetermined_texts:
					overlap_text = sen_text_FP2[:-1] + undetermined_text
					temp_result = check_overlapping(overlap_text, undetermined_text)
					check_result['FP2+FP1'] = temp_result
					if not temp_result:
						break
			# if the initial and final states of FP1 and FP2 are not consistent, the overlapping cannot realize.
			# The 2cF is redundant as long as one of the FPs is included in the main MEs' middle part.
			if NO_OVERLAP in check_result.values():
				undetermined_redundancy_pool.add(undetermined_2cF)

	filter_result = {'degenerated_seq': {}, 'sf_seq': {}, 'undetermined_faults': {}}
	for init_key in {'Init_0', 'Init_1', 'Init_-1'}:
		filter_result['degenerated_seq'][init_key] = degenerated_seq_pool[init_key] - degenerated_redundant_pool[
			init_key]
		filter_result['sf_seq'][init_key] = sf_seq_pool[init_key] - sf_redundancy_pool[init_key]
	filter_result['undetermined_faults'] = undetermined_faults - undetermined_redundancy_pool

	return filter_result


def nonCFds_constructor(degenerated_seq_pool, undetermined_fault_pool, sf_seq_pool):
	# since the sequences in unlinked_pool are unnecessary to abide any 2cF3 or 2cF2aa conditions, use 00/11 ME structure
	# to make sure no redundant sequence will be contained. Besides, the instant-detect rule still need to be obeyed for
	# nonCFds*nonCFds-type 2cFs

	# At the beginning of the construction, the sf pool should not be included, since it still need to be filtered by the
	# undetermined faults.

	# Add both FPs of each nonCFds*nonCFds 2cF to the candidate vertex pool, merge with the degenerated sequences.
	# Check the adjacent matrix after every time after sequence comes from nonCFds*nonCFds is chosen to build
	candidate_seq_pool = {'Init_0': copy.deepcopy(degenerated_seq_pool['Init_0']), 'Init_1': copy.deepcopy(degenerated_seq_pool['Init_1']),
						  'Init_-1': copy.deepcopy(degenerated_seq_pool['Init_-1'])}
	undetermined_seq_pool = {'Init_0': set(), 'Init_1': set()}
	create_undetermined_sequences(undetermined_fault_pool, undetermined_seq_pool)
	undetermined_fault_temp = {'Init_0': set(), 'Init_1': set()}

	for undetermined_key in undetermined_seq_pool.keys():
		for seq_obj in undetermined_seq_pool[undetermined_key]:
			undetermined_fault_temp[undetermined_key].update(seq_obj.belongings)
			merge_result = find_identical_objs(seq_obj, candidate_seq_pool[undetermined_key], {'ass_init', 'belongings', 'undetermined_tag', 'detect_tag', 'dr_tag', 'nest_tag'})
			if isinstance(merge_result, type(DIFFERENT)):
				# set a flag to clarify it only comes from undetermined fault pool
				setattr(seq_obj, 'undetermined_tag', True)
				candidate_seq_pool[undetermined_key].add(seq_obj)
			else:
				setattr(merge_result, 'belongings', seq_obj.belongings)
				setattr(seq_obj, 'undetermined_tag', False)
				merge_result.detect_tag |= seq_obj.detect_tag
				merge_result.dr_tag |= seq_obj.dr_tag
				merge_result.nest_tag = seq_obj.nest_tag

	vertex_pool_0 = define_vertices(candidate_seq_pool['Init_0'])
	vertex_pool_1 = define_vertices(candidate_seq_pool['Init_1'])
	vertex_scf_pool = define_vertices(candidate_seq_pool['Init_-1'])

	sf_vertex_pool_0 = define_vertices(sf_seq_pool['Init_0'])
	sf_vertex_pool_1 = define_vertices(sf_seq_pool['Init_1'])
	sf_aux_vertex_pool = define_vertices(sf_seq_pool['Init_-1'])

	# build "00" ME
	unlinked_me_00_text = construct_nonCFds_element(vertex_pool_0, vertex_scf_pool, sf_vertex_pool_0, sf_aux_vertex_pool, undetermined_fault_temp, 'Init_0')
	unlinked_me_00 = UnlinkedElementsBuilder.terminal_decorator(unlinked_me_00_text, 'Init_0')
	# build "11" ME
	unlinked_me_11_text = construct_nonCFds_element(vertex_pool_1, vertex_scf_pool, sf_vertex_pool_1, sf_aux_vertex_pool, undetermined_fault_temp, 'Init_1')
	unlinked_me_11 = UnlinkedElementsBuilder.terminal_decorator(unlinked_me_11_text, 'Init_1')

	# all 2cFs and SFs with certain initial states are covered in the nonCFds ME, only the remainders of sf_aux_pool
	# and vertex_scf_pool need to be covered in SCF ME. Additionally, since the sequences have been checked without any inclusion,
	# the two aux pools can be union directly
	return {'00_me': unlinked_me_00, '11_me': unlinked_me_11}, sf_aux_vertex_pool.union(vertex_scf_pool)


def scf_constructor(vertex_scf_pool):
	scf_me_candidates = set()

	def get_scf_me_candidates(init, vertex_candidate_pool):
		initial_vertex = UnlinkedElementsBuilder.get_vertex_winner(vertex_candidate_pool, set(), CoverageVertex({'coverage': [], 'diff': -1}), init)
		vertex_candidate_pool -= {initial_vertex}
		coverage_chain = initial_vertex.get_march_segment()

		while len(vertex_candidate_pool) > 0:
			build_result = build_coverage_chain(coverage_chain, -1, vertex_candidate_pool, UnlinkedElementsBuilder, set(), initial_vertex, init)
			vertex_candidate_pool -= set(build_result[0])
			coverage_chain += build_result[1]
		# make sure that the SF me also starts with read operation, since terminal decorator is not applied
		if coverage_chain[1] != 'r':
			coverage_chain = coverage_chain[0] + 'r' + coverage_chain

		return MarchElement(coverage_chain[1:])

	scf_me_candidates.add(get_scf_me_candidates('Init_0', copy.deepcopy(vertex_scf_pool)))
	scf_me_candidates.add(get_scf_me_candidates('Init_1', copy.deepcopy(vertex_scf_pool)))

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

	if len(seq_pool['undetermined_faults']) + sum(map(lambda p: len(seq_pool['sf_seq'][p]), seq_pool['sf_seq'].keys())) + sum(map(lambda p: len(seq_pool['degenerated_seq'][p]), seq_pool['degenerated_seq'].keys())) > 0:
		nonCFds_result = nonCFds_constructor(seq_pool['degenerated_seq'], seq_pool['undetermined_faults'], seq_pool['sf_seq'])
	if len(nonCFds_result[1]) > 0:
		scf_constructor(seq_pool['unlinked']['Init_-1'])
	pass

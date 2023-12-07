from linked_main_constructor import *


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
		initial_vertex = LinkedMainElementsBuilder.get_vertex_winner(vertex_candidate_pool, set(),
																	 CoverageVertex({'coverage': [], 'diff': -1}),
																	 'Init_-1')
		vertex_candidate_pool -= {initial_vertex}
		coverage_chain = initial_vertex.get_march_sequence()

		while len(vertex_candidate_pool) > 0:
			build_result = build_coverage_chain(coverage_chain, set(), -1, vertex_candidate_pool, set(), initial_vertex, 'Init_-1',
												LinkedMainElementsBuilder)
			vertex_candidate_pool -= build_result[0]
			coverage_chain += build_result[1]

		# the terminal of tail-cover ME need to be decorated to match the original terminal states of the ME where
		# tail-cover seq in. The initial state of the v-cell needs to be consistent with the original ME, since the AO is reversed.
		# for instance, assume a main ME has a tail-cover "1w0r0w1", a tail-cover ME "down,r1,w0,r0,w1" cannot cover
		# the fault <1w0r0w1;0/1/->*<1r1w1w0;1/0/->, since the tail-cover ME is "11" type, the initial state of v-cell is
		# different from the original "01" type ME. The consistency is natural in 2-operation fault list, while not in 3 and more
		# operation fault list.
		match coverage_chain[0] + coverage_chain[-1]:
			case '00':
				tail_cover_me = 'r1w' + coverage_chain
			case '11':
				tail_cover_me = 'r0w' + coverage_chain
			case _:
				if not coverage_chain[1:].startswith('r'):
					tail_cover_me = 'r' + coverage_chain
				else:
					tail_cover_me = coverage_chain[1:]

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
	violation_pool.append([get_violated_seq_texts(element_text_list), element_text_list[1], element_text_list[0]])
	element_text_list.reverse()
	violation_pool.append([get_violated_seq_texts(element_text_list), element_text_list[1], element_text_list[0]])

	if max(map(lambda t: len(t[0]), violation_pool)) > 0:
		return violation_pool
	else:
		return NOT_FOUND


def build_single_sensitization_chain(chain: str, vertex_pool: set, covered_vertex_pool: set):
	# an individual build function to make sure the sequences are sensitized only once in odd-sensitization ME
	# a set records the covered vertices by this build process
	for v_obj in vertex_pool:
		v_obj.diff = calculate_diff_value(chain, v_obj)

	vertex_candidates = sorted(vertex_pool, key=lambda v: v.diff)

	for vertex_candidate in vertex_candidates:
		candidate_segment = vertex_candidate.get_march_sequence()
		if vertex_candidate.diff < len(candidate_segment):
			candidate_appendix = candidate_segment[-vertex_candidate.diff:]
		else:
			candidate_appendix = 'w' + candidate_segment

		for covered_vertex in covered_vertex_pool:
			segment_under_check = covered_vertex.get_march_sequence()
			search_range = min(len(segment_under_check) - 2, len(chain))
			chain_under_check = chain[-search_range:] + candidate_appendix
			if segment_under_check in chain_under_check:
				break
		else:
			return {vertex_candidate}, candidate_appendix
	else:
		return {}, ''


def construct_odd_sensitization_elements(odd_violation, tail_cover):
	# for each sequence object in odd_violation, it should be regarded as a CFds, since it just violates the
	# odd-sensitization condition for CFds
	violation_pool = set()
	violation_me_candidates = []

	if isinstance(odd_violation, list):
		for odd_case in odd_violation:
			# list for multiple odd-sensitization MEs
			odd_mes = []
			violation_pool.clear()
			# if an odd violation has been in tail-cover ME, it will be sensitized individually,
			# so no need to add it in this ME
			odd_case[0] -= tail_cover

			for seq_obj in odd_case[0]:
				violation_seq = Sequence(seq_obj.__dict__.copy())
				violation_seq.detect_tag = False
				violation_vertex = CoverageVertex({'diff': -1, 'coverage': [violation_seq]})
				violation_pool.add(violation_vertex)

			vertex_candidate_pool = copy.deepcopy(violation_pool)
			covered_vertex_pool = set()

			max_range = max(set(map(lambda s: len(s.seq_text), odd_case[0]))) - 2
			search_range = min(max_range, len(odd_case[1]))
			coverage_chain = odd_case[1][-search_range:] + 'r' + odd_case[1][-1]

			# the initial coverage chain may also cover a sequence
			initial_vertex = set(filter(lambda v: v.coverage[0].seq_text in coverage_chain, vertex_candidate_pool))
			vertex_candidate_pool -= initial_vertex
			covered_vertex_pool.update(initial_vertex)

			while len(vertex_candidate_pool) > 0:

				build_result = build_single_sensitization_chain(coverage_chain, vertex_candidate_pool, covered_vertex_pool)
				if len(build_result[0]) > 0:
					# if there is sequence that satisfy the single-sensitization, append it normally
					vertex_candidate_pool -= build_result[0]
					covered_vertex_pool.update(build_result[0])
					coverage_chain += build_result[1]
				else:
					# if not, save and end this odd-sensitization ME
					odd_mes.append(coverage_chain[search_range:])
					covered_vertex_pool.clear()
					# initialize the coverage chain for next ME
					max_range = max(set(map(lambda v: len(v.get_march_sequence()), vertex_candidate_pool))) - 2
					search_range = min(max_range, len(odd_mes[-1]))
					coverage_chain = odd_mes[-1][-search_range:] + 'r' + odd_mes[-1][-1]

					initial_vertex = set(filter(lambda v: v.coverage[0].seq_text in coverage_chain, vertex_candidate_pool))
					vertex_candidate_pool -= initial_vertex
					covered_vertex_pool.update(initial_vertex)

			odd_mes.append(coverage_chain[search_range:])

			# each odd_sensitization candidate has to follow the tied precedent element, or the odd sensitization may be destroyed
			violation_me_candidates.append((odd_mes, odd_case[1], odd_case[2]))

		odd_sensitization = min(violation_me_candidates, key=lambda c: sum(map(lambda m: len(m), c[0])))
		return odd_sensitization
	else:
		return (NO_ELEMENT,)


def construct_ass_elements(main_elements, filtered_sequence_pool, original_sequence_pool):
	ass_elements = {'tail_cover_me': MarchElement(''), 'odd_sensitization_me': MarchElement('')}
	# check and build the ME for tail-cover first
	tail_cover = check_tail_cover(list(filter(lambda m: m.tail_tag, main_elements.values())), filtered_sequence_pool)
	tail_cover_me_text = construct_tail_cover_elements(tail_cover)

	if isinstance(tail_cover_me_text, str):
		tail_cover_me = MarchElement(tail_cover_me_text)
		tail_cover_me.ass_tag = True
		ass_elements['tail_cover_me'] = tail_cover_me

	# check and build the ME for odd sensitization violation
	odd_violation = check_odd_sensitization(main_elements, original_sequence_pool)
	construct_result = construct_odd_sensitization_elements(odd_violation, tail_cover)
	odd_sensitization_me_texts = construct_result[0]
	odd_sensitization_mes = []
	if isinstance(odd_sensitization_me_texts, list):
		for odd_sensitization_me_text in odd_sensitization_me_texts:
			odd_sensitization_me = MarchElement(odd_sensitization_me_text)
			odd_sensitization_me.ass_tag = True
			odd_sensitization_me.tied_element = []
			odd_sensitization_mes.append(odd_sensitization_me)

		odd_sensitization_mes.reverse()
		odd_sensitization_mes[0].tied_element.extend(odd_sensitization_mes[1:])

		for tied_element_text in construct_result[1:]:
			odd_sensitization_mes[0].tied_element.extend(list(filter(lambda m: m.content == tied_element_text, main_elements.values())))

		ass_elements['odd_sensitization_me'] = odd_sensitization_mes[0]

	return ass_elements


def linked_CFds_constructor(filtered_linked_seq_pool, original_linked_fault_pool):
	# For 01/10 type main ME, all linked seq need to be added for constructing; for 00/11 type, init_0 and init_1 pool
	# should be added for constructing separately. however, the different seq set in 00 and 11 ME may introduce the
	# violation to sensitization order. As a result, we only use the 01/10 type ME, which means that all sequences in
	# linked pool should be contained in both 2 ME, and the ass_init of them can be ignored
	filtered_union_pool = get_linked_CFds_union(filtered_linked_seq_pool['Init_0'], filtered_linked_seq_pool['Init_1'])
	vertex_union = define_vertices(filtered_union_pool)
	main_mes = construct_main_elements(vertex_union)

	linked_CFds_pool = copy.deepcopy(original_linked_fault_pool)
	original_linked_seq_pool = {'Init_0': set(), 'Init_1': set()}

	for linked_CFds in linked_CFds_pool:
		for comp_obj in linked_CFds.comps.values():
			fault_sequence = Sequence(get_sequence_properties(comp_obj))
			init_key = 'Init_' + fault_sequence.ass_init
			if isinstance(find_identical_objs(fault_sequence, original_linked_seq_pool[init_key], {}), type(DIFFERENT)):
				original_linked_seq_pool[init_key].add(copy.deepcopy(fault_sequence))

	original_union_pool = get_linked_CFds_union(original_linked_seq_pool['Init_0'], original_linked_seq_pool['Init_1'])
	ass_mes = construct_ass_elements(main_mes, filtered_union_pool, original_union_pool)

	return {'main_me': main_mes, 'ass_me': ass_mes}


if __name__ == '__main__':
	os.chdir("../")
	sys.path.append("src/")
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_2cF_pool = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'],
											 classified_pool['2cF_CFds']['unlinked'])
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], filtered_2cF_pool)
	flat_SF_pool = flatten_sf_pool(filtered_SF_pool)

	seq_pool = create_sequence_pool(flat_SF_pool, filtered_2cF_pool, classified_pool['2cF_CFds']['linked'])

	if len(seq_pool['linked']['Init_0']) + len(seq_pool['linked']['Init_1']) > 0:
		linked_CFds_constructor(seq_pool['linked'], classified_pool['2cF_CFds']['linked'])

	pass

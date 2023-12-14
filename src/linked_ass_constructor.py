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
	if isinstance(tail_cover, set):
		for seq_obj in tail_cover:
			tail_cover_seq = Sequence(seq_obj.__dict__.copy())
			# tail cover is for CFds
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
			build_result = build_coverage_chain(coverage_chain, -1, vertex_candidate_pool, set(), initial_vertex, 'Init_-1',
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


def check_odd_sensitization(main_elements, sequence_pool: set):
	# find all sequences that are sensitized even-number times in ME, which violate the odd-sensitization condition
	element_text_list = list(map(lambda e: e.content, main_elements.values()))
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


def check_vertices_covered_by_appendix(chain: str, appendix: str, vertex_pool: set):
	covered_vertices = set()
	for vertex in vertex_pool:
		segment_under_check = vertex.get_march_sequence()
		search_range = min(len(segment_under_check) - 2, len(chain))
		chain_under_check = chain[-search_range:] + appendix
		if segment_under_check in chain_under_check:
			covered_vertices.add(vertex)

	return covered_vertices


def construct_odd_sensitization_elements(odd_violation, tail_cover):
	# for each sequence object in odd_violation, it should be regarded as a CFds, since it just violates the
	# odd-sensitization condition for CFds
	violation_pool = set()
	violation_me_candidates = []

	if isinstance(odd_violation, list):
		for odd_case in odd_violation:
			# list for multiple odd-sensitization MEs
			odd_mes = []
			# in a specific order of main MEs, the odd-sensitization ME is empty
			if len(odd_case[0]) < 1:
				odd_mes.append('')
				violation_me_candidates.append((odd_mes, odd_case[1], odd_case[2]))
				continue

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
					covered_vertices = build_result[0].union(check_vertices_covered_by_appendix(coverage_chain, build_result[1], vertex_candidate_pool))
					covered_vertex_pool.update(covered_vertices)
					vertex_candidate_pool -= covered_vertices
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


def check_head_cover(main_elements: dict, chain: str, sequence_pool: set):
	# check whether the head of the ME destroys the sensitization of the first nonCF of middle part. The head operation
	# "r1,w0" or "r0,w1" may introduce new mis-sensitization sequences (sequences that are just sensitized but not detected),
	# which can result in the damage to the initial state of the following normal sequence in MP.
	head_cover = []
	head_cover_candidates = []
	victims = []
	expected_locations = []
	head_decorated_me = next(iter(filter(lambda m: m.head_tag is True, main_elements.values())))
	head_decorated_content = head_decorated_me.content[1] + head_decorated_me.content

	search_range = sorted(set(map(lambda s: len(s.seq_text), sequence_pool)), reverse=True)
	# only two possible head cover sequences, since the head consists of 2 operations
	for head_operation in range(2):
		for location in search_range:
			start_location = 2 * head_operation
			end_location = 2 * head_operation + location
			head_cover_candidate = next(iter(filter(lambda s: s.seq_text == head_decorated_content[start_location:end_location], sequence_pool)), '')
			if isinstance(head_cover_candidate, Sequence):
				head_cover_candidates.append((head_cover_candidate, head_operation))
				break
	if len(head_cover_candidates) < 1:
		return NOT_FOUND

	for candidate, index in head_cover_candidates:
		head_cover_text = candidate.seq_text
		# if there's read operation that follows the head_cover candidate, the cover is harmless, since the head cover
		# will be detected instantly
		if head_decorated_content.find(head_cover_text) == head_decorated_content.find(head_cover_text + 'r' + head_cover_text[-1]):
			continue
		expected_locations.append(chain.find(head_cover_text + 'r' + head_cover_text[-1]))
		for location in search_range:
			start_location = 2 * index + location - 1
			end_location = 2 * (index + location) - 1
			victim = next(iter(filter(lambda s: s.seq_text == head_decorated_content[start_location:end_location], sequence_pool)), '')
			if isinstance(victim, Sequence):
				victims.append(victim)
				break
	if len(victims) < 1:
		return NOT_FOUND

	head_information = list(zip(head_cover_candidates, victims, expected_locations))

	for information in head_information:
		victim_text = information[1].seq_text
		if chain[information[2] + 1 - len(victim_text):information[2] + 1] == victim_text:
			head_cover.append((information[0][0], information[1]))

	if len(head_cover) > 0:
		return head_cover
	else:
		return NOT_FOUND


def construct_head_cover_element(odd_element, main_elements, head_cover):
	if not isinstance(head_cover, list):
		return NO_ELEMENT
	# the head cover element should follow the odd_sensitization elements. If there's no odd elements, follow main elements
	vertex_candidate_pool = set()
	for seq_obj in map(lambda h: h[0], head_cover):
		head_cover_seq = Sequence(seq_obj.__dict__.copy())
		# head cover is for nonCFds
		head_cover_seq.detect_tag = True
		head_cover_vertex = CoverageVertex({'diff': -1, 'coverage': [head_cover_seq]})
		vertex_candidate_pool.add(head_cover_vertex)

	# determine the precedent element
	head_decorated_me = next(iter(filter(lambda m: m.head_tag is True, main_elements.values())))
	if isinstance(odd_element.tied_element, list):
		if len(odd_element.content) > 0:
			precedent = odd_element
		else:
			precedent = odd_element.tied_element[0]

		if precedent.content[-1] != head_decorated_me.content[1]:
			transition_me_dict = {'1': MarchElement('r0,w1'), '0': MarchElement('r1,w0')}
			transition_me = transition_me_dict[head_decorated_me.content[1]]
			transition_me.transition_flag = True
			transition_me.tied_element = [precedent]
			precedent = transition_me
	else:
		# if there's no precedent restriction, choose the main ME that final state equals to the initial state of head decorated ME
		precedent = next(iter(filter(lambda m: m.head_tag is False, main_elements.values())))

	max_range = max(set(map(lambda s: len(s[0].seq_text), head_cover))) - 2
	precedent_text = precedent.content
	search_range = min(max_range, len(precedent_text))

	chain_candidate = []
	# The initial state of the head_cover ME should be consistent with the initial state of the head-decorated main ME.
	# see the explanation in tail_cover part
	coverage_chain = precedent_text[-search_range:] + 'r' + precedent_text[-1]
	vertex_candidate_list = list(vertex_candidate_pool)

	def build_head_cover_chain(candidate_list, chain):
		head_cover_chain = chain
		max_check_range = max(set(map(lambda s: len(s[0].seq_text), head_cover)))
		check_range = min(max_check_range, len(chain))
		for vertex in candidate_list:
			diff = calculate_diff_value(head_cover_chain, vertex)
			segment = vertex.get_march_sequence()
			if diff == len(segment):
				head_cover_chain += 'w' + segment[0]

			victim = next(iter(set(filter(lambda i: i[0].seq_text == vertex.coverage[0].seq_text, head_cover))))[1]
			if victim.seq_text == head_cover_chain[-check_range:]:
				head_cover_chain += 'r' + head_cover_chain[-1]
				continue
			else:
				head_cover_chain += segment[1:]

		return head_cover_chain

	chain_candidate.append(build_head_cover_chain(vertex_candidate_list, coverage_chain)[search_range:])
	vertex_candidate_list.reverse()
	chain_candidate.append(build_head_cover_chain(vertex_candidate_list, coverage_chain)[search_range:])

	chain_winner = min(chain_candidate, key=lambda c: len(c))
	if chain_winner[-1] != head_decorated_me.content[-1]:
		chain_winner += 'w' + head_decorated_me.content[-1]

	head_cover_me = MarchElement(chain_winner)
	head_cover_me.ass_tag = True
	head_cover_me.tied_element = [precedent]
	return head_cover_me


def construct_ass_elements(main_elements, main_middle_part, filtered_sequence_pool, original_sequence_pool):
	ass_elements = {'head_cover_me': MarchElement(''), 'tail_cover_me': MarchElement(''), 'odd_sensitization_me': MarchElement('')}

	# check and build the ME for tail-cover first, it may be used in constructing odd-sensitization MEs
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
	# if the odd_sensitization_me_texts is a list, it means that the odd-sensitization ME exists, or not exist only under
	# one specific order of main MEs (the tied elements need to be arranged)
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

	# check and build the ME for head-cover
	head_cover = check_head_cover(main_elements, main_middle_part, filtered_sequence_pool)
	head_cover_me = construct_head_cover_element(odd_sensitization_mes[0], main_elements, head_cover)
	if isinstance(head_cover_me, MarchElement):
		ass_elements['head_cover_me'] = head_cover_me

	return ass_elements


def linked_CFds_constructor(filtered_linked_seq_pool, original_linked_fault_pool):
	# For 01/10 type main ME, all linked seq need to be added for constructing; for 00/11 type, init_0 and init_1 pool
	# should be added for constructing separately. however, the different seq set in 00 and 11 ME may introduce the
	# violation to sensitization order. As a result, we only use the 01/10 type ME, which means that all sequences in
	# linked pool should be contained in both 2 ME, and the ass_init of them can be ignored
	filtered_union_pool = get_linked_CFds_union(filtered_linked_seq_pool['Init_0'], filtered_linked_seq_pool['Init_1'])
	vertex_union = define_vertices(filtered_union_pool)
	main_result = construct_main_elements(vertex_union)
	main_mes = main_result[0]
	main_coverage_chain = main_result[1]

	linked_CFds_pool = copy.deepcopy(original_linked_fault_pool)
	original_linked_seq_pool = {'Init_0': set(), 'Init_1': set()}

	# for odd sensitization ME, the sequence of all number of operations need to be considered, since the inclusive rule
	# cannot make sure that a shorter sequence is only included once in longer sequence, so even the longer sequence is
	# included odd times, it does not mean that included shorter sequences are also included odd times.
	for linked_CFds in linked_CFds_pool:
		for comp_obj in linked_CFds.comps.values():
			fault_sequence = Sequence(get_sequence_properties(comp_obj))
			init_key = 'Init_' + fault_sequence.ass_init
			if isinstance(find_identical_objs(fault_sequence, original_linked_seq_pool[init_key], {}), type(DIFFERENT)):
				original_linked_seq_pool[init_key].add(copy.deepcopy(fault_sequence))

	original_union_pool = get_linked_CFds_union(original_linked_seq_pool['Init_0'], original_linked_seq_pool['Init_1'])
	ass_mes = construct_ass_elements(main_mes, main_coverage_chain, filtered_union_pool, original_union_pool)

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

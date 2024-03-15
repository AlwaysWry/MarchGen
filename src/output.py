# assign the obtained MEs and output them as standard format
from nonCFds_builder import *


def initial_based_assign(linked_me, scf_me, me_dict, precedent_list, element):
	# if current element is odd-sensitization element, it will be assigned with its two tied main elements in fixed order
	if isinstance(element.tied_element, list):
		precedent = element.tied_element[-1]
		precedent_list.extend(element.tied_element)
		return initial_based_assign(linked_me, scf_me, me_dict, precedent_list, precedent)

	match element.initial_state:
		case '0':
			feature_keys = ['00', '10', '11', '01']
			for key in feature_keys:
				iter_key = iter(me_dict[key])
				precedent = next(iter_key, '')
				if not isinstance(precedent, MarchElement):
					continue
				else:
					if (len(me_dict[key]) > 1) and (precedent in linked_me['main_me'].values()):
						precedent = next(iter_key)
					# if an ME for SFs are visited, the other scf_me candidate can be discarded
					if precedent in scf_me:
						useless_me = next(iter(scf_me - {precedent}))
						me_dict[useless_me.initial_state + useless_me.final_state].remove(useless_me)

					me_dict[key].remove(precedent)

					if key == '11' or key == '01':
						transition_text = 'r1w0'
						transition = MarchElement(transition_text)
						transition.address_order = precedent.address_order
						transition.transition_tag = True
						precedent_list.extend([transition])

					precedent_list.extend([precedent])

					return initial_based_assign(linked_me, scf_me, me_dict, precedent_list, precedent)
			return

		case '1':
			feature_keys = ['11', '01', '00', '10']
			for key in feature_keys:
				iter_key = iter(me_dict[key])
				precedent = next(iter_key, '')
				if not isinstance(precedent, MarchElement):
					continue
				else:
					# if an ME for SFs are visited, the other scf_me candidate can be discarded
					if precedent in scf_me:
						scf_me -= {precedent}
						useless_me = next(iter(scf_me))
						me_dict[useless_me.initial_state + useless_me.final_state].remove(useless_me)

					if (len(me_dict[key]) > 1) and (precedent in linked_me['main_me'].values()):
						precedent = next(iter_key)
					me_dict[key].remove(precedent)

					if key == '00' or key == '10':
						transition_text = 'r0w1'
						transition = MarchElement(transition_text)
						transition.address_order = precedent.address_order
						transition.transition_tag = True
						precedent_list.extend([transition])

					precedent_list.extend([precedent])

					return initial_based_assign(linked_me, scf_me, me_dict, precedent_list, precedent)
			return

		case _:
			return


def element_assigner(linked_CFds_me, nonCFds_me, scf_me):
	me_dict = {'00': [], '01': [], '10': [], '11': []}
	precedent_list = []
	assign_start_me = MarchElement('')
	substitute_me = MarchElement('')

	# the visit priority of MEs in the same state_key is unlinked_ME (because of no transition) > linked_ME > scf_ME(because it has candidates)
	for me in nonCFds_me.values():
		if len(me.content) > 0:
			state = me.initial_state + me.final_state
			me.address_order = 'up'
			me_dict[state].append(me)

	for me in linked_CFds_me['main_me'].values():
		if len(me.content) > 0:
			state = me.initial_state + me.final_state
			me.address_order = 'up'
			me_dict[state].append(me)

	if isinstance(linked_CFds_me['ass_me']['odd_sensitization_me'].tied_element, list):
		# if the odd sensitization element not exist, it could be under specific main element order. In this case, the
		# tied_element still need to use.
		tied_elements = linked_CFds_me['ass_me']['odd_sensitization_me'].tied_element

		if len(linked_CFds_me['ass_me']['odd_sensitization_me'].content) > 0:
			state = linked_CFds_me['ass_me']['odd_sensitization_me'].initial_state + linked_CFds_me['ass_me']['odd_sensitization_me'].final_state
			linked_CFds_me['ass_me']['odd_sensitization_me'].address_order = 'up'
			me_dict[state].append(linked_CFds_me['ass_me']['odd_sensitization_me'])

			for tied_element in tied_elements:
				if tied_element in linked_CFds_me['main_me'].values():
					tied_state = tied_element.initial_state + tied_element.final_state
					# remove the main element at other place, the element will be added with the odd sensitization ME
					me_dict[tied_state].remove(tied_element)
				else:
					tied_element.address_order = linked_CFds_me['ass_me']['odd_sensitization_me'].address_order
		else:
			# in this case, all ME in tied_elements are main MEs, so no need to add substitute_me into me_dict again
			substitute_me = tied_elements[0]
			substitute_me.tied_element = tied_elements[1:]
			for tied_element in tied_elements[1:]:
				tied_state = tied_element.initial_state + tied_element.final_state
				me_dict[tied_state].remove(tied_element)

	if len(linked_CFds_me['ass_me']['head_cover_me'].content) > 0:
		linked_CFds_me['ass_me']['head_cover_me'].address_order = 'up'
		state = linked_CFds_me['ass_me']['head_cover_me'].initial_state + linked_CFds_me['ass_me']['head_cover_me'].final_state
		me_dict[state].append(linked_CFds_me['ass_me']['head_cover_me'])
		tied_element = next(iter(linked_CFds_me['ass_me']['head_cover_me'].tied_element))
		# if the tied element is the transition ME, the transition ME not exist in me_dict,
		# so continue to check the tied element of the transition ME
		if tied_element.transition_tag:
			tied_element.address_order = linked_CFds_me['ass_me']['head_cover_me'].address_order
			sub_tied_element = next(iter(tied_element.tied_element))
			sub_tied_state = sub_tied_element.initial_state + sub_tied_element.final_state
			if sub_tied_element in me_dict[sub_tied_state]:
				me_dict[sub_tied_state].remove(sub_tied_element)
		# if not, just check the tied element itself
		else:
			tied_state = tied_element.initial_state + tied_element.final_state
			if tied_element in me_dict[tied_state]:
				me_dict[tied_state].remove(tied_element)

	for me in scf_me:
		# there are 2 candidate MEs in scf_me, just use one of them
		if len(me.content) > 0:
			state = me.initial_state + me.final_state
			me.address_order = 'up'
			me_dict[state].append(me)

	# since the tail ME has to be added separately (it has opposite AO), use this ME as the start_me. The build order is from
	# the bottom to the top
	if len(linked_CFds_me['ass_me']['tail_cover_me'].content) > 0:
		assign_start_me = linked_CFds_me['ass_me']['tail_cover_me']
		assign_start_me.address_order = 'down'
		for element in assign_start_me.tied_element:
			element.address_order = assign_start_me.address_order
		initial_based_assign(linked_CFds_me, scf_me, me_dict, precedent_list, assign_start_me)
	elif len(linked_CFds_me['ass_me']['head_cover_me'].content) > 0:
		assign_start_me = linked_CFds_me['ass_me']['head_cover_me']
		me_dict[assign_start_me.initial_state + assign_start_me.final_state].remove(assign_start_me)
		initial_based_assign(linked_CFds_me, scf_me, me_dict, precedent_list, assign_start_me)
	elif isinstance(linked_CFds_me['ass_me']['odd_sensitization_me'].tied_element, list):
		if len(linked_CFds_me['ass_me']['odd_sensitization_me'].content) > 0:
			assign_start_me = linked_CFds_me['ass_me']['odd_sensitization_me']
		else:
			assign_start_me = substitute_me
		# if odd-sensitization ME is chosen to start, it should be removed from me_dict first, since it is added before
		me_dict[assign_start_me.initial_state + assign_start_me.final_state].remove(assign_start_me)
		initial_based_assign(linked_CFds_me, scf_me, me_dict, precedent_list, assign_start_me)
	else:
		for state_key in ['00', '11', '01', '10']:
			if len(me_dict[state_key]) > 0:
				assign_start_me = next(iter(me_dict[state_key]))
				me_dict[state_key].remove(assign_start_me)
				if assign_start_me in scf_me:
					# if a scf_me candidate is chosen as the start_me, directly discard the other candidate
					useless_me = next(iter(scf_me - {assign_start_me}))
					me_dict[useless_me.initial_state + useless_me.final_state].remove(useless_me)
				break

		initial_based_assign(linked_CFds_me, scf_me, me_dict, precedent_list, assign_start_me)

	precedent_list.reverse()

	# only the CFds-detected ME, like main ME and odd-sensitization MEs, need to add a single read operation before the
	# tail-cover ME
	if (len(precedent_list) > 0) and (assign_start_me.address_order != precedent_list[0].address_order):
		try:
			insert_location = list(map(lambda o: o.address_order, precedent_list)).index(assign_start_me.address_order)
		except ValueError:
			# if no 'down' ME in precedent list, it means that the assign_start_me is the only 1 tail-cover ME, append it
			insert_location = -1

		check_location = -1 if insert_location == -1 else insert_location - 1
		if (precedent_list[check_location] not in scf_me) and (precedent_list[check_location].transition_tag is False):
			address_order_me = MarchElement('r' + precedent_list[-1].content[-1])
			address_order_me.address_order = precedent_list[0].address_order
			if check_location == -1:
				precedent_list.append(address_order_me)
			else:
				precedent_list.insert(insert_location, address_order_me)

	precedent_list += [assign_start_me]

	initial_me = MarchElement('w' + precedent_list[0].content[1])
	initial_me.address_order = 'any'
	final_me = MarchElement('r' + precedent_list[-1].content[-1])
	final_me.address_order = precedent_list[-1].address_order

	return [initial_me] + precedent_list + [final_me]


def output(march_test):
	output_file_name = 'results/march.m2c'
	text_list = []
	with open(output_file_name, 'w') as out:
		for me in march_test:
			me_text = me.address_order + ',' + me.content.replace('0', '0,').replace('1', '1,')[:-1]
			print(me_text)
			out.write(me_text + '\n')
			text_list.append(me_text)

	return text_list


if __name__ == '__main__':
	os.chdir("../")
	sys.path.append("src/")
	parsed_pool = parse_fault_pool(fault_list_file, default_fault_model_name)
	classified_pool = classify(parsed_pool)
	degenerated_2cFs = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']['unlinked'])
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], degenerated_2cFs)
	flat_SF_pool = flatten_sf_pool(filtered_SF_pool)

	seq_pool = create_sequence_pool(flat_SF_pool, degenerated_2cFs, classified_pool['2cF_CFds']['linked'], classified_pool['2cF_nonCFds_included']['nonCFds_nonCFds'])

	ME_dict = {'linked_CFds_ME': {'main_me': {'01_me': MarchElement(''), '10_me': MarchElement('')}, 'ass_me':
		{'tail_cover_me': MarchElement(''), 'odd_sensitization_me': MarchElement('')}}, 'nonCFds_ME':
				   {'00_me': MarchElement(''), '11_me': MarchElement('')}, 'scf_ME': MarchElement('')}

	if len(seq_pool['linked']['Init_0']) + len(seq_pool['linked']['Init_1']) > 0:
		ME_dict['linked_CFds_ME'] = linked_CFds_constructor(seq_pool['linked'], classified_pool['2cF_CFds']['linked'])

	if len(seq_pool['undetermined_faults']) + sum(map(lambda p: len(seq_pool['sf_seq'][p]), seq_pool['sf_seq'].keys())) + sum(map(lambda p: len(seq_pool['degenerated_seq'][p]), seq_pool['degenerated_seq'].keys())) > 0:
		nonCFds_result = nonCFds_constructor(seq_pool['degenerated_seq'], seq_pool['undetermined_faults'], seq_pool['sf_seq'])
	if len(nonCFds_result[1]) > 0:
		scf_constructor(seq_pool['unlinked']['Init_-1'])

	output(element_assigner(ME_dict['linked_CFds_ME'], ME_dict['nonCFds_ME'], ME_dict['scf_ME']))

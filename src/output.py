# assign the obtained MEs and output them as standard format
from unlinked_constructor import *


def initial_based_assign(linked_me, sf_me, me_dict, precedent_list, element):
	# if current element is odd-sensitization element, it will be assigned with its two tied main elements in fixed order
	if isinstance(element.tied_element, list):
		precedent = element.tied_element[-1]
		precedent_list.extend(element.tied_element)
		return initial_based_assign(linked_me, sf_me, me_dict, precedent_list, precedent)

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
					# if an ME for SFs are visited, the other sf_me candidate can be discarded
					if precedent in sf_me:
						useless_me = next(iter(sf_me - {precedent}))
						me_dict[useless_me.initial_state + useless_me.final_state].remove(useless_me)

					me_dict[key].remove(precedent)

					if key == '11' or key == '01':
						transition_text = 'r1w0'
						transition = MarchElement(transition_text)
						transition.address_order = precedent.address_order
						precedent_list.extend([transition])

					precedent_list.extend([precedent])

					return initial_based_assign(linked_me, sf_me, me_dict, precedent_list, precedent)
			return

		case '1':
			feature_keys = ['11', '01', '00', '10']
			for key in feature_keys:
				iter_key = iter(me_dict[key])
				precedent = next(iter_key, '')
				if not isinstance(precedent, MarchElement):
					continue
				else:
					# if an ME for SFs are visited, the other sf_me candidate can be discarded
					if precedent in sf_me:
						sf_me -= {precedent}
						useless_me = next(iter(sf_me))
						me_dict[useless_me.initial_state + useless_me.final_state].remove(useless_me)

					if (len(me_dict[key]) > 1) and (precedent in linked_me['main_me'].values()):
						precedent = next(iter_key)
					me_dict[key].remove(precedent)

					if key == '00' or key == '10':
						transition_text = 'r0w1'
						transition = MarchElement(transition_text)
						transition.address_order = precedent.address_order
						precedent_list.extend([transition])

					precedent_list.extend([precedent])

					return initial_based_assign(linked_me, sf_me, me_dict, precedent_list, precedent)
			return

		case _:
			return


def element_assigner(linked_me, unlinked_2cF_me, sf_me):
	me_dict = {'00': [], '01': [], '10': [], '11': []}
	precedent_list = []
	assign_start_me = MarchElement('')

	# the visit priority of MEs in the same state_key is unlinked_ME (because of no transition) > linked_ME > sf_ME(because it has candidates)
	for me in unlinked_2cF_me.values():
		if len(me.content) > 0:
			state = me.initial_state + me.final_state
			me.address_order = 'up'
			me_dict[state].append(me)

	for me in linked_me['main_me'].values():
		if len(me.content) > 0:
			state = me.initial_state + me.final_state
			me.address_order = 'up'
			me_dict[state].append(me)

	if len(linked_me['ass_me']['odd_sensitization_me'].content) > 0:
		if isinstance(linked_me['ass_me']['odd_sensitization_me'].tied_element, list):
			tied_elements = linked_me['ass_me']['odd_sensitization_me'].tied_element
			for tied_element in tied_elements:
				tied_state = tied_element.initial_state + tied_element.final_state
				# remove the main element at other place, the element will be added with the odd sensitization ME
				me_dict[tied_state].remove(tied_element)

		state = linked_me['ass_me']['odd_sensitization_me'].initial_state + linked_me['ass_me']['odd_sensitization_me'].final_state
		linked_me['ass_me']['odd_sensitization_me'].address_order = 'up'
		me_dict[state].append(linked_me['ass_me']['odd_sensitization_me'])

	for me in sf_me:
		# there are 2 candidate MEs in sf_me, just use one of them
		if len(me.content) > 0:
			state = me.initial_state + me.final_state
			me.address_order = 'up'
			me_dict[state].append(me)

	# since the tail ME has to be added separately (it has opposite AO), use this ME as the start_me. The build order is from
	# bottom to the top
	if len(linked_me['ass_me']['tail_cover_me'].content) > 0:
		assign_start_me = linked_me['ass_me']['tail_cover_me']
		assign_start_me.address_order = 'down'
		initial_based_assign(linked_me, sf_me, me_dict, precedent_list, assign_start_me)
	elif len(linked_me['ass_me']['odd_sensitization_me'].content) > 0:
		assign_start_me = linked_me['ass_me']['odd_sensitization_me']
		assign_start_me.address_order = 'up'
		# if odd-sensitization ME is chosen to start, it should be removed from me_dict first, since it is added before
		me_dict[assign_start_me.initial_state + assign_start_me.final_state].remove(assign_start_me)
		initial_based_assign(linked_me, sf_me, me_dict, precedent_list, assign_start_me)
	else:
		for state_key in ['00', '11', '01', '10']:
			if len(me_dict[state_key]) > 0:
				assign_start_me = next(iter(me_dict[state_key]))
				me_dict[state_key].remove(assign_start_me)
				if assign_start_me in sf_me:
					# if a sf_me candidate is chosen as the start_me, directly discard the other candidate
					useless_me = next(iter(sf_me - {assign_start_me}))
					me_dict[useless_me.initial_state + useless_me.final_state].remove(useless_me)
				break

		initial_based_assign(linked_me, sf_me, me_dict, precedent_list, assign_start_me)

	precedent_list.reverse()

	if (len(precedent_list) > 0) and (precedent_list[-1] in linked_me['main_me'].values()) and (assign_start_me is linked_me['ass_me']['tail_cover_me']):
		address_order_me = MarchElement('r' + precedent_list[-1].content[-1])
		address_order_me.address_order = precedent_list[-1].address_order
		precedent_list.append(address_order_me)

	precedent_list += [assign_start_me]

	initial_me = MarchElement('w' + precedent_list[0].content[1])
	initial_me.address_order = precedent_list[0].address_order
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
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_2cF_pool = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'],
											 classified_pool['2cF_CFds']['unlinked'])
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], filtered_2cF_pool)
	flat_SF_pool = flatten_sf_pool(filtered_SF_pool)

	seq_pool = create_sequence_pool(flat_SF_pool, filtered_2cF_pool, classified_pool['2cF_CFds']['linked'])

	ME_dict = {'linked_ME': {'main_me': {'01_me': MarchElement(''), '10_me': MarchElement('')}, 'ass_me':
		{'tail_cover_me': MarchElement(''), 'odd_sensitization_me': MarchElement('')}}, 'unlinked_2cF_ME':
		{'00_me': MarchElement(''), '11_me': MarchElement('')}, 'sf_ME': MarchElement('')}

	if len(seq_pool['linked']['Init_0']) + len(seq_pool['linked']['Init_1']) > 0:
		ME_dict['linked_ME'] = linked_CFds_constructor(seq_pool['linked'])

	if len(seq_pool['unlinked']['Init_0']) + len(seq_pool['unlinked']['Init_1']) > 0:
		ME_dict['unlinked_2cF_ME'] = unlinked_2cF_constructor(seq_pool['unlinked'])
	elif len(seq_pool['unlinked']['Init_-1']) > 0:
		ME_dict['sf_ME'] = sf_constructor(seq_pool['unlinked'])

	output(element_assigner(ME_dict['linked_ME'], ME_dict['unlinked_2cF_ME'], ME_dict['sf_ME']))

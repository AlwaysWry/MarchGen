# assign the obtained MEs and output them as standard format
from unlinked_constructor import *


def element_assigner(linked_me, unlinked_2cF_me, sf_me):
	me_dict = {'00': set(), '01': set(), '10': set(), '11': set()}
	precedent_list = []
	assign_start_me = MarchElement('')

	for me in linked_me['main_me'].values():
		if len(me.content) > 0:
			state = me.initial_state + me.final_state
			me.address_order = 'up'
			me_dict[state].add(me)

	for me in unlinked_2cF_me.values():
		if len(me.content) > 0:
			state = me.initial_state + me.final_state
			me.address_order = 'up'
			me_dict[state].add(me)

	if len(sf_me.content) > 0:
		state = sf_me.initial_state + sf_me.final_state
		sf_me.address_order = 'up'
		me_dict[state].add(sf_me)

	if len(linked_me['ass_me']['odd_sensitization_me'].content) > 0:
		if isinstance(linked_me['ass_me']['odd_sensitization_me'].tied_element, MarchElement):
			tied_element = linked_me['ass_me']['odd_sensitization_me'].tied_element
			tied_state = tied_element.initial_state + tied_element.final_state
			# remove the main element at other place, the element will be added with the odd sensitization ME
			me_dict[tied_state].remove(tied_element)

		state = linked_me['ass_me']['odd_sensitization_me'].initial_state + linked_me['ass_me']['odd_sensitization_me'].final_state
		linked_me['ass_me']['odd_sensitization_me'].address_order = 'up'
		me_dict[state].add(linked_me['ass_me']['odd_sensitization_me'])

	# since the tail ME has to be added separately (it has opposite AO), use this ME as the beginning
	def initial_based_assign(element):
		if isinstance(element.tied_element, MarchElement):
			precedent = element.tied_element
			precedent_list.extend([precedent])
			return initial_based_assign(precedent)
		match element.initial_state:
			case '0':
				iter_00 = iter(me_dict['00'])
				precedent = next(iter_00, '')
				if not isinstance(precedent, MarchElement):
					iter_10 = iter(me_dict['10'])
					precedent = next(iter_10, '')
					if not isinstance(precedent, MarchElement):
						iter_11 = iter(me_dict['11'])
						precedent = next(iter_11, '')
						if not isinstance(precedent, MarchElement):
							iter_01 = iter(me_dict['01'])
							precedent = next(iter_01, '')
							if not isinstance(precedent, MarchElement):
								return
							else:
								transition_text = 'r1w0'
								transition = MarchElement(transition_text)
								transition.address_order = precedent.address_order
								if (len(me_dict['01']) > 1) and (precedent in linked_me['main_me'].values()):
									precedent = next(iter_01)
								me_dict['01'].remove(precedent)
								precedent_list.extend([transition, precedent])
								return initial_based_assign(precedent)
						else:
							transition_text = 'r1w0'
							transition = MarchElement(transition_text)
							transition.address_order = precedent.address_order
							if (len(me_dict['11']) > 1) and (precedent in linked_me['main_me'].values()):
								precedent = next(iter_11)
							me_dict['11'].remove(precedent)
							precedent_list.extend([transition, precedent])
							return initial_based_assign(precedent)
					else:
						if (len(me_dict['10']) > 1) and (precedent in linked_me['main_me'].values()):
							precedent = next(iter_10)
						me_dict['10'].remove(precedent)
						precedent_list.extend([precedent])
						return initial_based_assign(precedent)
				else:
					if (len(me_dict['00']) > 1) and (precedent in linked_me['main_me'].values()):
						precedent = next(iter_00)
					me_dict['00'].remove(precedent)
					precedent_list.extend([precedent])
					return initial_based_assign(precedent)
			case '1':
				iter_11 = iter(me_dict['11'])
				precedent = next(iter_11, '')
				if not isinstance(precedent, MarchElement):
					iter_01 = iter(me_dict['01'])
					precedent = next(iter_01, '')
					if not isinstance(precedent, MarchElement):
						iter_00 = iter(me_dict['00'])
						precedent = next(iter_00, '')
						if not isinstance(precedent, MarchElement):
							iter_10 = iter(me_dict['10'])
							precedent = next(iter_10, '')
							if not isinstance(precedent, MarchElement):
								return
							else:
								transition_text = 'r0w1'
								transition = MarchElement(transition_text)
								transition.address_order = precedent.address_order
								if (len(me_dict['10']) > 1) and (precedent in linked_me['main_me'].values()):
									precedent = next(iter_10)
								me_dict['10'].remove(precedent)
								precedent_list.extend([transition, precedent])
								return initial_based_assign(precedent)
						else:
							transition_text = 'r0w1'
							transition = MarchElement(transition_text)
							transition.address_order = precedent.address_order
							if (len(me_dict['00']) > 1) and (precedent in linked_me['main_me'].values()):
								precedent = next(iter_00)
							me_dict['00'].remove(precedent)
							precedent_list.extend([transition, precedent])
							return initial_based_assign(precedent)
					else:
						if (len(me_dict['01']) > 1) and (precedent in linked_me['main_me'].values()):
							precedent = next(iter_01)
						me_dict['01'].remove(precedent)
						precedent_list.extend([precedent])
						return initial_based_assign(precedent)
				else:
					if (len(me_dict['11']) > 1) and (precedent in linked_me['main_me'].values()):
						precedent = next(iter_11)
					me_dict['11'].remove(precedent)
					precedent_list.extend([precedent])
					return initial_based_assign(precedent)
			case _:
				pass
		return

	if len(linked_me['ass_me']['tail_cover_me'].content) > 0:
		assign_start_me = linked_me['ass_me']['tail_cover_me']
		assign_start_me.address_order = 'down'
		initial_based_assign(assign_start_me)
	elif len(linked_me['ass_me']['odd_sensitization_me'].content) > 0:
		assign_start_me = linked_me['ass_me']['odd_sensitization_me']
		assign_start_me.address_order = 'up'
		initial_based_assign(assign_start_me)
	else:
		for state_key in ['00', '11', '01', '10']:
			if len(me_dict[state_key]) > 0:
				assign_start_me = next(iter(me_dict[state_key]))
				me_dict[state_key].remove(assign_start_me)
				break

		initial_based_assign(assign_start_me)

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

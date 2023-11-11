# a module of March test constructor based on the filtered fault pools
from sf_filter import *


class Sequence:
	"""structure of the basic element of March tests"""

	def __init__(self, props_dict):
		"""
		:param props_dict:
		seq_text: main content of a sequence object, use it to generate March test
		ass_init: the initial condition of the associate cell. For CFds, it is vInit, for nonCFds, it is aInit
		detect_tag: tag for detect operation of nonCFds
		nest_tag: tag for checking whether the sequence is available for nest sensitization

		"""
		self.seq_text = ''
		self.ass_init = ''
		self.dr_tag = ''
		self.detect_tag = ''
		self.nest_tag = ''
		self.__dict__.update(props_dict)


# End of class definitions


def get_sequence_properties(fault_obj):
	props_dict = {'seq_text': '', 'ass_init': '', 'detect_tag': ''}
	if fault_obj.CFdsFlag:
		props_dict['seq_text'] = fault_obj.aInit + fault_obj.Sen
		props_dict['ass_init'] = fault_obj.vInit
	else:
		props_dict['seq_text'] = fault_obj.vInit + fault_obj.Sen
		if fault_obj.aInit != '-':
			props_dict['ass_init'] = fault_obj.aInit
		else:
			props_dict['ass_init'] = '-1'

	props_dict['dr_tag'] = (fault_obj.rdFlag == -1)
	props_dict['detect_tag'] = not bool(fault_obj.CFdsFlag)
	props_dict['nest_tag'] = fault_obj.nestSenFlag

	return props_dict


def flatten_sf_pool(sf_pool):
	flattened_sf_pool = set()
	for init_sub in copy.deepcopy(sf_pool).values():
		for op_num_sub in init_sub.values():
			for comp_obj in op_num_sub:
				flattened_sf_pool.add(comp_obj.comps['comp1'])

	return flattened_sf_pool


def find_inclusive_sequences(seq_obj, candidate_pool):
	for candidate in candidate_pool:
		if (seq_obj is not candidate) and (seq_obj.seq_text in candidate.seq_text):
			return REDUNDANT

	return NOT_REDUNDANT


def filter_redundant_linked_sequences(linked_seq_pool):
	redundant_seq = set()
	for init_key in linked_seq_pool.keys():
		for seq in linked_seq_pool[init_key]:
			if not seq.detect_tag:
				if find_inclusive_sequences(seq, linked_seq_pool[init_key]):
					redundant_seq.add(seq)
			else:
				detect_seq = copy.deepcopy(seq)
				detect_seq.seq_text = seq.seq_text + 'r' + seq.seq_text[-1]
				if find_inclusive_sequences(detect_seq, linked_seq_pool[init_key]):
					redundant_seq.add(seq)

		linked_seq_pool[init_key] -= redundant_seq
	return


def create_sequence_pool(sf_pool, unlinked_2cF_pool, linked_CFds_pool):
	linked_seq_pool = {'Init_0': set(), 'Init_1': set()}
	# the 'Init_-1' class is for the not redundant nonCFs, it can be transformed into either of main MEs
	unlinked_seq_pool = {'Init_0': set(), 'Init_1': set(), 'Init_-1': set()}

	unlinked_pool = sf_pool.union(unlinked_2cF_pool)

	for linked_CFds in linked_CFds_pool:
		for comp_obj in linked_CFds.comps.values():
			fault_sequence = Sequence(get_sequence_properties(comp_obj))
			init_key = 'Init_' + fault_sequence.ass_init
			if isinstance(find_identical_objs(fault_sequence, linked_seq_pool[init_key], {}), type(DIFFERENT)):
				linked_seq_pool[init_key].add(copy.deepcopy(fault_sequence))

	for unlinked_fault in unlinked_pool:
		fault_sequence = Sequence(get_sequence_properties(unlinked_fault))
		init_key = 'Init_' + fault_sequence.ass_init

		# if the nonCF still exists after former filter phases, keep it into the sequence pool directly
		# if init_key == 'Init_-1':
		#	unlinked_seq_pool[init_key].add(copy.deepcopy(fault_sequence))
		#	continue

		# compare with the linked CFds pool. if the same sequence exists, merge into linked CFds pool by merging the detect_tag,
		# dr_tag and the nest_tag of nonCFds into CFds-sequence objects
		if init_key != 'Init_-1':
			find_result = find_identical_objs(fault_sequence, linked_seq_pool[init_key], {'detect_tag', 'dr_tag', 'nest_tag'})
			if not isinstance(find_result, type(DIFFERENT)):
				if not unlinked_fault.CFdsFlag:
					# only nonCFds with same sequence needs to merge the detect_tag and nest_tag
					find_result.detect_tag |= fault_sequence.detect_tag
					find_result.dr_tag |= fault_sequence.dr_tag
					find_result.nest_tag = fault_sequence.nest_tag

				continue

		# if there's no identical seq_text in linked CFds pool, compare with current unlinked sequence pool. Similarly,
		# the detect_tag, dr_tag and the nest_tag should be merged
		find_result = find_identical_objs(fault_sequence, unlinked_seq_pool[init_key], {'detect_tag', 'dr_tag', 'nest_tag'})
		if not isinstance(find_result, type(DIFFERENT)):
			if not unlinked_fault.CFdsFlag:
				find_result.detect_tag |= fault_sequence.detect_tag
				find_result.dr_tag |= fault_sequence.dr_tag
				find_result.nest_tag = fault_sequence.nest_tag
		else:
			unlinked_seq_pool[init_key].add(copy.deepcopy(fault_sequence))

	filter_redundant_linked_sequences(linked_seq_pool)

	return {'linked': linked_seq_pool, 'unlinked': unlinked_seq_pool}


if __name__ == '__main__':
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_unlinked_pool = (filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']['unlinked']))
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], filtered_unlinked_pool)

	for pool in create_sequence_pool(flatten_sf_pool(filtered_SF_pool),
									 filtered_unlinked_pool, classified_pool['2cF_CFds']['linked']).values():
		for sub_pool in pool.values():
			for sequence in sub_pool:
				print([sequence.nest_tag, sequence.seq_text])
			print("\n")

# a module of March test constructor based on the filtered fault pools
from _2cF_filter import *


class Sequence:
	"""structure of the basic element of March tests"""

	def __init__(self, props_dict):
		"""
		:param props_dict:
		seq_text: main content of a sequence object, use it to generate March test
		ass_init: the initial condition of the associate cell. For CFds, it is vInit, for nonCFds, it is aInit
		CFdsFlag: tag for CFds
		"""
		self.seq_text = ''
		self.ass_init = ''
		self.detect_tag = ''
		self.__dict__.update(props_dict)

# End of class definitions


def get_sequence_properties(fault_obj):
	props_dict = {'seq_text': '', 'ass_init': '', 'detect_tag': ''}
	if fault_obj.CFdsFlag:
		props_dict['seq_text'] = fault_obj.aInit + fault_obj.Sen
		props_dict['ass_init'] = fault_obj.vInit
	else:
		props_dict['seq_text'] = fault_obj.vInit + fault_obj.Sen
		props_dict['ass_init'] = fault_obj.aInit

	props_dict['detect_tag'] = not bool(fault_obj.CFdsFlag)

	return props_dict


def flatten_sf_pool(sf_pool):
	flattened_sf_pool = set()
	for init_sub in copy.deepcopy(sf_pool).values():
		for op_num_sub in init_sub.values():
			for comp_obj in op_num_sub:
				flattened_sf_pool.add(comp_obj.comps['comp1'])

	return flattened_sf_pool


def find_inclusive_sequences(seq, candidate_pool):
	for candidate in candidate_pool:
		if (seq is not candidate) and (seq.seq_text in candidate.seq_text):
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


def create_sequence_pool(sf_pool, unlinked_2cF_pool, linked_pool):
	linked_seq_pool = {'Init_0': set(), 'Init_1': set()}
	# the 'Init_-1' class is for the not redundant nonCFs, it can be transformed into either of main MEs
	unlinked_seq_pool = {'Init_0': set(), 'Init_1': set(), 'Init_-1': set()}

	unlinked_pool = sf_pool.union(unlinked_2cF_pool)

	for linked_CFds in linked_pool:
		for comp_obj in linked_CFds.comps.values():
			fault_sequence = Sequence(get_sequence_properties(comp_obj))
			init_key = 'Init_' + fault_sequence.ass_init
			if isinstance(find_identical_objs(fault_sequence, linked_seq_pool[init_key], {}), int):
				linked_seq_pool[init_key].add(copy.deepcopy(fault_sequence))

	for unlinked_fault in unlinked_pool:
		fault_sequence = Sequence(get_sequence_properties(unlinked_fault))
		init_key = 'Init_' + fault_sequence.ass_init

		# if the nonCF still exists after former filter phases, keep it into the sequence pool directly
		if init_key == 'Init_-1':
			unlinked_seq_pool[init_key].add(copy.deepcopy(fault_sequence))
			continue

		# compare with the CFds pool. if the same sequence exists, merge into CFds pool, but change the detect_tag to valid
		find_result = find_identical_objs(fault_sequence, linked_seq_pool[init_key], {'detect_tag'})
		if find_result == DIFFERENT:
			unlinked_seq_pool[init_key].add(copy.deepcopy(fault_sequence))
		# only nonCFds with same sequence needs to validate the detect_tag
		elif not unlinked_fault.CFdsFlag:
			find_result.detect_tag |= fault_sequence.detect_tag

	filter_redundant_linked_sequences(linked_seq_pool)

	return [linked_seq_pool, unlinked_seq_pool]


if __name__ == '__main__':
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'])
	filtered_unlinked_pool = (filter_redundant_2cF(filtered_SF_pool, classified_pool['2cF_nonCFds_included'],
												   classified_pool['2cF_CFds']))

	for pool in create_sequence_pool(flatten_sf_pool(filtered_SF_pool),
									 filtered_unlinked_pool, classified_pool['2cF_CFds']['linked']):
		for sub_pool in pool.values():
			for sequence in sub_pool:
				print(sequence.seq_text)
			print("\n")

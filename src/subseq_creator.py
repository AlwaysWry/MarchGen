# a module of March test constructor based on the filtered fault pools
from sf_filter import *


class Sequence:
	"""structure of the basic element of March tests"""

	def __init__(self, props_dict):
		"""
		:param props_dict:
		seq_text: main content of a sequence object, use it to generate March test
		ass_init: the initial condition of the associate cell. For CFds, it is vInit, for nonCFds, it is aInit
		dr_tag: tag for checking whether the sequence belongs to a CFdr fault
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


def find_inclusive_seq(target_seq, candidate_seq_pool, ignored_keys):
	if target_seq.detect_tag:
		detect_seq_text = target_seq.seq_text + 'r' + target_seq.seq_text[-1]
	else:
		detect_seq_text = target_seq.seq_text

	candidate_seqs = set(filter(lambda s: (detect_seq_text in s.seq_text) or (s.seq_text == target_seq.seq_text), candidate_seq_pool))
	return find_identical_objs(target_seq, candidate_seqs, ignored_keys)


def merge_undetermined_2cFs(_2cF_pool, _2cF_cover, linked_seq_pool, degenerated_seq_pool):
	ignored_keys = {'detect_tag', 'dr_tag', 'nest_tag'}
	redundant_undetermined_2cFs = set()
	redundant_undetermined_cover = set()
	# for the nonCFds*nonCFds, the inclusive rule is not allowed, since the included sensitization sequence of the nonCFds
	# fault cannot be ensured without overlapping by diff-based search. Since diff-based search only includes the sequences
	# of linked CFds*CFds faults or the faults in degeneration set. The difference-based build method just for making sure that
	# the merged nonCFds*nonCFds in the linked ME and unlinked ME avoid overlapping. For those nonCFds*nonCFds faults that cannot
	# be merged, still need other methods to detect.
	for nonCFds_2cF in _2cF_pool:
		find_result = {'linked_CFds': [], 'degenerated': []}
		for comp_obj in nonCFds_2cF.comps.values():
			target_seq = Sequence(get_sequence_properties(comp_obj))
			init_key = 'Init_' + target_seq.ass_init
			find_result['linked_CFds'].append(find_identical_objs(target_seq, linked_seq_pool[init_key], ignored_keys))
			find_result['degenerated'].append(find_identical_objs(target_seq, degenerated_seq_pool[init_key], ignored_keys))
		# only when the 2 FPs of the 2cF in the linked_seq_pool or the degeneration_seq_pool at the same time, i.e. no FP
		# is discarded, the merging is successful.
		if DIFFERENT not in find_result['linked_CFds'] or DIFFERENT not in find_result['degenerated']:
			redundant_undetermined_2cFs.add(nonCFds_2cF)

	# remove the merged 2cFs from the fault pool, and remove the corresponding self-filtered composite from the fault cover set.
	# Since once the
	_2cF_pool -= redundant_undetermined_2cFs
	for redundant_2cF in redundant_undetermined_2cFs:
		for comp_obj in redundant_2cF.comps.values():
			find_result = find_identical_objs(comp_obj, _2cF_cover, set())
			if not isinstance(find_result, type(DIFFERENT)):
				redundant_undetermined_cover.add(find_result)
	_2cF_cover -= redundant_undetermined_cover
	return


def filter_redundant_degenerated_sequences(composite, target_seq, candidate_seq_pool, ignored_properties):
	find_result = DIFFERENT
	seq_init_key = 'Init_' + target_seq.ass_init
	for pool_init_key in candidate_seq_pool.keys():
		# only for nonCF, all pool_init_key should be visited
		if (seq_init_key != 'Init_-1') and (pool_init_key != seq_init_key):
			continue
		# the seq_text is unnecessary to be identical, since the inclusion/read-as-root rules will be applied
		find_result = find_inclusive_seq(target_seq, candidate_seq_pool[pool_init_key], ignored_properties.union({'seq_text'}))
		if not isinstance(find_result, type(DIFFERENT)):
			break

	if not isinstance(find_result, type(DIFFERENT)):
		if (find_result.seq_text == target_seq.seq_text) and (not composite.CFdsFlag):
			# only nonCFds with same sequence needs to merge the detect_tag and nest_tag
			find_result.detect_tag |= target_seq.detect_tag
			find_result.dr_tag |= target_seq.dr_tag
			find_result.nest_tag = target_seq.nest_tag

		return REDUNDANT
	else:
		return NOT_REDUNDANT


def filter_redundant_linked_CFds_sequences(linked_seq_pool):
	def find_inclusive_sequences(seq_obj, candidate_pool):
		for candidate in candidate_pool:
			if (seq_obj is not candidate) and (seq_obj.seq_text in candidate.seq_text):
				return REDUNDANT

		return NOT_REDUNDANT

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


def create_sequence_pool(sf_pool, degenerated_2cF_pool, linked_CFds_pool, undetermined_2cF_pool, undetermined_2cF_cover ):
	linked_seq_pool = {'Init_0': set(), 'Init_1': set()}
	# the 'Init_-1' class is for the not redundant nonCFs, it can be transformed into either of main MEs
	degenerated_seq_pool = {'Init_0': set(), 'Init_1': set(), 'Init_-1': set()}
	degenerated_CF_pool = set(filter(lambda f: f.aInit != '-', degenerated_2cF_pool))
	degenerated_nonCF_pool = degenerated_2cF_pool - degenerated_CF_pool

	# generate sequence objects for linked CFds*CFds faults
	for linked_CFds in linked_CFds_pool:
		for comp_obj in linked_CFds.comps.values():
			target_sequence = Sequence(get_sequence_properties(comp_obj))
			init_key = 'Init_' + target_sequence.ass_init
			if isinstance(find_identical_objs(target_sequence, linked_seq_pool[init_key], {}), type(DIFFERENT)):
				linked_seq_pool[init_key].add(copy.deepcopy(target_sequence))

	for comp_obj in degenerated_CF_pool:
		target_sequence = Sequence(get_sequence_properties(comp_obj))
		init_key = 'Init_' + target_sequence.ass_init

		# compare with the linked CFds pool. if the same sequence or the sequence that includes the composite's sequence
		# exists, merge into linked CFds pool by merging the detect_tag, dr_tag and the nest_tag of nonCFds into CFds-sequence objects
		if filter_redundant_degenerated_sequences(comp_obj, target_sequence, linked_seq_pool,
												  {'detect_tag', 'dr_tag', 'nest_tag'}):
			continue

		# if the filter failed in linked CFds pool, compare with current unlinked sequence pool. Similarly,
		# the detect_tag, dr_tag and the nest_tag should be merged
		if not filter_redundant_degenerated_sequences(comp_obj, target_sequence, degenerated_seq_pool,
													  {'detect_tag', 'dr_tag', 'nest_tag'}):
			degenerated_seq_pool[init_key].add(copy.deepcopy(target_sequence))

	# check the redundancy of nonCFs after the check of CFs
	for comp_obj in degenerated_nonCF_pool:
		target_sequence = Sequence(get_sequence_properties(comp_obj))
		init_key = 'Init_' + target_sequence.ass_init

		if filter_redundant_degenerated_sequences(comp_obj, target_sequence, linked_seq_pool,
												  {'ass_init', 'detect_tag', 'dr_tag', 'nest_tag'}):
			continue
		if not filter_redundant_degenerated_sequences(comp_obj, target_sequence, degenerated_seq_pool,
													  {'ass_init', 'detect_tag', 'dr_tag', 'nest_tag'}):
			degenerated_seq_pool[init_key].add(copy.deepcopy(target_sequence))

	# check the same_init nonCFds*nonCFds faults

	filter_redundant_linked_CFds_sequences(linked_seq_pool)

	return {'linked': linked_seq_pool, 'unlinked': degenerated_seq_pool}


if __name__ == '__main__':
	os.chdir("../")
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filter_result = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'],
										 classified_pool['2cF_CFds']['unlinked'])
	degenerated_2cFs = filter_result[0]
	undetermined_2cFs = filter_result[1]
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], degenerated_2cFs)

	for pool in create_sequence_pool(flatten_sf_pool(filtered_SF_pool), degenerated_2cFs,
									 classified_pool['2cF_nonCFds_included']['nonCFds_nonCFds']['same_init'],
									 classified_pool['2cF_CFds']['linked']).values():
		for sub_pool in pool.values():
			for sequence in sub_pool:
				print([sequence.nest_tag, sequence.seq_text])
			print("\n")

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

	candidate_seqs = set(
		filter(lambda s: (detect_seq_text in s.seq_text) or (s.seq_text == target_seq.seq_text), candidate_seq_pool))
	return find_identical_objs(target_seq, candidate_seqs, ignored_keys)


def merge_undetermined_2cFs(_2cF_pool, linked_seq_pool, degenerated_seq_pool):
	linked_ignored_keys = {'ass_init', 'detect_tag', 'dr_tag', 'nest_tag'}
	degenerated_ignored_keys = {'detect_tag', 'dr_tag', 'nest_tag'}
	redundant_undetermined_2cFs = set()
	# for the nonCFds*nonCFds, the inclusive/read-as-root rules are not allowed, since the included sensitization sequence of the nonCFds
	# fault cannot be ensured without overlapping by diff-based search. Since diff-based search only includes the sequences
	# of linked CFds*CFds faults or the faults in degeneration set. The difference-based build method just for making sure that
	# the merged nonCFds*nonCFds in the linked ME and unlinked ME avoid overlapping. For those nonCFds*nonCFds faults that cannot
	# be merged, still need other methods to detect.
	for nonCFds_2cF in _2cF_pool:
		find_result = {'linked_CFds': [], 'degenerated': []}
		seq_temp = []
		for comp_obj in nonCFds_2cF.comps.values():
			degenerated_find_result = {'Init_0': DIFFERENT, 'Init_1': DIFFERENT, 'Init_-1': DIFFERENT}
			target_seq = Sequence(get_sequence_properties(comp_obj))
			seq_temp.append(target_seq)
			# for linked CFds sequences, Init_0 and Init_1 will be merged before the main ME is generated, since all sequences
			# in Init_0 and Init_1 need to be covered in the same middle part. As a result, the nonCFds can be merged without caring about ass_init
			linked_result = DIFFERENT
			for linked_init_key in linked_seq_pool.keys():
				linked_result = find_identical_objs(target_seq, linked_seq_pool[linked_init_key], linked_ignored_keys)
				if not isinstance(linked_result, type(DIFFERENT)):
					break
			find_result['linked_CFds'].append(linked_result)
			# for degenerated sequences, the sequences in Init_0 and Init_1 are covered in different MEs, so the nonCFds
			# should be merged according to the ass_init
			init_key = 'Init_' + target_seq.ass_init
			if init_key == 'Init_-1':
				for pool_init_key in degenerated_seq_pool.keys():
					temp_result = find_identical_objs(target_seq, degenerated_seq_pool[pool_init_key], degenerated_ignored_keys)
					if isinstance(temp_result, Sequence):
						degenerated_find_result[pool_init_key] = temp_result
				find_result['degenerated'].append(degenerated_find_result)
			else:
				temp_result = find_identical_objs(target_seq, degenerated_seq_pool[init_key], degenerated_ignored_keys)
				degenerated_find_result[init_key] = temp_result
				find_result['degenerated'].append(degenerated_find_result)

		# check the difference in linked_seq-pool. Only when the 2 FPs of the 2cF in the linked_seq_pool at the same time,
		# i.e. no FP is discarded, the merging is successful.
		if DIFFERENT not in find_result['linked_CFds']:
			for found_seq, temp in zip(find_result['linked_CFds'], seq_temp):
				found_seq.detect_tag |= temp.detect_tag
				found_seq.dr_tag |= temp.dr_tag
				found_seq.nest_tag = temp.nest_tag

			redundant_undetermined_2cFs.add(nonCFds_2cF)
			continue
		# if the nonCFds cannot be merged into the linked_seq_pool, try degenerated_seq_pool
		aInit_comp1 = nonCFds_2cF.comps['comp1'].aInit
		aInit_comp2 = nonCFds_2cF.comps['comp2'].aInit
		find_result_init_0 = list(map(lambda r: r['Init_0'], find_result['degenerated']))
		find_result_init_1 = list(map(lambda r: r['Init_1'], find_result['degenerated']))
		# if the aInits of 2 FPs in nonCFds*nonCFds are the same (or be don't care state '-'), both FPs can be sensitized
		# in once application on v-cell, so the merging is successful unless both FPs are covered by the same degenerated ME
		if (aInit_comp1 == '-') or (aInit_comp2 == '-') or (aInit_comp1 == aInit_comp2):
			feasible_result = next(iter(filter(lambda r: DIFFERENT not in r, [find_result_init_0, find_result_init_1])), [])
			if len(feasible_result) > 0:
				for found_seq, temp in zip(feasible_result, seq_temp):
					found_seq.detect_tag |= temp.detect_tag
					found_seq.dr_tag |= temp.dr_tag
					found_seq.nest_tag = temp.nest_tag
				redundant_undetermined_2cFs.add(nonCFds_2cF)
		# if the aInits are different, since degenerated MEs are 00/11 type, the states of a-cells cannot change during
		# the ME applications, so the 2 FPs cannot be sensitized in one application. As long as one of FPs is covered by the
		# degenerated_seq_pool, the nonCFds*nonCFds can be determined as redundant.
		else:
			degenerated_result_0 = list(filter(lambda r: isinstance(r, Sequence), find_result_init_0))
			degenerated_result_1 = list(filter(lambda r: isinstance(r, Sequence), find_result_init_1))
			feasible_result = next(iter(filter(lambda r: len(r) > 0, [degenerated_result_0, degenerated_result_1])), [])
			if len(feasible_result) > 0:
				for found_seq, temp in zip(feasible_result, seq_temp):
					if isinstance(found_seq, Sequence):
						found_seq.detect_tag |= temp.detect_tag
						found_seq.dr_tag |= temp.dr_tag
						found_seq.nest_tag = temp.nest_tag
						redundant_undetermined_2cFs.add(nonCFds_2cF)
				redundant_undetermined_2cFs.add(nonCFds_2cF)

	return _2cF_pool - redundant_undetermined_2cFs


def filter_redundant_other_sequences(composite, target_seq, candidate_seq_pool, ignored_properties, merged_flag=False):
	find_result = DIFFERENT
	seq_init_key = 'Init_' + target_seq.ass_init
	for pool_init_key in candidate_seq_pool.keys():
		# merged_flag represents the filter is carried out by linked_CFds_seq_pool.
		# If merged_flag is False, then only for nonCF, all pool_init_keys need to be visited
		if (not merged_flag) and (seq_init_key != 'Init_-1') and (pool_init_key != seq_init_key):
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


def create_undetermined_sequences(undetermined_fault_pool, undetermined_seq_pool):
	def save_seq_belongings(fault_obj, init_key):
		for comp_obj in fault_obj.comps.values():
			target_sequence = Sequence(get_sequence_properties(comp_obj))
			# since the FP can be a CF or a nonCF, the dr_tag should be ignored (only CF is possible to be a CFdr)
			find_result = find_identical_objs(target_sequence, undetermined_seq_pool[init_key], {'ass_init', 'dr_tag', 'belongings'})
			if isinstance(find_result, type(DIFFERENT)):
				setattr(target_sequence, 'belongings', set())
				target_sequence.belongings.add(fault_obj)
				undetermined_seq_pool[init_key].add(target_sequence)
			else:
				# the same sequence in undetermined_seq_pool is definitely already with the 'belongings' property
				find_result.belongings.add(fault_obj)

	for remainder in undetermined_fault_pool:
		init_information = remainder.comps['comp1'].aInit + remainder.comps['comp2'].aInit
		match init_information:
			case '00' | '0-' | '-0':
				save_seq_belongings(remainder, 'Init_0')
			case '11' | '1-' | '-1':
				save_seq_belongings(remainder, 'Init_1')
			case '01' | '10':
				target_comp = next(iter(remainder.comps.values()))
				target_seq = Sequence(get_sequence_properties(target_comp))
				target_key = 'Init_' + target_seq.ass_init
				target_find_result = find_identical_objs(target_seq, undetermined_seq_pool[target_key], {'dr_tag', 'belongings'})
				if isinstance(target_find_result, type(DIFFERENT)):
					setattr(target_seq, 'belongings', set())
					target_seq.belongings.add(remainder)
					undetermined_seq_pool[target_key].add(target_seq)
				else:
					target_find_result.belongings.add(remainder)
			case '--':
				save_seq_belongings(remainder, min(('Init_0', 'Init_1'), key=lambda k: len(undetermined_seq_pool[k])))
			case _:
				pass


def create_sequence_pool(sf_pool, degenerated_2cF_pool, linked_CFds_pool, undetermined_2cF_pool):
	linked_CFds_seq_pool = {'Init_0': set(), 'Init_1': set()}
	# the 'Init_-1' class is for the not redundant nonCFs, it can be transformed into either of main MEs
	degenerated_seq_pool = {'Init_0': set(), 'Init_1': set(), 'Init_-1': set()}
	degenerated_CF_pool = set(filter(lambda f: f.aInit != '-', degenerated_2cF_pool))
	degenerated_nonCF_pool = degenerated_2cF_pool - degenerated_CF_pool

	sf_seq_pool = {'Init_0': set(), 'Init_1': set(), 'Init_-1': set()}
	sf_CF_pool = set(filter(lambda f: f.aInit != '-', sf_pool))
	sf_nonCF_pool = sf_pool - sf_CF_pool

	def create_sequences(fault_pool, merged_seq_pool, target_seq_pool, ignored_properties):
		for fault_obj in fault_pool:
			target_seq = Sequence(get_sequence_properties(fault_obj))

			# compare with the linked CFds pool. if the same sequence or the sequence that includes the composite's sequence
			# exists, merge into linked CFds pool by merging the detect_tag, dr_tag and the nest_tag of nonCFds into CFds-sequence objects.
			# Note that for the filter by linked_CFds_seq_pool, the ass_init also needs to be ignored, since the sequences in linked_CFds_seq_pool
			# will be added in both linked CFds*CFds MEs, which cover the non-CFds-included faults with both 2 possible ass_inits.
			# So no matter the target_seq is CF or nonCF, it is allowed to be merged into any of the init class.
			merged_flag = True
			if filter_redundant_other_sequences(fault_obj, target_seq, merged_seq_pool, ignored_properties.union({'ass_init'}), merged_flag):
				continue
			# if the filter failed in linked CFds pool, compare with current unlinked sequence pool. Similarly,
			# the detect_tag, dr_tag and the nest_tag should be merged, and for nonCFs, the ass_init also needs to be ignored
			if not filter_redundant_other_sequences(fault_obj, target_seq, target_seq_pool, ignored_properties):
				target_seq_pool['Init_' + target_seq.ass_init].add(copy.deepcopy(target_seq))

		return

	# generate sequence objects for linked CFds*CFds faults
	for linked_CFds in linked_CFds_pool:
		for comp_obj in linked_CFds.comps.values():
			target_sequence = Sequence(get_sequence_properties(comp_obj))
			init_key = 'Init_' + target_sequence.ass_init
			if isinstance(find_identical_objs(target_sequence, linked_CFds_seq_pool[init_key], {}), type(DIFFERENT)):
				linked_CFds_seq_pool[init_key].add(copy.deepcopy(target_sequence))

	ignored_props_CF = {'detect_tag', 'dr_tag', 'nest_tag'}
	ignored_props_nonCF = {'ass_init', 'detect_tag', 'dr_tag', 'nest_tag'}
	# create and filter redundant sequences from the degenerated and simple fault pool respectively
	create_sequences(degenerated_CF_pool, linked_CFds_seq_pool, degenerated_seq_pool, ignored_props_CF)
	create_sequences(degenerated_nonCF_pool, linked_CFds_seq_pool, degenerated_seq_pool, ignored_props_nonCF)
	create_sequences(sf_CF_pool, linked_CFds_seq_pool, sf_seq_pool, ignored_props_CF)
	create_sequences(sf_nonCF_pool, linked_CFds_seq_pool, sf_seq_pool, ignored_props_nonCF)

	# merge a part of the nonCFds*nonCFds faults
	undetermined_remainder = merge_undetermined_2cFs(undetermined_2cF_pool, linked_CFds_seq_pool, degenerated_seq_pool)

	# self-filter the linked CFds*CFds sequences finally, after all filters for degenerated and SF sequences that based
	# on the linked CFds sequences finish
	filter_redundant_linked_CFds_sequences(linked_CFds_seq_pool)

	return {'linked_CFds_seq': linked_CFds_seq_pool, 'degenerated_seq': degenerated_seq_pool, 'undetermined_faults': undetermined_remainder, 'sf_seq': sf_seq_pool}


if __name__ == '__main__':
	os.chdir("../")
	parsed_pool = parse_fault_pool(fault_list_file, default_fault_model_name)
	classified_pool = classify(parsed_pool)
	degenerated_2cFs = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']['unlinked'])
	filtered_SF_pool = filter_redundant_SF(classified_pool['SF'], degenerated_2cFs)

	for pool in create_sequence_pool(flatten_sf_pool(filtered_SF_pool), degenerated_2cFs, classified_pool['2cF_CFds']['linked'],
									 classified_pool['2cF_nonCFds_included']['nonCFds_nonCFds']).values():
		for sub_pool in pool.values():
			for sequence in sub_pool:
				print([sequence.nest_tag, sequence.seq_text])
			print("\n")

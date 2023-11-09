# sf_filter.py is a module to remove the redundant simple faults for March generation.
# for SF list, the principle of removal is op_num principle and inclusion principle
# for 2cF-nonCFds and unlinked 2cF-CFds list, use minimum weight vertices cover method to simplify the fault list
from _2cF_filter import *


def find_all_occurrences(seq):
	feature = re.compile('r', re.IGNORECASE)
	positions = []

	for match in feature.finditer(seq):
		positions.append(match.start())

	return positions


def generate_nonCFds_search_candidates(seq):
	# for check the redundancy of nonCFds, apply read-root rule to the search sequence to
	# generate the candidates for comparing
	indexes = find_all_occurrences(seq)
	candidate_list = []

	if len(indexes) > 0:
		for index in indexes:
			sub_seq = seq[:index]
			for pointer in range(0, len(sub_seq) - 2, 2):
				candidate_list.append(sub_seq[pointer:])

	return candidate_list


def generate_CFds_search_candidates(seq):
	# for check the redundancy of CFds, use the sensitization sequence as the search sequence directly
	candidate_list = [seq]
	return candidate_list


def generate_inclusive_search_set(classified_fault_pool):
	# 'CFds' key and 'nonCFds' key is for checking the redundancy of CFds and nonCFds respectively, since they have
	# different redundancy checking methods.
	set_dict = {'CFds': {}, 'nonCFds': {}}

	# inclusive rules for dynamic composites
	def dynamic_faults_inclusion(c_obj):
		if arbit_SF_CFds(c_obj):
			search_seq = c_obj.comps['comp1'].aInit + c_obj.comps['comp1'].Sen
			# The candidate set for CFds can also be generated from existing nonCFds faults, since the
			# sensitization sequences in the ME (arranged according to the aInit condition of nonCFds)
			# can always satisfy the vInit condition of CFds. The nonCFds candidate set case is vice versa.
			candidate_CFds_list = generate_CFds_search_candidates(search_seq)
			candidate_nonCFds_list = generate_nonCFds_search_candidates(search_seq)
		else:
			# the possible search sequence for nonCFds also includes a detect operation
			search_seq = c_obj.comps['comp1'].vInit + c_obj.comps['comp1'].Sen + 'r' + c_obj.comps['comp1'].Sen[-1]
			candidate_CFds_list = generate_CFds_search_candidates(search_seq)
			candidate_nonCFds_list = generate_nonCFds_search_candidates(search_seq)

		return candidate_CFds_list, candidate_nonCFds_list

	# inclusive rules for static composites
	def static_faults_inclusion(c_obj):
		if arbit_SF_CFds(c_obj):
			# a CFds with single operation cannot include any other CFds or nonCFds
			candidate_CFds_list = []
			candidate_nonCFds_list = []
		else:
			# a nonCFds with single operation can only include CFds, the only nonCFds it can include is itself
			search_seq = c_obj.comps['comp1'].vInit + c_obj.comps['comp1'].Sen + 'r' + c_obj.comps['comp1'].Sen[-1]
			candidate_CFds_list = generate_CFds_search_candidates(search_seq)
			candidate_nonCFds_list = []

		return candidate_CFds_list, candidate_nonCFds_list

	for init_key in classified_fault_pool.keys():
		for type_key in set_dict.keys():
			if set_dict[type_key].get(init_key) is None:
				set_dict[type_key][init_key] = {}

		for op_num_key in classified_fault_pool[init_key].keys():
			candidate_CFds_set = set()
			candidate_nonCFds_set = set()

			if op_num_key > '#O_1':
				for comp_obj in classified_fault_pool[init_key][op_num_key]:
					dynamic_inclusion = dynamic_faults_inclusion(comp_obj)
					candidate_CFds_set.update(dynamic_inclusion[0])
					candidate_nonCFds_set.update(dynamic_inclusion[1])
			else:
				for comp_obj in classified_fault_pool[init_key][op_num_key]:
					static_inclusion = static_faults_inclusion(comp_obj)
					candidate_CFds_set.update(static_inclusion[0])
					candidate_nonCFds_set.update(static_inclusion[1])

			set_dict['CFds'][init_key][op_num_key] = candidate_CFds_set.copy()
			set_dict['nonCFds'][init_key][op_num_key] = candidate_nonCFds_set.copy()

	return set_dict


def check_nonCFds_redundancy(fault, candidate_dict, init):
	fault_op_num = fault.SenOpsNum
	fault_op_num_key = '#O_' + str(fault_op_num)
	match_seq = fault.vInit + fault.Sen
	# for nonCFds, still need to consider nest sensitization
	nest_match_seq = ''
	if fault.nestSenFlag != 'invalid':
		if fault.nestSenFlag == 'donor':
			nest_match_seq = match_seq + fault.Sen
		elif fault.vInit == '1':
			nest_match_seq = '0' + 2 * fault.Sen
		else:
			nest_match_seq = '1' + 2 * fault.Sen

	if init == 'Init_-1':
		# 1) if the nonCFds is a nonCF, its sensitization sequence can be matched by the candidate sets
		# 	 regardless the initial conditions
		for init_key in candidate_dict.keys():
			op_num_keys = sorted(candidate_dict[init_key].keys())
			# nonCF is impossible be covered by other nonCFs in the same #O class (otherwise the nonCF will be identical with
			# the current fault), it can be possibly covered by other nonCFs with more operations
			if init_key == init:
				op_num_keys = op_num_keys[op_num_keys.index(fault_op_num_key) + 1:]
				for op_num_key in op_num_keys:
					if (match_seq in candidate_dict[init_key][op_num_key]) or (
							nest_match_seq in candidate_dict[init_key][op_num_key]):
						return REDUNDANT
			# or be covered by other CFs that contain the same sensitization sequence or more operations
			else:
				for op_num_key in op_num_keys[op_num_keys.index(fault_op_num_key):]:
					if (match_seq in candidate_dict[init_key][op_num_key]) or (
							nest_match_seq in candidate_dict[init_key][op_num_key]):
						return REDUNDANT
	else:
		# 2) if the nonCFds is a CF, check whether the faults with target sequence
		# 	 exist in >fault_op_num subclasses of corresponding init class of candidate_dict
		op_num_keys = sorted(candidate_dict[init].keys())
		op_num_keys = op_num_keys[op_num_keys.index(fault_op_num_key) + 1:]
		for op_num_key in op_num_keys:
			if (match_seq in candidate_dict[init][op_num_key]) or (nest_match_seq in candidate_dict[init][op_num_key]):
				return REDUNDANT

	return NOT_REDUNDANT


def check_CFds_redundancy(fault, candidate_dict, init):
	fault_op_num = fault.SenOpsNum
	fault_op_num_key = '#O_' + str(fault_op_num)
	op_num_keys = sorted(candidate_dict[init].keys())
	# For checking the redundancy of a CFds, apply inclusive rule to every seq in candidate set to check redundancy.
	# the following 3 terms need to be considered:

	# 1) Check whether the nonCFds that has the target sequence exists in fault_op_num class of candidate_dict,
	# 	 i.e. search the candidate includes target sequence while is longer than the target

	# 2) Check whether the target sequence in the fault_op_num-1 class of candidate_dict, since the additional 1 detect
	#    operation of nonCFds in that class may cover some CFds here. However, it requires that every nonCFds has to apply
	#    the detect operation after sensitization sequence, and cannot emit when it at the tail of an ME.

	# 3) Check whether the faults that have the target sequence exist in >fault_op_num classes of candidate_dict

	search_range = op_num_keys[op_num_keys.index(fault_op_num_key) - 1:]
	match_seq = fault.aInit + fault.Sen

	for op_num_key in search_range:
		if op_num_key == '#O_' + str(fault_op_num):
			for seq in candidate_dict[init][op_num_key]:
				if (match_seq in seq) and (len(match_seq) < len(seq)):
					return REDUNDANT
			continue

		for seq in candidate_dict[init][op_num_key]:
			if match_seq in seq:
				return REDUNDANT

	return NOT_REDUNDANT


def generate_2cF_set_dict(_2cF_pool):
	flat_2cF_pool = set()
	inner_sf_pool = {'Init_0': {}, 'Init_1': {}, 'Init_-1': {}}
	inner_sf = TwoComposite()
	for _2cF_obj in _2cF_pool:
		inner_sf.get_Comp_objects(_2cF_obj, _2cF_obj)
		flat_2cF_pool.add(copy.deepcopy(inner_sf))

	for inner_obj in flat_2cF_pool:
		classify_SF(inner_sf_pool, inner_obj)

	candidate_set_dict = generate_inclusive_search_set(inner_sf_pool)
	return candidate_set_dict


def filter_redundant_SF(classified_fault_pool, _2cF_pool):
	# the filter search of SF can only happen in the same init_sub_list, since a CF can be covered by another CF
	# only if they have the same sensitization initial conditions, while nonCFs not have such limitations
	print("Filtering redundant simple faults...\n")

	filtered_fault_pool = copy.deepcopy(classified_fault_pool)
	candidate_set_dict_SF = generate_inclusive_search_set(filtered_fault_pool)
	candidate_set_dict_2cF = generate_2cF_set_dict(_2cF_pool)
	redundant_fault_pool = {}

	for init in filtered_fault_pool.keys():
		redundant_fault_pool[init] = {}
		for op_num in filtered_fault_pool[init].keys():
			redundant_fault_pool[init][op_num] = set()
			for comp_obj in filtered_fault_pool[init][op_num]:
				if arbit_SF_CFds(comp_obj):
					sf_redundancy = check_CFds_redundancy(comp_obj.comps['comp1'], candidate_set_dict_SF['CFds'], init)
					_2cF_redundancy = check_CFds_redundancy(comp_obj.comps['comp1'], candidate_set_dict_2cF['CFds'], init)
				else:
					sf_redundancy = check_nonCFds_redundancy(comp_obj.comps['comp1'], candidate_set_dict_SF['nonCFds'], init)
					_2cF_redundancy = check_nonCFds_redundancy(comp_obj.comps['comp1'], candidate_set_dict_2cF['nonCFds'], init)

				# check if the same faults exist in 2cF pool
				identical_flag = find_identical_objs(comp_obj.comps['comp1'], _2cF_pool, {'aCell'})

				if sf_redundancy or _2cF_redundancy or isinstance(identical_flag, type(comp_obj.comps['comp1'])):
					redundant_fault_pool[init][op_num].add(comp_obj)

	for init in filtered_fault_pool.keys():
		for op_num in filtered_fault_pool[init].keys():
			filtered_fault_pool[init][op_num] = filtered_fault_pool[init][op_num] - redundant_fault_pool[init][op_num]

	print("Simple faults are filtered.\n")
	return filtered_fault_pool


if __name__ == '__main__':
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_2cF_pool = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']['unlinked'])
	filter_redundant_SF(copy.deepcopy(classified_pool['SF']), filtered_2cF_pool)

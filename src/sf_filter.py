# sf_filter.py is a module to remove the redundant simple faults for March generation.
# for SF list, the principle of removal is op_num principle and inclusion principle
# for 2cF-nonCFds and unlinked 2cF-CFds list, use minimum weight vertices cover method to simplify the fault list

from basic import fault_parser as ps
import classifier as cf
import copy
import re

REDUNDANT = True
NOT_REDUNDANT = False


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
	# for check the redundancy of CFds, use the sensitization sequence directly
	candidate_list = [seq]
	return candidate_list


def generate_fault_search_set(classified_fault_pool):
	# 'CFds' key and 'nonCFds' key is for checking the redundancy of CFds and nonCFds respectively, since they have
	# different redundancy checking methods.
	set_dict = {'CFds': {}, 'nonCFds': {}}

	for init_key in classified_fault_pool.keys():

		for type_key in set_dict.keys():
			if set_dict[type_key].get(init_key) is None:
				set_dict[type_key][init_key] = {}

		for op_num_key in classified_fault_pool[init_key].keys():
			candidate_CFds_set = set()
			candidate_nonCFds_set = set()

			if op_num_key > '#O_1':
				for comp_obj in classified_fault_pool[init_key][op_num_key]:
					if cf.arbit_SF_CFds(comp_obj):
						search_seq = comp_obj.comps['comp1'].aInit + comp_obj.comps['comp1'].Sen
						# The candidate set for CFds can also be generated from existing nonCFds faults, since the
						# sensitization sequences in the ME (arranged according to the aInit condition of nonCFds)
						# can always satisfy the vInit condition of CFds. The nonCFds candidate set case is vice versa.
						candidate_CFds_list = generate_CFds_search_candidates(search_seq)
						candidate_nonCFds_list = generate_nonCFds_search_candidates(search_seq)
					else:
						# the possible search sequence for nonCFds also includes a detect operation
						search_seq = comp_obj.comps['comp1'].vInit + comp_obj.comps['comp1'].Sen + 'r' + comp_obj.comps['comp1'].Sen[-1]
						candidate_CFds_list = generate_CFds_search_candidates(search_seq)
						candidate_nonCFds_list = generate_nonCFds_search_candidates(search_seq)

					candidate_CFds_set.update(candidate_CFds_list)
					candidate_nonCFds_set.update(candidate_nonCFds_list)

			set_dict['CFds'][init_key][op_num_key] = candidate_CFds_set.copy()
			set_dict['nonCFds'][init_key][op_num_key] = candidate_nonCFds_set.copy()

	return set_dict


def check_nonCFds_redundancy(fault, candidate_dict, init):
	fault_op_num = fault.comps['comp1'].SenOpsNum
	match_seq = fault.comps['comp1'].vInit + fault.comps['comp1'].Sen

	if init == 'Init_-1':
		# if the fault is a non-CF fault, its sensitization sequence can be matched by the candidate sets
		# regardless the initial conditions
		for init_key in candidate_dict.keys():
			op_num_keys = sorted(candidate_dict[init_key].keys())
			# nonCF cannot be removed by itself, but can be covered by other CFs with the same sensitization sequence,
			# or other nonCFs with different operation number
			if init_key == init:
				for op_num_key in op_num_keys[fault_op_num:]:
					if match_seq in candidate_dict[init_key][op_num_key]:
						return REDUNDANT
			else:
				for op_num_key in op_num_keys:
					if match_seq in candidate_dict[init_key][op_num_key]:
						return REDUNDANT
	else:
		op_num_keys = sorted(candidate_dict[init].keys())
		for op_num_key in op_num_keys[fault_op_num:]:
			if match_seq in candidate_dict[init][op_num_key]:
				return REDUNDANT

	return NOT_REDUNDANT


def check_CFds_redundancy(fault, candidate_dict, init):
	fault_op_num = fault.comps['comp1'].SenOpsNum
	op_num_keys = sorted(candidate_dict[init].keys())
	# for CFds, besides the sets that op_num larger than the SenOpsNum, the set whose op_num is 1 less than SenOpsNum
	# also needs to be checked. Since the additional 1 detect operation of nonCFds in that class may cover some CFds here
	search_range = op_num_keys[max(0, fault_op_num - 2):]
	search_range.remove('#O_' + str(fault_op_num))

	match_seq = fault.comps['comp1'].aInit + fault.comps['comp1'].Sen

	for op_num_key in search_range:
		for seq in candidate_dict[init][op_num_key]:
			if match_seq in seq:
				return REDUNDANT

	return NOT_REDUNDANT


def filter_redundant_SF(classified_fault_pool):
	# the filter search of SF can only happen in the same init_sub_list, since a CF can be covered by another CF
	# only if they have the same sensitization initial conditions, while nonCFs not have such limitations
	print("Filtering redundant simple faults...\n")

	filtered_fault_pool = copy.deepcopy(classified_fault_pool)
	candidate_set_dict = generate_fault_search_set(filtered_fault_pool)
	redundant_fault_pool = {}

	for init in filtered_fault_pool.keys():
		redundant_fault_pool[init] = {}
		for op_num in filtered_fault_pool[init].keys():
			redundant_fault_pool[init][op_num] = set()
			for comp_obj in filtered_fault_pool[init][op_num]:
				if cf.arbit_SF_CFds(comp_obj):
					if check_CFds_redundancy(comp_obj, candidate_set_dict['CFds'], init):
						redundant_fault_pool[init][op_num].add(comp_obj)
				else:
					if check_nonCFds_redundancy(comp_obj, candidate_set_dict['nonCFds'], init):
						redundant_fault_pool[init][op_num].add(comp_obj)

	for init in filtered_fault_pool.keys():
		for op_num in filtered_fault_pool[init].keys():
			filtered_fault_pool[init][op_num] = filtered_fault_pool[init][op_num] - redundant_fault_pool[init][op_num]

	return filtered_fault_pool


if __name__ == '__main__':
	parsed_pool = cf.parse_fault_pool(ps.fault_list_file, ps.fault_model_name)
	classify_result = cf.classify(parsed_pool)
	sf_pool = copy.deepcopy(classify_result['SF'])
	filter_redundant_SF(sf_pool)

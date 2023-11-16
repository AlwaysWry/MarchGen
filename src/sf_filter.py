# sf_filter.py is a module to remove the redundant simple faults for March generation.
# for SF list, the principle of removal is op_num principle and inclusion principle
# for 2cF-nonCFds and unlinked 2cF-CFds list, use minimum weight vertices cover method to simplify the fault list
from _2cF_filter import *


def generate_2cF_SF_set_dict(_2cF_pool):
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
	print("\n***Filtering redundant simple faults...\n")

	filtered_fault_pool = copy.deepcopy(classified_fault_pool)
	candidate_set_dict_SF = generate_inclusive_search_set(filtered_fault_pool)
	# faults in 2cF pool get more strict test condition than SF, so use 2cF_pool to filter the current SF pool
	candidate_set_dict_2cF = generate_2cF_SF_set_dict(_2cF_pool)
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

				# check if the same faults exist in 2cF pool, it covers the missing case in sf_filter, the 2cF in 2cF_pool
				# with the same sensitization sequence as the SF.
				identical_flag = find_identical_objs(comp_obj.comps['comp1'], _2cF_pool, {'aCell'})

				if sf_redundancy or _2cF_redundancy or isinstance(identical_flag, type(comp_obj.comps['comp1'])):
					redundant_fault_pool[init][op_num].add(comp_obj)

	for init in filtered_fault_pool.keys():
		for op_num in filtered_fault_pool[init].keys():
			filtered_fault_pool[init][op_num] = filtered_fault_pool[init][op_num] - redundant_fault_pool[init][op_num]

	# print("Simple faults are filtered.\n")
	return filtered_fault_pool


if __name__ == '__main__':
	os.chdir("../")
	parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
	classified_pool = classify(parsed_pool)
	filtered_2cF_pool = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']['unlinked'])
	filter_redundant_SF(copy.deepcopy(classified_pool['SF']), filtered_2cF_pool)

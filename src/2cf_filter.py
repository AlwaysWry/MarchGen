# 2cf_filter.py is a module of filter the redundant 2-composite nonCFds included and the unlinked CFds*CFds faults.
# The filter strategies include the sf-based method and the minimum weight vertices coverage method.
from basic import fault_parser as ps
import classifier as cf
import sf_filter as sff

IDENTICAL = True
DIFFERENT = False
REDUNDANT = True
NOT_REDUNDANT = False


def find_identical_comps(fault_obj, sub_sf_pool, ignored_keys):
	record = []
	for sf_obj in sub_sf_pool:
		for key, value in fault_obj.__dict__.items():
			if key not in ignored_keys and sf_obj.comps['comp1'].__dict__[key] != value:
				record.append(False)
				break
			else:
				record.append(True)

		if False not in record:
			return IDENTICAL
		record.clear()

	return DIFFERENT


def check_2cF_redundancy_by_SF(sf_pool, fault_obj, init):
	op_num_key = '#O_' + str(fault_obj.SenOpsNum)
	if init == '-':
		for init_key in sf_pool.keys():
			if find_identical_comps(fault_obj, sf_pool[init_key][op_num_key], {fault_obj.aInit, fault_obj.aCell}):
				return REDUNDANT
	else:
		init_key = 'Init_' + init
		if find_identical_comps(fault_obj, sf_pool[init_key][op_num_key], set()):
			return REDUNDANT

	return NOT_REDUNDANT


def remove_2cF_based_on_SF(sf_pool, _2cF_nonCFds_pool):
	redundancy = []
	redundant_2cF_pool = set()
	for _2cF_obj in _2cF_nonCFds_pool:
		for comp in _2cF_obj.comps.values():
			if comp.CFdsFlag:
				init = comp.vInit
			else:
				init = comp.aInit
			redundancy.append(check_2cF_redundancy_by_SF(sf_pool, comp, init))

		if REDUNDANT in redundancy:
			redundant_2cF_pool.add(_2cF_obj)
		redundancy.clear()

	return _2cF_nonCFds_pool - redundant_2cF_pool


def remove_2cF_based_on_MWVC():
	pass


def filter_redundant_2cF(_2cF_nonCFds_pool, _2cF_CFds_pool):
	pass


parsed_pool = cf.parse_fault_pool(ps.fault_list_file, ps.fault_model_name)
classified_pool = cf.classify(parsed_pool)
filtered_SF_pool = sff.filter_redundant_SF(classified_pool['SF'])
simplified_2cF_nonCFds_pool = remove_2cF_based_on_SF(classified_pool['SF'], classified_pool['2cF_nonCFds_included'])
for fault in simplified_2cF_nonCFds_pool:
	print(fault.fault_primitive)

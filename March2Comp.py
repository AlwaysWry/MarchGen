# main script of the March2Comp program
from src.unlinked_constructor import *


def March2Comp(fault_list, fault_model):
	# parse the input fault list
	parsed_faults = parse_fault_pool(fault_list, fault_model)
	if len(parsed_faults) == 0:
		print("Fault list is empty!\n")
		return
	# classify
	classified_faults = classify(parsed_faults)
	# filter 2cF pool
	filtered_2cFs = filter_redundant_2cF(classified_faults['2cF_nonCFds_included'], classified_faults['2cF_CFds']['unlinked'])
	# filter SF pool
	filtered_SFs = filter_redundant_SF(classified_faults['SF'], filtered_2cFs)
	flat_SF_set = flatten_sf_pool(filtered_SFs)
	# create sequence objects
	sequence_pool = create_sequence_pool(flat_SF_set, filtered_2cFs, classified_faults['2cF_CFds']['linked'])
	# build MEs for linked CFds


	pass


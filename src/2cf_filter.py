# 2cf_filter.py is a module of filter the redundant 2-composite nonCFds included and the unlinked CFds*CFds faults.
# The filter strategies include the sf-based method and the minimum weight vertices coverage method.
from basic import fault_parser as ps
import classifier as cf
import sf_filter as sff
import sys

DIFFERENT = -1
REDUNDANT = True
NOT_REDUNDANT = False


def find_identical_comps(fault_obj, candidate_pool, ignored_keys):
	record = []
	for comp_obj in candidate_pool:
		# check the type of the input fault, for the function can be used in different situations
		if type(comp_obj).__name__ == 'SimpleFault':
			candidate = comp_obj
		else:
			candidate = comp_obj.comps['comp1']

		for key, value in fault_obj.__dict__.items():
			if key not in ignored_keys and candidate.__dict__[key] != value:
				record.append(False)
				break
			else:
				record.append(True)

		if False not in record:
			return comp_obj
		record.clear()

	return DIFFERENT


def check_2cF_redundancy_by_SF(sf_pool, fault_obj, init):
	op_num_key = '#O_' + str(fault_obj.SenOpsNum)
	if init == '-':
		for init_key in sf_pool.keys():
			if find_identical_comps(fault_obj, sf_pool[init_key][op_num_key], {'aInit', 'aCell'}) != DIFFERENT:
				return REDUNDANT
	else:
		init_key = 'Init_' + init
		if find_identical_comps(fault_obj, sf_pool[init_key][op_num_key], {'aCell'}) != DIFFERENT:
			return REDUNDANT

	return NOT_REDUNDANT


def remove_2cF_based_on_SF(sf_pool, _2cF_unlinked_pool):
	redundancy = []
	redundant_2cF_pool = set()
	for _2cF_obj in _2cF_unlinked_pool:
		for comp in _2cF_obj.comps.values():
			if comp.CFdsFlag:
				init = comp.vInit
			else:
				init = comp.aInit
			redundancy.append(check_2cF_redundancy_by_SF(sf_pool, comp, init))

		if REDUNDANT in redundancy:
			redundant_2cF_pool.add(_2cF_obj)
		redundancy.clear()

	return _2cF_unlinked_pool - redundant_2cF_pool


def get_fault_operation_num(fault_obj):
	if fault_obj.CFdsFlag:
		weight = fault_obj.SenOpsNum
	elif fault_obj.nestSenFlag == 'donor':
		# for the nonCFds while can be nest-sensitized faults, the middle detect operation can be emitted
		weight = fault_obj.SenOpsNum
	else:
		weight = fault_obj.SenOpsNum + 1

	return weight


def output_graph_file_QUICK_VC(vertices, edges):
	with open('../resources/unlinked_2cF.m2c', 'w') as graph:
		graph.write(str(len(edges)) + ' ' + str(len(vertices)) + '\n')
		for vertex in vertices[:-1]:
			weight = list(vertex.values())
			graph.write(str(weight[0]) + ' ')
		weight = list(vertices[-1].values())
		graph.write(str(weight[0]) + '\n')

		for edge in edges[:-1]:
			terminal_1 = edge[0] + 1
			terminal_2 = edge[1] + 1
			graph.write(str(terminal_1) + ' ' + str(terminal_2) + '\n')
		graph.write(str(edges[-1][0] + 1) + ' ' + str(edges[-1][1] + 1))

	return


def output_graph_file_DYNWVC2(vertices, edges):
	with open('../resources/unlinked_2cF.m2c', 'w') as graph:
		graph.write('p edge ' + str(len(vertices)) + ' ' + str(len(edges)) + '\n')
		for vertex in vertices:
			index = list(vertex.keys())
			weight = list(vertex.values())
			graph.write('v ' + str(index[0]) + ' ' + str(weight[0]) + '\n')

		for edge in edges[:-1]:
			terminal_1 = edge[0] + 1
			terminal_2 = edge[1] + 1
			graph.write('e ' + str(terminal_1) + ' ' + str(terminal_2) + '\n')
		graph.write('e ' + str(edges[-1][0] + 1) + ' ' + str(edges[-1][1] + 1))

	return


def build_unlinked_2cF_graph(_2cF_unlinked_pool):
	vertices_map = []
	vertices = []
	edges = []
	for _2cF_obj in _2cF_unlinked_pool:
		edge_info = []
		for comp in _2cF_obj.comps.values():
			if comp.aInit == '-':
				ignore_keys = {'aInit', 'aCell'}
			else:
				ignore_keys = {'aCell'}

			find_result = find_identical_comps(comp, vertices_map, ignore_keys)
			if find_result == DIFFERENT:
				vertices_map.append(comp)
				vertices.append({len(vertices_map) - 1: get_fault_operation_num(comp)})
				edge_info.append(vertices_map.index(comp))
			else:
				edge_info.append(vertices_map.index(find_result))

		edge_info.sort()
		edges.append(edge_info)
		edges.sort()

	if sys.platform.startswith('linux'):
		output_graph_file_QUICK_VC(vertices, edges)
	else:
		output_graph_file_DYNWVC2(vertices, edges)

	return


def remove_2cF_based_on_MWVC():

	pass


def filter_redundant_2cF(_2cF_nonCFds_pool, _2cF_CFds_pool):
	pass


parsed_pool = cf.parse_fault_pool(ps.fault_list_file, ps.fault_model_name)
classified_pool = cf.classify(parsed_pool)
filtered_SF_pool = sff.filter_redundant_SF(classified_pool['SF'])
# simplified_2cF_nonCFds_pool = remove_2cF_based_on_SF(classified_pool['SF'], classified_pool['2cF_nonCFds_included'] |
# 													 classified_pool['2cF_CFds']['unlinked'])
build_unlinked_2cF_graph(classified_pool['2cF_nonCFds_included'] | classified_pool['2cF_CFds']['unlinked'])

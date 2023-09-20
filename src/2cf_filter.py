# 2cf_filter.py is a module of filter the redundant 2-composite nonCFds included and the unlinked CFds*CFds faults.
# The filter strategies include the sf-based method and the minimum weight vertices coverage method.
from basic import fault_parser as ps
import classifier as cf
import sf_filter as sff
import sys

if sys.platform.startswith('linux'):
	import quickVC_solver
else:
	import dynWVC2_solver

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
	# if the simple fault pool is empty, return the original 2cF unlinked pool directly
	for init_index, init_key in enumerate(sf_pool.keys(), 1):
		if any(sf_pool[init_key]):
			break
		elif init_index < len(sf_pool.keys()):
			continue
		else:
			return _2cF_unlinked_pool

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
	filename = '../results/unlinked_2cF.m2c'
	with open(filename, 'w') as graph:
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

	return filename


def output_graph_file_DYNWVC2(vertices, edges):
	filename = '../results/unlinked_2cF.m2c'
	with open(filename, 'w') as graph:
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

	return filename


def build_unlinked_2cF_graph(_2cF_unlinked_pool, vertices_map, CFdr_map):
	vertices = []
	edges = []
	for _2cF_obj in _2cF_unlinked_pool:
		edge_info = []
		CFdr_flag = 0
		for comp in _2cF_obj.comps.values():
			# if a composite is a CFdr, it must be detected, so the 2cF can't be involved
			if (comp.rdFlag == -1) and (comp.CFdsFlag == 0):
				CFdr_map.append(comp)
				CFdr_flag = 1
				break

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

		# if a composite is a CFdr, it must be detected, so the 2cF can't be involved
		if CFdr_flag:
			continue
		edge_info.sort()
		edges.append(edge_info)
		edges.sort()

	if sys.platform.startswith('linux'):
		graph_file = output_graph_file_QUICK_VC(vertices, edges)
	else:
		graph_file = output_graph_file_DYNWVC2(vertices, edges)

	return graph_file


def remove_2cF_based_on_MWVC(graph_file):
	if sys.platform.startswith('linux'):
		quickVC_solver.quickVC_solver(graph_file)
	else:
		dynWVC2_solver.dynWVC2_solver(graph_file)

	return


def filter_redundant_2cF(sf_pool, _2cF_nonCFds_pool, _2cF_CFds_pool):
	print("Filtering redundant 2-composite faults...\n")

	unlinked_2cF_pool = _2cF_nonCFds_pool | _2cF_CFds_pool['unlinked']
	simplified_unlinked_2cF_pool = remove_2cF_based_on_SF(sf_pool, unlinked_2cF_pool)
	_2cF_cover = set()

	if len(simplified_unlinked_2cF_pool) > 0:
		vertices_map = []
		CFdr_map = []
		print("Invoking MWVC solver...\n")
		graph_file = build_unlinked_2cF_graph(simplified_unlinked_2cF_pool, vertices_map, CFdr_map)
		remove_2cF_based_on_MWVC(graph_file)

		with open("../results/mwvc.log", "r") as result:
			for vertex in result.readlines():
				vertex.strip()
				_2cF_cover.add(vertices_map[int(vertex) - 1])

		for CFdr in CFdr_map:
			if CFdr.aInit == '-':
				ignore_keys = {'aInit', 'aCell'}
			else:
				ignore_keys = {'aCell'}
			if find_identical_comps(CFdr, _2cF_cover, ignore_keys) == DIFFERENT:
				_2cF_cover.add(CFdr)

	return _2cF_cover


if __name__ == '__main__':
	try:
		parsed_pool = cf.parse_fault_pool(ps.fault_list_file, ps.fault_model_name)
		classified_pool = cf.classify(parsed_pool)
		filtered_SF_pool = sff.filter_redundant_SF(classified_pool['SF'])
		for fault in filter_redundant_2cF(classified_pool['SF'], classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']):
			print(fault.text)
	except TypeError:
		pass

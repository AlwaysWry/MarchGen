# _2cF_filter.py is a module of filter the redundant 2-composite nonCFds included and the unlinked CFds*CFds faults.
# The filter strategies include the sf-based method and the minimum weight vertices coverage method.
import traceback
from classifier import *
import sys

REDUNDANT = True
NOT_REDUNDANT = False

# apply different MWVC solver for different OS
if sys.platform.startswith('linux'):
	import quickVC_solver
else:
	import dynWVC2_solver

DIFFERENT = -1


def find_identical_objs(obj, candidate_pool, ignored_keys):
	record = []
	for comp_obj in candidate_pool:
		# check the type of the input fault, for the function can be used in different situations
		if type(comp_obj).__name__ == 'TwoComposite':
			# if the candidate is TwoComposite type, use its comp one by one
			candidate = comp_obj.comps.values()
		else:
			# if the candidate is other type, use its properties directly
			candidate = [comp_obj]

		for sub_candidate in candidate:
			for key, value in obj.__dict__.items():
				if (key not in ignored_keys) and (sub_candidate.__dict__[key] != value):
					record.append(False)
					break
				else:
					record.append(True)

			if False not in record:
				return comp_obj
			record.clear()

	return DIFFERENT


def get_vertex_weight(fault_obj):
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

			find_result = find_identical_objs(comp, vertices_map, ignore_keys)
			if find_result == DIFFERENT:
				vertices_map.append(comp)
				vertices.append({len(vertices_map) - 1: get_vertex_weight(comp)})
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


def remove_unlinked_2cF_by_MWVC(graph_file):
	if sys.platform.startswith('linux'):
		quickVC_solver.quickVC_solver(graph_file)
	else:
		dynWVC2_solver.dynWVC2_solver(graph_file)

	return


def remove_unlinked_2cF_by_linked_2cF_CFds(simplified_unlinked_2cF_cover, linked_2cF_CFds_pool):
	redundant_comp = set()
	for comp in simplified_unlinked_2cF_cover:
		if isinstance(find_identical_objs(comp, linked_2cF_CFds_pool, {'aCell'}), int):
			redundant_comp.add(comp)

	return simplified_unlinked_2cF_cover - redundant_comp


def filter_redundant_2cF(_2cF_nonCFds_pool, _2cF_CFds_pool):
	print("Filtering redundant 2-composite faults...\n")

	unlinked_2cF_pool = _2cF_nonCFds_pool | _2cF_CFds_pool
	unlinked_2cF_cover = set()

	if len(unlinked_2cF_pool) > 0:
		vertices_map = []
		CFdr_map = []
		print("Building unlinked 2cF graph...\n")
		graph_file = build_unlinked_2cF_graph(unlinked_2cF_pool, vertices_map, CFdr_map)
		print("Invoking MWVC solver...\n")
		remove_unlinked_2cF_by_MWVC(graph_file)

		with open("../results/mwvc.log", "r") as result:
			for vertex in result.readlines():
				vertex.strip()
				unlinked_2cF_cover.add(vertices_map[int(vertex) - 1])

		for CFdr in CFdr_map:
			if CFdr.aInit == '-':
				ignore_keys = {'aInit', 'aCell'}
			else:
				ignore_keys = {'aCell'}
			if isinstance(find_identical_objs(CFdr, unlinked_2cF_cover, ignore_keys), int):
				unlinked_2cF_cover.add(CFdr)

	print("\n2-composite faults are filtered.\n")
	return unlinked_2cF_cover


if __name__ == '__main__':
	try:
		parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
		classified_pool = classify(parsed_pool)
		# TODO: consider the case that 2cF pool is empty at top level
		filtered_2cF_pool = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']['unlinked'])
	except TypeError:
		print("fail")
		traceback.print_exc()

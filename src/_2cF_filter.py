# _2cF_filter.py is a module of filter the redundant 2-composite nonCFds included and the unlinked CFds*CFds faults.
# The filter strategies include the sf-based method and the minimum weight vertices coverage method.
import traceback
from classifier import *

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
	filename = 'results/unlinked_2cF_graph.m2c'
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
	filename = 'results/unlinked_2cF_graph.m2c'
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


def build_degenerate_2cF_graph(_2cF_unlinked_pool, vertices_map):
	vertices = []
	edges = []
	for _2cF_obj in _2cF_unlinked_pool:
		edge_info = []
		for comp in _2cF_obj.comps.values():
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

		edge_info.sort()
		edges.append(edge_info)

	edges.sort()

	if (len(edges) > 0) and (len(vertices) > 0):
		if sys.platform.startswith('linux'):
			graph_file = output_graph_file_QUICK_VC(vertices, edges)
		else:
			graph_file = output_graph_file_DYNWVC2(vertices, edges)
		return graph_file
	else:
		return False


def remove_degenerate_2cF_by_MWVC(graph_file):
	if sys.platform.startswith('linux'):
		quickVC_solver.quickVC_solver(graph_file)
	else:
		dynWVC2_solver.dynWVC2_solver(graph_file)

	return


def find_all_read_occurrences(seq):
	feature = re.compile('r', re.IGNORECASE)
	positions = []

	for match in feature.finditer(seq):
		positions.append(match.start())

	return positions


def generate_nonCFds_search_candidates(seq):
	# for check the redundancy of nonCFds, apply read-root rule to the search sequence to
	# generate the candidates for comparing
	indexes = find_all_read_occurrences(seq)
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
			# the possible search sequence for nonCFds should also include a detect operation, but it should be emitted
			# because the nest sensitization can eliminate the intermediate read operation, and the nest information cannot be
			# known until the ME is generated. As a result, the search sequence keeps the sensitization sequence itself.
			search_seq = c_obj.comps['comp1'].vInit + c_obj.comps['comp1'].Sen
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
			search_seq = c_obj.comps['comp1'].vInit + c_obj.comps['comp1'].Sen
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

			# skip it if candidate_dict is empty
			if len(op_num_keys) == 0:
				continue

			# nonCF is impossible be covered by other nonCFs in the same #O class (otherwise the nonCF will be identical with
			# the current fault), it can be possibly covered by other nonCFs with more operations
			if init_key == init:
				op_num_keys = sorted(filter(lambda k: k > fault_op_num_key, op_num_keys))
				for op_num_key in op_num_keys:
					if (match_seq in candidate_dict[init_key][op_num_key]) or (
							nest_match_seq in candidate_dict[init_key][op_num_key]):
						return REDUNDANT
			# or be covered by other CFs that contain the same sensitization sequence or more operations
			else:
				op_num_keys = sorted(filter(lambda k: k >= fault_op_num_key, op_num_keys))
				for op_num_key in op_num_keys:
					if (match_seq in candidate_dict[init_key][op_num_key]) or (
							nest_match_seq in candidate_dict[init_key][op_num_key]):
						return REDUNDANT
	else:
		# 2) if the nonCFds is a CF, check whether the faults with target sequence
		# 	 exist in >fault_op_num subclasses of corresponding init class of candidate_dict
		op_num_keys = sorted(candidate_dict[init].keys())

		# skip it if candidate_dict is empty
		if len(op_num_keys) == 0:
			return NOT_REDUNDANT

		op_num_keys = sorted(filter(lambda k: k > fault_op_num_key, op_num_keys))
		for op_num_key in op_num_keys:
			if (match_seq in candidate_dict[init][op_num_key]) or (nest_match_seq in candidate_dict[init][op_num_key]):
				return REDUNDANT

	return NOT_REDUNDANT


def check_CFds_redundancy(fault, candidate_dict, init):
	fault_op_num = fault.SenOpsNum
	op_num_keys = sorted(candidate_dict[init].keys())

	# skip it if candidate_dict is empty
	if len(op_num_keys) == 0:
		return NOT_REDUNDANT
	# For checking the redundancy of a CFds, apply inclusive rule to every seq in candidate set to check redundancy.
	# the following 3 terms need to be considered:

	# 1) Check whether the nonCFds that has the target sequence exists in fault_op_num class of candidate_dict,
	# 	 i.e. search the candidate includes target sequence while is longer than the target

	# 2) Check whether the target sequence in the fault_op_num-1 class of candidate_dict, since the additional 1 detect
	#    operation of nonCFds in that class may cover some CFds here. However, it requires that every nonCFds has to apply
	#    the detect operation after sensitization sequence, and cannot emit when it at the tail of an ME.

	# 3) Check whether the faults that have the target sequence exist in >fault_op_num classes of candidate_dict

	op_num_keys = sorted(filter(lambda k: k >= '#O_' + str(fault_op_num - 1), op_num_keys))
	match_seq = fault.aInit + fault.Sen

	for op_num_key in op_num_keys:
		if op_num_key == '#O_' + str(fault_op_num):
			for seq in candidate_dict[init][op_num_key]:
				if (match_seq in seq) and (len(match_seq) < len(seq)):
					return REDUNDANT
			continue

		for seq in candidate_dict[init][op_num_key]:
			if match_seq in seq:
				return REDUNDANT

	return NOT_REDUNDANT


def generate_2cF_2cF_set_dict(sf_set, _2cF_set):
	inner_sf_pool = {'Init_0': {}, 'Init_1': {}, 'Init_-1': {}}
	inner_sf = TwoComposite()
	for sf_obj in sf_set:
		inner_sf.get_Comp_objects(sf_obj, sf_obj)
		_2cF_set.add(copy.deepcopy(inner_sf))

	for inner_obj in _2cF_set:
		classify_SF(inner_sf_pool, inner_obj)

	return generate_inclusive_search_set(inner_sf_pool)


def self_filter(_2cF_cover):
	# the 2cFs can be filtered according to the inclusive rule
	redundant_fault_pool = set()
	filtered_fault_pool = set()
	candidate_set_dict = generate_2cF_2cF_set_dict(_2cF_cover, filtered_fault_pool)

	for fault in filtered_fault_pool:
		comp = fault.comps['comp1']
		if arbit_SF_CFds(fault):
			init = 'Init_' + comp.vInit
			redundancy = check_CFds_redundancy(comp, candidate_set_dict['CFds'], init)
		else:
			if comp.aInit != '-':
				init = 'Init_' + comp.aInit
			else:
				init = 'Init_-1'
			redundancy = check_nonCFds_redundancy(comp, candidate_set_dict['nonCFds'], init)

		if redundancy:
			redundant_fault_pool.add(fault)

	return set(map(lambda f: f.comps['comp1'], filtered_fault_pool - redundant_fault_pool))


def filter_redundant_2cF(_2cF_nonCFds_pool, unlinked_2cF_CFds_pool):
	print("***Filtering redundant 2-composite faults...\n")

	# all nonCFds*nonCFds faults cannot be degenerated, since the discarded nonCFds FP can be mis-sensitized if it is
	# out of the diff-based search range. For the nonCFds*nonCFds with different a-cell requirements (such as <1;...>*<0;...>),
	# both of the FPs can be sensitized under [a1, v, a2] or [a2, v, a1]. and for the faults with identical a-cell requirements,
	# the FPs can be sensitized under [a1, a2, v] etc.
	degenerate_2cF_pool = _2cF_nonCFds_pool['nonCFds_CFds'] | unlinked_2cF_CFds_pool
	degenerate_2cF_cover = set()

	if len(degenerate_2cF_pool) > 1:
		vertices_map = []
		print("Building unlinked 2cF graph...")
		graph_file = build_degenerate_2cF_graph(degenerate_2cF_pool, vertices_map)

		if isinstance(graph_file, str):
			print("Invoking MWVC solver...")
			remove_degenerate_2cF_by_MWVC(graph_file)

			with open("results/mwvclog", "r") as result:
				for vertex in result.readlines():
					vertex.strip()
					degenerate_2cF_cover.add(vertices_map[int(vertex) - 1])
	else:
		# if there is only 1 2cF, use its 2 vertices as the final cover directly
		for degenerated_2cF in degenerate_2cF_pool:
			degenerate_2cF_cover.update(set(degenerated_2cF.comps.values()))

	# filter the 2cF cover after MWVC, according to inclusive rule
	if len(degenerate_2cF_cover) > 1:
		degenerate_2cF_cover = self_filter(degenerate_2cF_cover)

	# the nonCFds*nonCFds faults can only carry out self-filter
	# TODO: really can carry out self-filter?
	undetermined_2cF_pool = _2cF_nonCFds_pool['nonCFds_nonCFds']
	undetermined_2cF_cover = set()

	for _2cF_obj in undetermined_2cF_pool:
		for comp in _2cF_obj.comps.values():
			ignore_keys = {'aCell'}

			find_result = find_identical_objs(comp, undetermined_2cF_cover, ignore_keys)
			if find_result == DIFFERENT:
				undetermined_2cF_cover.add(comp)

	if len(undetermined_2cF_cover) > 1:
		undetermined_2cF_cover = self_filter(undetermined_2cF_cover)

	return degenerate_2cF_cover, undetermined_2cF_cover


if __name__ == '__main__':
	os.chdir("../")
	try:
		parsed_pool = parse_fault_pool(fault_list_file, fault_model_name)
		classified_pool = classify(parsed_pool)
		filtered_2cF_pool = filter_redundant_2cF(classified_pool['2cF_nonCFds_included'], classified_pool['2cF_CFds']['unlinked'])[0]

	except TypeError:
		print("fail")
		traceback.print_exc()

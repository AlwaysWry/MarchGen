# A parse and classify module of input simple faults and 2cFs
import copy
import sys

sys.path.append("../")

from basic.fault_parser import *

# MACRO definitions
CLASSIFY_ERROR = -1
IS_SF = True
NOT_SF = False
IS_SF_CFDS = True
NOT_SF_CFDS = False
IS_NONCFDS_INCLUDED = True
NOT_NONCFDS_INCLUDED = False


class TwoComposite:
	"""a 2-composite fault contains two simple faults"""
	link_flag = ''
	fp_text = ''

	def __init__(self):
		self.comps = {}

	def get_FP_text(self, text):
		self.fp_text = text

	def get_Comp_objects(self, obj1, obj2):
		self.comps.update({'comp1': obj1, 'comp2': obj2})

	def get_Link_conditions(self):
		# condition M'
		if self.comps['comp1'].vFault == self.comps['comp2'].vFault:
			self.link_flag = 0

		# for nonCFds included 2cF, the condition C' needs to be met Note that the compatibility for a-cell is
		# useless. Since the state of a-cell can always be updated correctly, which means that if a tester attempts to
		# use a March sequence to sensitize the two composites, the actual behaviour shown by a-cells is consistent
		# with the tester's expectation. On the other hand, the behaviour of v-cell can be influenced by the fault
		# values, the real state of v-cell and the tester-expected state can be different. Remember that decide
		# whether a 2cF is an LF, the original method is using the March test that designed according to a tester's
		# expectation who DOES NOT KNOW the existence of linked faults. Masking only happens when the final v-cell state of
		# the first-sensitized fault is the same as the initial v-cell state of the second-sensitized fault, because that
		# is the situation of reality in March test.

		elif self.comps['comp1'].CFdsFlag & self.comps['comp2'].CFdsFlag:
			self.link_flag = 1
		else:
			# for nonCFds, any of fault composites cannot be rd or ir fault, or it doesn't satisfy condition D
			if (self.comps['comp1'].rdFlag == 1) or (self.comps['comp2'].rdFlag == 1):
				self.link_flag = 0
			elif (self.comps['comp1'].rdFlag == 2) or (self.comps['comp2'].rdFlag == 2):
				self.link_flag = 0
			elif (self.comps['comp1'].vFault != self.comps['comp2'].vInit) and \
					(self.comps['comp2'].vFault != self.comps['comp1'].vInit):
				self.link_flag = 0
			else:
				self.link_flag = 1
		return


# End of class definition


def parse_fault_pool(fault_pool, fault_model):
	fault_obj_list = get_fault_primitive(fault_pool, fault_model)
	parsed_pool = []
	fault_comps = TwoComposite()
	for obj in fault_obj_list:
		fault_comps.get_FP_text(obj[0])
		fault_comps.get_Comp_objects(obj[1], obj[2])
		fault_comps.get_Link_conditions()
		parsed_pool.append(copy.deepcopy(fault_comps))
	return parsed_pool


def arbit_SF(comp_obj):
	if comp_obj.comps['comp1'] is comp_obj.comps['comp2']:
		return IS_SF
	else:
		return NOT_SF


def arbit_SF_CFds(comp_obj):
	if comp_obj.comps['comp1'].CFdsFlag:
		return IS_SF_CFDS
	else:
		return NOT_SF_CFDS


def arbit_2cF_nonCFds_included(comp_obj):
	if comp_obj.comps['comp1'].CFdsFlag & comp_obj.comps['comp2'].CFdsFlag:
		return NOT_NONCFDS_INCLUDED
	else:
		return IS_NONCFDS_INCLUDED


def arbit_linked_2cF_CFds(comp_obj):
	return comp_obj.link_flag


def classify_based_on_Init(comp_obj):
	ENCODE_0 = 0
	ENCODE_1 = 1
	ENCODE__ = -1

	aInit_dict = {'0': ENCODE_0, '1': ENCODE_1, '-': ENCODE__}
	vInit_dict = {'0': ENCODE_0, '1': ENCODE_1}

	if arbit_SF(comp_obj):
		# for SF, comp1 and comp2 are the same object, use comp1
		if arbit_SF_CFds(comp_obj):
			return vInit_dict.get(comp_obj.comps['comp1'].vInit, CLASSIFY_ERROR)
		else:
			return aInit_dict.get(comp_obj.comps['comp1'].aInit, CLASSIFY_ERROR)
	elif arbit_2cF_nonCFds_included(comp_obj):
		# classify method for filter 2cF by SF
		init_result = tuple()
		for comp_key in comp_obj.comps.keys():
			if comp_obj.comps[comp_key].CFdsFlag:
				init_result += (vInit_dict.get(comp_obj.comps[comp_key].vInit, CLASSIFY_ERROR),)
			else:
				init_result += (aInit_dict.get(comp_obj.comps[comp_key].aInit, CLASSIFY_ERROR),)
		return init_result
	else:
		# The 2 fault composites in a 2cF-CFds need to be classified individually, return a tuple to store the
		# classify results
		return (vInit_dict.get(comp_obj.comps['comp1'].vInit, CLASSIFY_ERROR), vInit_dict.get(
			comp_obj.comps['comp2'].vInit, CLASSIFY_ERROR))


def classify_based_on_SenOpsNum(fault_pool, comp_obj):
	init_encode = classify_based_on_Init(comp_obj)
	if fault_pool['Init_' + str(init_encode)].get('#O_' + str(comp_obj.comps['comp1'].SenOpsNum)) is None:
		fault_pool['Init_' + str(init_encode)]['#O_' + str(comp_obj.comps['comp1'].SenOpsNum)] = set()

	fault_pool['Init_' + str(init_encode)]['#O_' + str(comp_obj.comps['comp1'].SenOpsNum)].add(comp_obj)
	return


def classify_SF(sf_pool, comp_obj):
	classify_based_on_SenOpsNum(sf_pool, comp_obj)
	pass
	return


def classify_2cF_nonCFds_included(_2cF_nonCFds_pool, comp_obj):
	_2cF_nonCFds_pool.add(comp_obj)
	return


def classify_2cF_CFds(_2cF_CFds_pool, comp_obj):
	if arbit_linked_2cF_CFds(comp_obj):
		_2cF_CFds_pool['linked'].add(comp_obj)
	else:
		_2cF_CFds_pool['unlinked'].add(comp_obj)

	return


def classify(unclassified_fault_pool):
	sf_pool = {'Init_0': {}, 'Init_1': {}, 'Init_-1': {}}
	_2cF_nonCFds_pool = set()
	_2cF_CFds_pool = {'linked': set(), 'unlinked': set()}
	for comp_obj in unclassified_fault_pool:
		if arbit_SF(comp_obj):
			classify_SF(sf_pool, comp_obj)
		elif arbit_2cF_nonCFds_included(comp_obj):
			classify_2cF_nonCFds_included(_2cF_nonCFds_pool, comp_obj)
		else:
			classify_2cF_CFds(_2cF_CFds_pool, comp_obj)

	return {'SF': sf_pool, '2cF_nonCFds_included': _2cF_nonCFds_pool, '2cF_CFds': _2cF_CFds_pool}


if __name__ == '__main__':
	classify_result = classify(parse_fault_pool(fault_list_file, fault_model_name))

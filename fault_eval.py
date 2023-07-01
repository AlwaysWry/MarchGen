# a March test evaluation tool #

import sys
import os


# define fault primitive class
class FP:
	"""fault primitive (FP) is a standard format of describing a functional memory
	fault. It contains involved cells, sensitize operations, fault value, etc."""

	def __init__(self, a, v, s, f, cfds, rv):
		"""
		:param a: sensitize initial state of a-cell
		:param v: sensitize initial state of v-cell
		:param s: sensitize sequence
		:param f: fault state
		:param cfds: CFds flag
		:param rv: read-verify operation
		"""
		self.aInit = a
		self.vInit = v
		# TODO: for dynamic faults, Sen is more than 1 ops. Need to add a num_ops property
		self.Sen = s
		self.vFault = f
		self.CFdsFlag = cfds
		self.rvFlag = rv

	def dbg_Get_aInitState(self):
		return self.aInit

	def dbg_Get_vInitState(self):
		return self.vInit

	def dbg_Get_SenOp(self):
		return self.Sen

	def dbg_Get_vFault(self):
		return self.vFault

	def dbg_Get_rvFlag(self):
		return self.rvFlag


def init_cell_state(March):
	if len(March[0].split(',')) == 2:
		# cell_state contains current states of [a1 a2 v] respectively
		if '0' in March[0]:
			cell_state = ['0' for n in range(3)]
		elif '1' in March[0]:
			cell_state = ['1' for n in range(3)]
		else:
			return error
		print(cell_state)
		return cell_state
	else:
		return error


def eval_2comp(FP1, FP2, March):

	cell_state = init_cell_state(March)
	# define cell traverse order
	traverse_up_order = [['a1', 'a2', 'v'], ['a1', 'v', 'a2'], ['a2', 'a1', 'v']]
	traverse_down_order = [['v', 'a2', 'a1'], ['a2', 'v', 'a1'], ['v', 'a1', 'a2']]

	# traverse march elements
	for m, element in enumerate(March[1:]):
		print("element: %s" % element)
		Ops = element.split(',')

		for order_index in range(3):
			if Ops[0] == 'up':
				traverse_cell_order = traverse_up_order[order_index]
			else:
				traverse_cell_order = traverse_down_order[order_index]
			print("current traversing cell order: %s" % traverse_cell_order)

			for i in traverse_cell_order:
				print("Traversing %s, and operation element M%d is (%s)" % (i, m + 1, element))
				for n, op in enumerate(Ops[1:]):  # operations of each element
					print("Operation: %s" % op)
					# a-cell case
					if i != 'v':
						# write operation case
						if 'w' in op:
							# Sensitization operation
							print("FP1.CFdsFlag:%s cell[2]:%s FP1.vInit:%s op:%s FP1.Sen:%s FP1.Init:%s cell[0]:%s" % (
								FP1.CFdsFlag, cell_state[2], FP1.vInit, op, FP1.Sen, FP1.aInit, cell_state[0]))

							if (FP1.CFdsFlag == 1) and (cell_state[2] == FP1.vInit) and (op == FP1.Sen) and (
									FP1.aInit == cell_state[0]) and (i == 'a1'):
								cell_state[2] = FP1.vFault
							elif (FP2.CFdsFlag == 1) and (cell_state[2] == FP2.vInit) and (op == FP2.Sen) and (
									FP2.aInit == cell_state[1]) and (i == 'a2'):
								cell_state[2] = FP2.vFault

							if i == 'a1':
								cell_state[0] = op[1]
							elif i == 'a2':
								cell_state[1] = op[1]

							print("%s: %s, %s, %s" % (i, cell_state[0], cell_state[1], cell_state[2]))

						# read operation case
						elif 'r' in op:
							print("FP1.CFdsFlag: %s cell[2]:%s FP1.vInit:%s op:%s FP1.Sen:%s FP1.Init:%s cell[0]:%s" % (
								FP1.CFdsFlag, cell_state[2], FP1.vInit, op, FP1.Sen, FP1.aInit, cell_state[0]))
							print("%s: %s %s %s" % (i, cell_state[0], cell_state[1], cell_state[2]))
							if (i == 'a1' and cell_state[0] != op[1]) or (i == 'a2' and cell_state[1] != op[1]):
								return Detected
							else:
								if (FP1.CFdsFlag == 1) and (cell_state[2] == FP1.vInit) and (op == FP1.Sen) and (
										FP1.aInit == cell_state[0]) and (i == 'a1'):
									cell_state[2] = FP1.vFault
								elif (FP2.CFdsFlag == 1) and (cell_state[2] == FP2.vInit) and (op == FP2.Sen) and (
										FP2.aInit == cell_state[1]) and (i == 'a2'):
									cell_state[2] = FP2.vFault
								print(
									"CFds read on %s cell: %s %s %s" % (i, cell_state[0], cell_state[1], cell_state[2]))
					# v-cell case
					else:
						# write operation case
						if 'w' in op:
							# Sensitization operation
							if (FP1.CFdsFlag == 0) and (FP1.aInit == cell_state[0]) and (op == FP1.Sen) and (
									cell_state[2] == FP1.vInit):
								cell_state[2] = FP1.vFault
							elif (FP2.CFdsFlag == 0) and (FP2.aInit == cell_state[0]) and (op == FP2.Sen) and (
									cell_state[2] == FP2.vInit):
								cell_state[2] = FP2.vFault
							else:
								cell_state[2] = op[1]

						# read operation case
						elif 'r' in op:
							if (FP1.CFdsFlag == 0) and (FP1.aInit == cell_state[0]) and (op == FP1.Sen) and (
									cell_state[2] == FP1.vInit):
								if FP1.rdFlag == 1:
									return Detected
								else:
									cell_state[2] = FP1.vFault
							elif (FP2.CFdsFlag == 0) and (FP2.aInit == cell_state[0]) and (op == FP2.Sen) and (
									cell_state[2] == FP2.vInit):
								if FP2.rdFlag == 1:
									return Detected
								else:
									cell_state[2] = FP2.vFault
							elif cell_state[2] != op[1]:
								return Detected

						print("%s: %s %s %s" % (i, cell_state[0], cell_state[1], cell_state[2]))

		return Undetected

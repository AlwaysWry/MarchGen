# a March test evaluation tool #

import sys
import os


# define fault primitive class
class FP:

	"""fault primitive (FP) is a standard format of describing a functional memory
	fault. It contains involved cells, sensitize operations, fault value, etc."""
	
	def __init__(self, a, v, s, f, cfds, rv):
		self.aInit = a
		self.vInit = v
		self.SenOp = s
		self.vFault = f
		self.CFdsFlag = cfds
		self.rvFlag = rv

	def dbg_GetAInitState(self):
		return self.aInit

	def dbg_GetVInitState(self):
		return self.vInit

	def dbg_GetSenOp(self):
		return self.SenOp

	def dbg_GetvFault(self):
		return self.vFault

	def dbg_GetrvFlag(self):
		return self.rvFlag


def eval2comp(LF1, LF2, March):

	if len(March[0].split(',')) == 2:
		if '0' in March[0]:
			cell = ['0' for n in range(3)]
		elif '1' in March[0]:
			cell = ['1' for n in range(3)]
		else:
			return error
		print(cell)
	else:
		return error

	# define cell traverse order
	traverse_up_order = [['a1', 'a2', 'v'], ['a1', 'v', 'a2'], ['a2', 'a1', 'v']]
	traverse_down_order = [['v', 'a2', 'a1'], ['a2', 'v', 'a1'], ['v', 'a1', 'a2']]
	order_index = 0

	# traverse March elements
	for m, element in enumerate(March[1:]):    
		print("element: %s" % (element))
		Ops = element.split(',')
		if Ops[0] == 'up':
			traverse_cell_order = traverse_up_order[order_index]

		else:
			traverse_cell_order = traverse_down_order[order_index]

		order_index = order_index + 1

		for i in traverse_cell_order:              
			print("=======Traversing %s cell=======, and operatation element M%d is (%s)" % (i, m + 1, element))
			for n, op in enumerate(Ops[1:]):    # operations of each element
				print("Operation: %s" % (op))
				if i != 'v':
					if ('w' in op):
						print("LF1.CFdsFlag:%s cell[2]:%s LF1.vInit:%s op:%s LF1.Sen:%s LF1.Init:%s cell[0]:%s" %(LF1.CFdsFlag,cell[2],LF1.vInit,op,LF1.Sen,LF1.aInit,cell[0]))
						if (LF1.CFdsFlag == 1) and (cell[2] == LF1.vInit) and (op == LF1.Sen) and (LF1.aInit == cell[0]) and (i==0):
							cell[2] = LF1.vFault
						elif(LF2.CFdsFlag == 1) and (cell[2] == LF2.vInit) and (op == LF2.Sen) and (LF2.aInit == cell[1]) and (i ==1):
							cell[2] = LF2.vFault
						cell[i] = op[1]
						print("%d: %s %s %s" %(i, cell[0], cell[1], cell[2]))

					elif ('r' in op):
						print("LF1.CFdsFlag:%s cell[2]:%s LF1.vInit:%s op:%s LF1.Sen:%s LF1.Init:%s cell[0]:%s" %(LF1.CFdsFlag,cell[2],LF1.vInit,op,LF1.Sen,LF1.aInit,cell[0]))
						print("%d: %s %s %s" %(i, cell[0], cell[1], cell[2]))
						if cell[i] != op[1]:
							return Detected
						else:
							if(LF1.CFdsFlag == 1) and (cell[2] == LF1.vInit) and (op == LF1.Sen) and (LF1.aInit == cell[0]) and (i==0):
								cell[2] = LF1.vFault
							elif(LF2.CFdsFlag == 1) and (cell[2] == LF2.vInit) and (op == LF2.Sen) and (LF2.aInit == cell[1]) and (i==1):
								cell[2] = LF2.vFault
							print("CFds read on %d cell: %s %s %s" %(i, cell[0], cell[1], cell[2]))


				else:            #traversing v-cell
					if ('w' in op):
						#Sensatation operation
						if(LF1.CFdsFlag == 0) and (LF1.aInit == cell[0]) and (op == LF1.Sen) and (cell[2] == LF1.vInit):
							cell[2] = LF1.vFault
						elif(LF2.CFdsFlag == 0) and (LF2.aInit == cell[0]) and (op == LF2.Sen) and (cell[2] == LF2.vInit):
							cell[2] = LF2.vFault
						else: cell[2] = op[1]
							
					elif('r' in op):
						if(LF1.CFdsFlag == 0) and (LF1.aInit == cell[0]) and (op == LF1.Sen) and (cell[2] == LF1.vInit):
							if(LF1.rdFlag == 1): return Detected
							else: cell[2] = LF1.vFault
						elif(LF2.CFdsFlag == 0) and (LF2.aInit == cell[0]) and (op == LF2.Sen) and (cell[2] == LF2.vInit):
							if(LF2.rdFlag == 1): return Detected
							else: cell[2] = LF2.vFault
						elif cell[2] != op[1]:
							return Detected


					print("%d: %s %s %s" %(i, cell[0], cell[1], cell[2]))

	return Undetected

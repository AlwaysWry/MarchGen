import re
import evaluator as ev


def get_March_algorithm(filename):
	file = open(filename, 'r')
	march = []
	for ME in file.readlines():
		if ME.strip().startswith('#'):
			continue
		elif ME.strip():
			march.append(ME.strip())

	file.close()
	print("March test is loaded successfully.\n")
	return march


def get_fault_properties(fault_comps, model):
	fault_props = []
	for fp_index, fp in enumerate(fault_comps, 1):
		empty_dict = {'aCell': '', 'aInit': '',
					  'vCell': '', 'vInit': '',
					  'SenOpsNum': '', 'Sen': '',
					  'vFault': '',
					  'CFdsFlag': '', 'rdFlag': ''}

		vCell = model[0]
		# check if FP is CF
		if ';' in fp:
			if fp_index == 1:
				aCell = model[1]
			elif fp_index == 2:
				aCell = model[2]
			else:
				aCell = ''
				print("Illegal fault %s Found!\n" % fp)
				return fault_props

			fault_prop = fp.partition(';')
			a_prop = fault_prop[0]
			v_prop = fault_prop[2].split('/')
			aInit = a_prop[0]
			vInit = v_prop[0][0]
			vFault = v_prop[1]

			# check if FP is CFds
			if len(a_prop) > 1:
				CFdsFlag = 1
				Sen = a_prop[1:]
			else:
				CFdsFlag = 0
				Sen = v_prop[0][1:]

			# check rd/ir fault
			if (Sen[-2] == 'r') and (Sen[-1] != v_prop[-1]):
				if vFault == Sen[-1]:
					# ir
					rdFlag = 2
				else:
					# rd
					rdFlag = 1
			# check drd fault
			elif (Sen[-2] == 'r') and (Sen[-1] == v_prop[-1]):
				rdFlag = -1
			else:
				rdFlag = 0

		else:
			aCell = vCell
			fault_prop = fp.split('/')
			a_prop = fault_prop[0]
			v_prop = a_prop
			aInit = '-'
			vInit = v_prop[0]
			vFault = fault_prop[1]

			CFdsFlag = 0
			Sen = a_prop[1:]

			# check rd/ir fault
			if (Sen[-2] == 'r') and (Sen[-1] != fault_prop[-1]):
				if vFault == Sen[-1]:
					# ir
					rdFlag = 2
				else:
					# rd
					rdFlag = 1
			# check drd fault
			elif (Sen[-2] == 'r') and (Sen[-1] == fault_prop[-1]):
				rdFlag = -1
			else:
				rdFlag = 0

		SenOpsNum = Sen.count('r') + Sen.count('w')

		props_list = [aCell, aInit, vCell, vInit, SenOpsNum, Sen, vFault, CFdsFlag, rdFlag]
		fault_props_dict = dict(zip(empty_dict.keys(), props_list))
		fault_props.append(fault_props_dict)

	return fault_props


def get_fault_primitive(filename, modelname):
	# all_lines = ['<1;0r0/1/0>*<1;1w1/0/->']
	with open(filename, 'r') as fobj:
		all_lines = fobj.readlines()

	fobj_list = []

	# choose 2cF model
	model_dict = {'2cF_2aa': ['v', 'a1', 'a1'], '2cF_3': ['v', 'a1', 'a2']}

	# The content of lf is like: <S/F/R>*<S/F/R>
	for llf in all_lines:
		lf = llf.strip()
		if '*' in lf:
			# if contains "*", it is a 2-composite fault. extract FP1 and FP2 as ['S/F/R', 'S/F/R']
			fault_comps = re.findall(r"(?<=<).*?(?=>)", lf)
		else:
			fault_comps = [lf[1:-1]]

		fault_props = get_fault_properties(fault_comps, model_dict[modelname])
		if len(fault_props) == 0:
			return fobj_list

		# Put FP1 and FP2 into object FP
		FP1 = ev.FP()
		FP1.__dict__.update(fault_props[0])
		if len(fault_props) > 1:
			FP2 = ev.FP()
			FP2.__dict__.update(fault_props[1])
		else:
			# if it is a simple fault, make FP1 and FP2 the same object
			FP2 = FP1

		fobj_list.append([lf, FP1, FP2])

	return fobj_list

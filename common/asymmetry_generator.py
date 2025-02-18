# a script for generating asymmetry fault sets

def get_sensitization_sequence(op_num):
	def get_subsequence(seq_string):
		return {seq_string + 'w0', seq_string + 'w1', seq_string + 'r' + seq_string[-1]}

	if op_num < 1:
		return {}
	else:
		seq_set_init = {'0w0', '0w1', '1w0', '1w1', '0r0', '1r1'}
		match op_num:
			case 1:
				return seq_set_init
			case _:
				seq_set_temp = seq_set_init.copy()
				seq_set = set()
				while op_num > 1:
					seq_set.clear()
					for seq in seq_set_temp:
						seq_set.update(get_subsequence(seq))
					seq_set_temp.clear()
					seq_set_temp.update(seq_set)
					op_num -= 1

		return seq_set


def get_asymmetry_sequences(seq_set):
	def get_symmetric_sequence(seq):
		seq_temp = seq.replace('0', '-')
		seq_temp = seq_temp.replace('1', '0')
		return seq_temp.replace('-', '1')

	asym_set = set()
	for seq in seq_set:
		seq_sym = get_symmetric_sequence(seq)
		if seq in asym_set:
			continue
		elif seq_sym not in asym_set:
			asym_set.add(seq_sym)

	return asym_set


def get_FR(seq):
	opp_dict = {'0': '1', '1': '0'}
	if seq[-2] == 'r':
		return [(seq[-1], opp_dict[seq[-1]]), (opp_dict[seq[-1]], seq[-1]), (opp_dict[seq[-1]], opp_dict[seq[-1]])]
	else:
		return [(opp_dict[seq[-1]], '-')]


def get_asymmetry_set(seq_set, original_set, wf):
	asym_fp = []
	for seq in seq_set:
		fr_lst = get_FR(seq)
		for fr in fr_lst:

			asym_fp.append('<' + '0;' + seq + '/' + fr[0] + '/' + fr[1] + '>')
			asym_fp.append('<' + '1;' + seq + '/' + fr[0] + '/' + fr[1] + '>')
		asym_fp.append('<' + seq + ';0' + '/1/->')
		asym_fp.append('<' + seq + ';1' + '/0/->')
	for seq in seq_set:
		fr_lst = get_FR(seq)
		for fr in fr_lst:
			asym_fp.append('<' + seq + '/' + fr[0] + '/' + fr[1] + '>')
	for fp in asym_fp:
		wf.write(fp + '\n')


def get_symmetric_fault(fp):
	if ';' in fp:
		fault_prop = fp[1:-1].partition(';')
		a_prop = fault_prop[0]
		v_prop = fault_prop[2].split('/')[0]
		rest_prop = fault_prop[2].split('/')[1:]
	else:
		a_prop = '-'
		v_prop = fp[1:-1].split('/')[0]
		rest_prop = fp[1:-1].split('/')[1:]

	sym_a_prop = a_prop.replace('0', 'x').replace('1', '0').replace('x', '1')
	sym_v_prop = v_prop.replace('0', 'x').replace('1', '0').replace('x', '1')
	sym_rest_prop = [rest.replace('0', 'x').replace('1', '0').replace('x', '1') for rest in rest_prop]
	if a_prop == '-':
		sym_fp = '<' + sym_v_prop + '/' + sym_rest_prop[0] + '/' + sym_rest_prop[1] + '>'
	else:
		sym_fp = '<' + sym_a_prop + ';' + sym_v_prop + '/' + sym_rest_prop[0] + '/' + sym_rest_prop[1] + '>'

	print(sym_fp)
	return sym_fp


if __name__ == '__main__':
	wfilename = '..\\resources\\fault_lists\\asym_1'
	wf = open(wfilename, 'w')
	print(len(get_sensitization_sequence(int(wfilename[-1]))))
	print(s := get_sensitization_sequence(int(wfilename[-1])))
	print(len(get_asymmetry_sequences(s)))
	print(ss := get_asymmetry_sequences(s))
	get_asymmetry_set(ss, s, wf)

	# rfilename = '..\\resources\\fault_lists\\simple_dynamic'
	# rf = open(rfilename, 'r')
	# all_lines = set(fp.strip() for fp in rf.readlines())
	# asym_set = set(filter(lambda f: ';0/' in f or ';1/' in f, all_lines))
	# asym_set.update(set(filter(lambda f: ';' in f and 'w0/' not in f and 'w1/' not in f, all_lines)))
	# asym_set = set(filter(lambda f: '<0;' in f or ';0/' in f, asym_set))
	# # # for line in all_lines:
	# # # 	fp = line.strip()
	# # # 	sym_fp = get_symmetric_fault(fp)
	# # # 	if (sym_fp not in asym_set) and (fp not in asym_set):
	# # # 		asym_set.add(fp)
	# for asym_fp in asym_set:
	# 	wf.write(asym_fp + '\n')
	# rf.close()

	wf.close()


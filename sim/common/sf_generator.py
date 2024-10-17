def get_FR(seq):
	opp_dict = {'0': '1', '1': '0'}
	if seq[-2] == 'r':
		return [(seq[-1], opp_dict[seq[-1]]), (opp_dict[seq[-1]], seq[-1]), (opp_dict[seq[-1]], opp_dict[seq[-1]])]
	else:
		return [(opp_dict[seq[-1]], '-')]


if __name__ == '__main__':
	rfilename = "C:\\Users\\AlwaysWry\\Desktop\\test.txt"
	wfilename = "ss_4dynamic"
	rf = open(rfilename, 'r')
	all_lines = rf.readlines()
	rf.close()
	wf = open(wfilename, 'w')

	opp_dict = {'0': '1', '1': '0'}
	
	seq_lst = []
	for llf in all_lines:
		lf = llf.strip()
		seq_lst.append(lf.split('/')[0][1:])
	for seq in set(seq_lst):
		new_seq = [seq + 'w0', seq + 'w1', seq + 'r' + seq[-1]]
		new_fp = []
		for s in new_seq:
			fr_lst = get_FR(s)
			for fr in fr_lst:
				new_fp.append('<' + s + '/' + fr[0] + '/' + fr[1] + '>')
				new_fp.append('<' + '0;' + s + '/' + fr[0] + '/' + fr[1] + '>')
				new_fp.append('<' + '1;' + s + '/' + fr[0] + '/' + fr[1] + '>')
			new_fp.append('<' + s + ';0' + '/1/->')
			new_fp.append('<' + s + ';1' + '/0/->')
		for fp in new_fp:
			wf.write(fp + '\n')

	wf.close()

	wrf = open(wfilename, 'r')
	all_lines = wrf.readlines()
	print(all_lines)
	if len(all_lines) == len(set(all_lines)):
		print('success')


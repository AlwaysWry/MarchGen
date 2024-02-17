import random

with open('./fault_lists/3_complete', 'r') as r_obj:
	all_lines = r_obj.readlines()
	random_lines = random.sample(all_lines, 18888)

w_obj = open('./fault_lists/random_test', 'w')
for line in random_lines:
	w_obj.write(line.strip())
	w_obj.write('\n')
import os
import random


def random_generator(fault_list, num, part):
	try:
		os.mkdir('./resources/fault_lists/random_tests')
	except FileExistsError:
		pass

	with open(fault_list, 'r') as r_obj:
		all_lines = r_obj.readlines()
		for index in range(10):
			random_lines = random.sample(all_lines, num)
			w_obj = open('./resources/fault_lists/random_tests/' + fault_list_file + '_' + str(part) + '_' + str(index), 'w')
			for line in random_lines:
				w_obj.write(line.strip())
				w_obj.write('\n')
			w_obj.close()


if __name__ == '__main__':
	fault_list_dir = './resources/fault_lists/2cF_2aa/'
	fault_list_file = 'universal_#O_1'
	with open(fault_list_dir + fault_list_file, 'r') as f:
		line_count = sum(1 for line in f)

	random_size = [int(line_count * r / 10) for r in range(1, 11)]

	for list_part, fault_num in enumerate(random_size, 1):
		random_generator(fault_list_dir + fault_list_file, fault_num, list_part)

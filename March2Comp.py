# main script of the March2Comp program
import sys
import time
import datetime

sys.path.append("src/")
sys.path.append("sim/")

from src.output import *
from sim.sim2comp import *


def March2Comp(fault_list, fault_model, fp):
	print("*********************************************************************")
	print("* March2Comp: A March test generator for dynamic 2-composite faults *")
	print("*                             v1.0                                  *")
	print("*********************************************************************")
	start_time = time.time()
	# parse the input fault list
	parsed_faults = parse_fault_pool(fault_list, fault_model)
	if len(parsed_faults) == 0:
		print("Fault list is empty!\n")
		return
	else:
		pass
	# print("Fault list is successfully loaded.\n")
	fp.write("Total fault number: %d\n" % len(parsed_faults))
	# classify
	classified_faults = classify(parsed_faults)
	# print("Classification finished.\n")
	# filter 2cF pool
	filtered_2cFs = filter_redundant_2cF(classified_faults['2cF_nonCFds_included'],
										 classified_faults['2cF_CFds']['unlinked'])
	# filter SF pool
	filtered_SFs = filter_redundant_SF(classified_faults['SF'], filtered_2cFs)
	flat_SFs = flatten_sf_pool(filtered_SFs)
	# create sequence objects
	print("***Generating sensitization units...\n")
	sequence_pool = create_sequence_pool(flat_SFs, filtered_2cFs, classified_faults['2cF_CFds']['linked'])
	# print("Sensitization units are generated.\n")
	# build MEs for linked CFds
	print("***Building march elements...\n")
	me_dict = {'linked_ME': {'main_me': {'01_me': MarchElement(''), '10_me': MarchElement('')}, 'ass_me':
		{'tail_cover_me': MarchElement(''), 'odd_sensitization_me': MarchElement('')}}, 'unlinked_2cF_ME':
				   {'00_me': MarchElement(''), '11_me': MarchElement('')}, 'sf_ME': MarchElement('')}

	if len(sequence_pool['linked']['Init_0']) + len(sequence_pool['linked']['Init_1']) > 0:
		print("Building linked CFds*CFds ME...")
		me_dict['linked_ME'] = linked_CFds_constructor(sequence_pool['linked'])

	if len(sequence_pool['unlinked']['Init_0']) + len(sequence_pool['unlinked']['Init_1']) > 0:
		print("Building unlinked 2cF ME...")
		me_dict['unlinked_2cF_ME'] = unlinked_2cF_constructor(sequence_pool['unlinked'])

	if len(sequence_pool['unlinked']['Init_-1']) > 0:
		print("Building SF ME...")
		me_dict['sf_ME'] = sf_constructor(sequence_pool['unlinked'])

	# print("All elements are built.\n")

	print("\n***Implementing element rearrangement...\n")
	march_element_list = element_assigner(me_dict['linked_ME'], me_dict['unlinked_2cF_ME'], me_dict['sf_ME'])
	print("***March test are successfully generated:\n")

	result = output(march_element_list)
	end_time = time.time()

	fp.write("\nGenerated March test:\n")
	for me in result:
		fp.write("%s\n" % me)
	fp.write("\nMarch element number: %d\n" % len(result))
	fp.write(
		"March test complexity: %dn\n" % sum(map(lambda t: t.count('r') + t.count('w') - t.count('down'), result)))
	fp.write("Elapsed generation time: %lf s.\n" % (end_time - start_time))

	return result


if __name__ == '__main__':
	march_test_file = 'results/march.m2c'
	test_logs_file = 'results/testlog'
	# fault_list_file = 'resources/fault_lists/' + 'complete_with_novel'
	fault_list_file = sys.argv[1]

	with open("results/generation_report.txt", 'w') as report:
		report.write("***********************************************************************************\n")
		report.write("*            March2Comp: A March test generator for 2-composite faults            *\n")
		report.write("*                                    v1.0                                         *\n")
		report.write("*                            Author: Sunrui Zhang                                 *\n")
		report.write("***********************************************************************************\n")
		report.write("\n-------------------------------GENERATION REPORT-----------------------------------\n")
		report.write("\nGeneration starts at %s\n" % datetime.datetime.now())
		report.write("Fault list file: " + fault_list_file + "\n")

		generation_result = March2Comp(fault_list_file, fault_model_name, report)
		march_test = generation_result[0]

		print("\n***Calculating fault coverage...\n")
		print("**Coverage result under 2cF_3 model:")
		fault_model_name = '2cF_3'
		_2cF_3_coverage, _2cF_3_undetected = sim2Comp(march_test_file, test_logs_file, fault_list_file,
													  fault_model_name)
		print("**Coverage result under 2cF_2aa model:")
		fault_model_name = '2cF_2aa'
		_2cF_2aa_coverage, _2cF_2aa_undetected = sim2Comp(march_test_file, test_logs_file, fault_list_file,
														  fault_model_name)

		table_t = '|{:^80s}|\n'
		table = '|{:^4s}{:^18s}{:^4s}|{:^4s}{:^18s}{:^4s}|{:^4s}{:^18s}{:^4s}|\n'
		table_h1 = '|{:^4s}{:>18s}{:^4s}|{:^26s}|{:^26s}|\n'
		table_h2 = '|{:^26s}|{:^4s}{:^18s}{:^4s}|{:^4s}{:^18s}{:^4s}|\n'
		table_h3 = '|{:^4s}{:<18s}{:^4s}|{:^26s}|{:^26s}|\n'
		table_f_1_f_2 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.1f}{:<1s}{:^10s}|{:^10s}{:^.1f}{:<1s}{:^10s}|\n'
		table_2f_1_f_2 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.2f}{:<1s}{:^10s}|{:^10s}{:^.1f}{:<1s}{:^10s}|\n'
		table_2f_2_f_1 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.1f}{:<1s}{:^10s}|{:^10s}{:^.2f}{:<1s}{:^10s}|\n'
		table_2f_1_2f_2 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.2f}{:<1s}{:^10s}|{:^10s}{:^.2f}{:<1s}{:^10s}|\n'
		table_3f_1_f_2 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.3f}{:<1s}{:^10s}|{:^10s}{:^.1f}{:<1s}{:^10s}|\n'
		table_3f_2_f_1 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.1f}{:<1s}{:^10s}|{:^10s}{:^.3f}{:<1s}{:^10s}|\n'
		table_3f_1_2f_2 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.3f}{:<1s}{:^10s}|{:^10s}{:^.2f}{:<1s}{:^10s}|\n'
		table_3f_2_2f_1 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.2f}{:<1s}{:^10s}|{:^10s}{:^.3f}{:<1s}{:^10s}|\n'
		table_3f_1_3f_2 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.3f}{:<1s}{:^10s}|{:^10s}{:^.3f}{:<1s}{:^10s}|\n'

		report.write("\n|--------------------------------------------------------------------------------|\n")
		report.write(table_t.format("FAULT SIMULATION RESULT"))
		report.write("|--------------------------------------------------------------------------------|\n")
		report.write(table_h1.format("", "Fault Models", "", "", ""))
		report.write(table_h2.format("", "", "2cF_3 model", "", "", "2cF_2aa model", ""))
		report.write(table_h3.format("", "Items", "", "", ""))
		report.write("|--------------------------|--------------------------|--------------------------|\n")
		report.write(table.format("", "undetected faults*", "", "", str(len(_2cF_3_undetected)), "", "",
								  str(len(_2cF_2aa_undetected)), ""))
		report.write("|--------------------------|--------------------------|--------------------------|\n")

		if (_2cF_3_coverage >= 100) and (_2cF_2aa_coverage >= 100):
			report.write(
				table_f_1_f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage, "%",
									 ""))
		elif (_2cF_3_coverage < 100) and (_2cF_2aa_coverage < 100):
			if (_2cF_3_coverage < 10) and (_2cF_2aa_coverage < 10):
				report.write(table_3f_1_3f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "",
													_2cF_2aa_coverage, "%", ""))
			elif _2cF_3_coverage < 10:
				report.write(table_3f_1_2f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "",
													_2cF_2aa_coverage, "%", ""))
			elif _2cF_2aa_coverage < 10:
				report.write(table_3f_2_2f_1.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "",
													_2cF_2aa_coverage, "%", ""))
			else:
				report.write(table_2f_1_2f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "",
													_2cF_2aa_coverage, "%", ""))
		elif _2cF_3_coverage < 100:
			if _2cF_3_coverage < 10:
				report.write(
					table_3f_1_f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage,
										  "%", ""))
			else:
				report.write(
					table_2f_1_f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage,
										  "%", ""))
		else:
			if _2cF_2aa_coverage < 10:
				report.write(
					table_3f_2_f_1.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage,
										  "%", ""))
			else:
				report.write(
					table_2f_2_f_1.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage,
										  "%", ""))
		report.write("|--------------------------------------------------------------------------------|\n")
		report.write("*a 2-composite fault are regarded as detected only when it can be detected under all possible cell orders.\n")

		if len(_2cF_3_undetected) > 0:
			report.write("\nUndetected faults in 2cF_3 model:\n")
			for undetected in _2cF_3_undetected:
				report.write(undetected + '\n')

		if len(_2cF_2aa_undetected) > 0:
			report.write("\nUndetected faults in 2cF_2aa model:\n")
			for undetected in _2cF_2aa_undetected:
				report.write(undetected + '\n')

		print("\n***Check elaborate information in \"results/generation_report.txt\".\n")
		report.write("\nSee elaborate test logs in file \"results/testlog\".\n")
		report.write("\n-----------------------------------------------------------------------------------\n")

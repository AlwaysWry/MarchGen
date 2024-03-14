# main script of the March2Comp program
import sys
import time

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
	if not isinstance(parsed_faults, list):
		return
	elif len(parsed_faults) == 0:
		print("Fault list is empty! No March test is generated.\n")
		return
	# print("Fault list is successfully loaded.\n")
	fp.write("Total fault number: %d\n" % len(parsed_faults))

	# classify
	classified_faults = classify(parsed_faults)
	# print("Classification finished.\n")
	classification_time = time.time()

	# filter 2cF pool
	degeneration_result = filter_redundant_2cF(classified_faults['2cF_nonCFds_included'], classified_faults['2cF_CFds']['unlinked'])

	# filter SF pool
	filtered_SFs = filter_redundant_SF(classified_faults['SF'], degeneration_result)
	flat_SFs = flatten_sf_pool(filtered_SFs)
	degeneration_time = time.time()

	# create sequence objects
	print("***Generating sensitization units...\n")
	sequence_pool = create_sequence_pool(flat_SFs, degeneration_result, classified_faults['2cF_CFds']['linked'], classified_faults['2cF_nonCFds_included']['nonCFds_nonCFds'])
	# print("Sensitization units are generated.\n")
	bundle_create_time = time.time()

	# build MEs for linked CFds
	print("***Building march elements...\n")
	me_dict = {'linked_CFds_ME': {'main_me': {'01_me': MarchElement(''), '10_me': MarchElement('')}, 'ass_me':
		{'head_cover_me': MarchElement(''), 'tail_cover_me': MarchElement(''), 'odd_sensitization_me': MarchElement('')}}, 'nonCFds_ME':
				   {'00_me': MarchElement(''), '11_me': MarchElement('')}, 'scf_ME': {MarchElement('')}}
	scf_vertex_pool = set()

	if len(sequence_pool['linked_CFds_seq']['Init_0']) + len(sequence_pool['linked_CFds_seq']['Init_1']) > 0:
		print("Building linked CFds ME...")
		main_construct_result = linked_CFds_constructor(sequence_pool['linked_CFds_seq'],
														classified_faults['2cF_CFds']['linked'])
		me_dict['linked_CFds_ME'] = {'main_me': main_construct_result['main_me'],
									 'ass_me': main_construct_result['ass_me']}
		# if main MEs exist, try to carry out inter ME filter to remove the faults/sequences covered by the main MEs
		unlinked_sequence_pool = inter_ME_filter(me_dict['linked_CFds_ME']['main_me'],
												 main_construct_result['main_me_middle_part'],
												 sequence_pool['degenerated_seq'], sequence_pool['sf_seq'],
												 sequence_pool['undetermined_faults'])
		sequence_pool['degenerated_seq'] = unlinked_sequence_pool['degenerated_seq']
		sequence_pool['sf_seq'] = unlinked_sequence_pool['sf_seq']
		sequence_pool['undetermined_faults'] = unlinked_sequence_pool['undetermined_faults']

	degenerated_size = sum(
		map(lambda k: len(sequence_pool['degenerated_seq'][k]), sequence_pool['degenerated_seq'].keys()))
	sf_size = sum(map(lambda k: len(sequence_pool['sf_seq'][k]), sequence_pool['sf_seq'].keys()))
	if degenerated_size + len(sequence_pool['undetermined_faults']) + sf_size > 0:
		print("Building nonCFds 2cF ME...")
		nonCFds_construct_result = nonCFds_constructor(sequence_pool['degenerated_seq'],
													   sequence_pool['undetermined_faults'], sequence_pool['sf_seq'])
		me_dict['nonCFds_ME'] = nonCFds_construct_result[0]
		scf_vertex_pool = nonCFds_construct_result[1]

	if len(scf_vertex_pool) > 0:
		print("Building SCF ME...")
		me_dict['scf_ME'] = scf_constructor(scf_vertex_pool)

	# print("All elements are built.\n")

	print("\n***Implementing element rearrangement...\n")
	march_element_list = element_assigner(me_dict['linked_CFds_ME'], me_dict['nonCFds_ME'], me_dict['scf_ME'])
	print("***March test are successfully generated:\n")

	result = output(march_element_list)
	end_time = time.time()

	print("\nClassification spends %lf s." % (classification_time - start_time))
	print("Degeneration spends %lf s." % (degeneration_time - start_time))
	print("Total elapsed generation time is %lf s.\n" % (end_time - start_time))

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
	test_logs_file_2cF3 = 'results/testlog_2cF3'
	test_logs_file_2cF_2aa = 'results/testlog_2cF_2aa'
	# fault_list_file = 'resources/fault_lists/' + '2_complete_with_novel'
	fault_list_file = sys.argv[1]
	fault_model_name = default_fault_model_name
	if len(sys.argv) > 2:
		fault_model_name = sys.argv[2]

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
		if generation_result is None:
			sys.exit()

		march_test = generation_result[0]

		print("\n***Calculating fault coverage...\n")
		if fault_model_name == '2cF_2aa':
			print("**Coverage result under 2cF_2aa model:")
			_2cF_2aa_coverage, _2cF_2aa_undetected = sim2Comp(march_test_file, test_logs_file_2cF_2aa, fault_list_file,
														  '2cF_2aa')
		else:
			_2cF_2aa_coverage = -1
			_2cF_2aa_undetected = tuple()

		print("**Coverage result under 2cF_3 model:")
		_2cF_3_coverage, _2cF_3_undetected = sim2Comp(march_test_file, test_logs_file_2cF3, fault_list_file,
													  '2cF_3')


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
		table_3f_2_f_1 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.1f}{:<1s}{:^10s}|{:^10s}{:^.3f}{:<1s}{:^10s}|\n' if fault_model_name == '2cF_2aa' else \
			'|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.1f}{:<1s}{:^10s}|{:^12s}{:d}{:<1s}{:^11s}|\n'
		table_3f_1_2f_2 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.3f}{:<1s}{:^10s}|{:^10s}{:^.2f}{:<1s}{:^10s}|\n'
		table_3f_2_2f_1 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.2f}{:<1s}{:^10s}|{:^10s}{:^.3f}{:<1s}{:^10s}|\n' if fault_model_name == '2cF_2aa' else \
			'|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.2f}{:<1s}{:^10s}|{:^12s}{:^d}{:<1s}{:^11s}|\n'
		table_3f_1_3f_2 = '|{:^4s}{:^18s}{:^4s}|{:^10s}{:^.3f}{:<1s}{:^10s}|{:^10s}{:^.3f}{:<1s}{:^10s}|\n'

		report.write("\n|--------------------------------------------------------------------------------|\n")
		report.write(table_t.format("FAULT SIMULATION RESULT"))
		report.write("|--------------------------------------------------------------------------------|\n")
		report.write(table_h1.format("", "Fault Models", "", "", ""))
		report.write(table_h2.format("", "", "2cF_3 model", "", "", "2cF_2aa model", ""))
		report.write(table_h3.format("", "Items", "", "", ""))
		report.write("|--------------------------|--------------------------|--------------------------|\n")
		report.write(table.format("", "undetected faults*", "", "", str(len(_2cF_3_undetected)), "", "",
								  str(len(_2cF_2aa_undetected)) if fault_model_name == '2cF_2aa' else '-', ""))
		report.write("|--------------------------|--------------------------|--------------------------|\n")

		if (_2cF_3_coverage >= 100) and (_2cF_2aa_coverage >= 100):
			report.write(
				table_f_1_f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage, "%",
									 ""))
		elif (_2cF_3_coverage < 100) and (_2cF_2aa_coverage < 100):
			if (_2cF_3_coverage < 10) and (_2cF_2aa_coverage < 10):
				report.write(table_3f_1_3f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "",
													_2cF_2aa_coverage, "%" if fault_model_name == '2cF_2aa' else '', ""))
			elif _2cF_3_coverage < 10:
				report.write(table_3f_1_2f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "",
													_2cF_2aa_coverage, "%" if fault_model_name == '2cF_2aa' else '', ""))
			elif _2cF_2aa_coverage < 10:
				report.write(table_3f_2_2f_1.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "",
													_2cF_2aa_coverage, "%" if fault_model_name == '2cF_2aa' else '', ""))
			else:
				report.write(table_2f_1_2f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "",
													_2cF_2aa_coverage, "%" if fault_model_name == '2cF_2aa' else '', ""))
		elif _2cF_3_coverage < 100:
			if _2cF_3_coverage < 10:
				report.write(
					table_3f_1_f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage,
										  "%" if fault_model_name == '2cF_2aa' else '', ""))
			else:
				report.write(
					table_2f_1_f_2.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage,
										  "%" if fault_model_name == '2cF_2aa' else '', ""))
		else:
			if _2cF_2aa_coverage < 10:
				report.write(
					table_3f_2_f_1.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage,
										  "%" if fault_model_name == '2cF_2aa' else '', ""))
			else:
				report.write(
					table_2f_2_f_1.format("", "fault coverage", "", "", _2cF_3_coverage, "%", "", "", _2cF_2aa_coverage,
										  "%" if fault_model_name == '2cF_2aa' else '', ""))
		report.write("|--------------------------------------------------------------------------------|\n")
		report.write(
			"*a 2-composite fault are regarded as detected only when it can be detected under all possible cell orders.\n")

		if len(_2cF_3_undetected) > 0:
			report.write("\nUndetected faults in 2cF_3 model:\n")
			for undetected in _2cF_3_undetected:
				report.write(undetected + '\n')
		else:
			report.write("\nAll faults in 2cF_3 model are detected.\n")

		if len(_2cF_2aa_undetected) > 0:
			report.write("\nUndetected faults in 2cF_2aa model:\n")
			for undetected in _2cF_2aa_undetected:
				report.write(undetected + '\n')
		elif fault_model_name == '2cF_2aa':
			report.write("\nAll faults in 2cF_2aa model are detected.\n")

		print("\n***Check elaborate information in \"results/generation_report.txt\".\n")
		if fault_model_name == '2cF_3':
			report.write("\nSee elaborate test logs in file \"results/testlog_2cF3\".\n")
		else:
			report.write("\nSee elaborate test logs in file \"results/testlog_2cF3\" and \"results/testlog_2cF_2aa\".\n")

		report.write("\n-----------------------------------------------------------------------------------\n")

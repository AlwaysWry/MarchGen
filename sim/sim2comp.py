# ---------------------------- #
# a March test evaluation tool #
# ---------------------------- #
import copy
import datetime
import os
import sys
import concurrent.futures as cf
from functools import partial

sys.path.append("../")

from common import fault_parser as ps
import evaluator as ev

# re-define if necessary

# fault_list_file = '../resources/fault_lists/' + '2_complete'
# fault_model_name = '2cF_3'


def atomic_sim(sim_info):
    # logfile.write("\n********\nDetecting fault %s:\n" % lf)
    return ev.eval_2comp(sim_info[0][1], sim_info[0][2], sim_info[1], ev.PROFOUND)
    # sim_flag = result[0]
    # if sim_flag == ev.ERROR:
    #     return ev.ERROR
    # elif sim_flag == ev.UNDETECTED:
    #     return ev.UNDETECTED
    #     # logfile.write("!!! %s is NOT detected !!!\n" % lf)
    #     # undetected_fault.append(lf)
    # else:
    #     return ev.DETECTED
    #     # logfile.write("%s is detected.\n" % lf)


def sim2Comp(test_file, logs_file, fault_list, fault_model):
    logfile = open(logs_file, 'w')
    # logfile = sys.stdout

    march = ps.get_March_algorithm(test_file)
    if not isinstance(march, list):
        return ev.ERROR

    fobj_list = ps.get_fault_primitive(fault_list, fault_model)
    if len(fobj_list) == 0:
        print("\nEmpty or illegal fault list!\n")
        return ev.ERROR

    print("Fault list is successfully loaded.")
    print("Applying March test...")

    undetected_fault = []

    # multi-processor program to accelerate simulation
    sim_info_list = [(f, march) for f in fobj_list]
    with cf.ProcessPoolExecutor() as ppe:
        eval_results = ppe.map(atomic_sim, sim_info_list)

    for fault_obj, eval_result in zip(fobj_list, eval_results):
        logfile.write(f"\n********\nDetecting fault {fault_obj[0]}:\n")
        eval_flag, eval_record = eval_result

        for order_key in eval_record.keys():
            for me_key in eval_record[order_key].keys():
                logfile.write(f"  evaluating element \"{me_key}\" under {order_key}\n")
                match eval_record[order_key][me_key][0]:
                    case 'success':
                        logfile.write(f"    current fault is detected at operation {eval_record[order_key][me_key][1]}.\n")
                    case 'fail':
                        logfile.write("    current fault is NOT detected.\n")

        match eval_flag:
            case ev.ERROR:
                return ev.ERROR
            case ev.UNDETECTED:
                logfile.write(f"!!! {fault_obj[0]} is NOT detected !!!\n")
                undetected_fault.append(fault_obj[0])
            case ev.DETECTED:
                logfile.write(f"{fault_obj[0]} is detected.\n")

    if len(undetected_fault) > 0:
        print("%d linked faults cannot be detected by this March sequence: " % (len(undetected_fault)))
        fault_coverage = ((len(fobj_list) - len(undetected_fault)) / len(fobj_list)) * 100
        print("Fault coverage is %.2f%%" % fault_coverage)
        [print(fault) for fault in undetected_fault]
    else:
        print("Congratulations! All listed faults can be detected by this March sequence. \n")
        fault_coverage = 100
        # [print(element) for element in march]

    if type(logfile) == str:
        logfile.close()

    return fault_coverage, undetected_fault


if __name__ == '__main__':
    os.chdir("../")
    march_test_file = sys.argv[1]
    fault_list_file = sys.argv[2]
    test_logs_file_2cF3 = 'results/testlog_2cF3'
    test_logs_file_2cF_2aa = 'results/testlog_2cF_2aa'

    with open("results/simulation_report.txt", 'w') as report:
        report.write("***********************************************************************************\n")
        report.write("*            March2Comp: A March test generator for 2-composite faults            *\n")
        report.write("*                                    v1.0                                         *\n")
        report.write("*                            Author: Sunrui Zhang                                 *\n")
        report.write("***********************************************************************************\n")
        report.write("\n-------------------------------SIMULATION REPORT-----------------------------------\n")
        report.write("\nGeneration starts at %s\n" % datetime.datetime.now())
        report.write("Fault list file: " + fault_list_file + "\n")
        report.write("March test file: " + march_test_file + "\n")

        print("\n***Calculating fault coverage...\n")
        print("**Coverage result under 2cF_3 model:")
        _2cF_3_coverage, _2cF_3_undetected = sim2Comp(march_test_file, test_logs_file_2cF3, fault_list_file, '2cF_3')
        print("**Coverage result under 2cF_2aa model:")
        _2cF_2aa_coverage, _2cF_2aa_undetected = sim2Comp(march_test_file, test_logs_file_2cF_2aa, fault_list_file, '2cF_2aa')

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

        print("\n***Check elaborate information in \"results/simulation_report.txt\".\n")
        report.write("\nSee elaborate test logs in file \"results/testlog_2cF3\" and \"results/testlog_2cF_2aa\".\n")
        report.write("\n-----------------------------------------------------------------------------------\n")

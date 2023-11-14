# ---------------------------- #
# a March test evaluation tool #
# ---------------------------- #
import os
import sys
sys.path.append("../")

from basic import fault_parser as ps
import evaluator as ev

# re-define if necessary

# fault_list_file = '../resources/fault_lists/' + 'complete'
# fault_model_name = '2cF_3'


def sim2Comp(test_file, logs_file, fault_list_file, fault_model_name):
    logfile = open(logs_file, 'w')
    # logfile = sys.stdout

    march = ps.get_March_algorithm(test_file)

    fobj_list = ps.get_fault_primitive(fault_list_file, fault_model_name)
    if len(fobj_list) == 0:
        print("\nEmpty or illegal fault list!\n")
        return ev.ERROR

    # print("Fault list is successfully loaded.")
    print("Applying March test...")

    undetected_fault = []
    for fobj in fobj_list:
        lf = fobj[0]
        FP1 = fobj[1]
        FP2 = fobj[2]

        logfile.write("\n********\nDetecting fault %s:\n" % lf)
        eval_result = ev.eval_2comp(FP1, FP2, march, ev.PROFOUND, logfile)
        if eval_result == ev.ERROR:
            return ev.ERROR
        elif eval_result == ev.UNDETECTED:
            logfile.write("!!! %s is NOT detected !!!\n" % lf)
            undetected_fault.append(lf)
        else:
            pass
            logfile.write("%s is detected.\n" % lf)

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
    march_test_file = 'resources/march_tests'
    test_logs_file = 'results/testlog'
    sim2Comp(march_test_file, test_logs_file, ps.fault_list_file, ps.fault_model_name)

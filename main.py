# ---------------------------- #
# a March test evaluation tool #
# ---------------------------- #

import sys
import grasp as gr
import evaluator as ev

fault_list_file = 'resources\\fault_lists\\' + 'dyn2comp'
march_test_file = 'resources\\march_tests'
test_logs_file = 'testlog'
fault_model_name = '2cF_3'


def main():
    logfile = open(test_logs_file, 'w')
    # logfile = sys.stdout

    march = gr.get_March_algorithm(march_test_file)

    fobj_list = gr.get_fault_primitive(fault_list_file, fault_model_name)
    if len(fobj_list) == 0:
        print("Empty or illegal fault list!\n")
        return ev.ERROR

    print("Fault list is loaded successfully.\n")
    print("Applying March test...\n")

    Undetected_fault = []
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
            Undetected_fault.append(lf)
        else:
            pass
            logfile.write("%s is detected.\n" % lf)

    if len(Undetected_fault) > 0:
        print("%d linked faults cannot be detected by this March sequence: " % (len(Undetected_fault)))
        [print(fault) for fault in Undetected_fault]
    else:
        print("Congratulations! All listed faults can be detected by this March sequence: \n")
        [print(element) for element in march]

    logfile.close()
    print("\nSee elaborate test logs in file \"testlog\".")


if __name__ == '__main__':
    main()

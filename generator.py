# A script for generating all realistic 2-composite faults
# (includes up to 2-operation dynamic faults)

import main as main
import itertools as it

s_dyn = main.get_fault_primitive("resources/fault_lists/simple_dynamic")
ss_dyn = main.get_fault_primitive("resources/fault_lists/ss_dynamic")
s_stat = main.get_fault_primitive("resources/fault_lists/simple_static")
ss_stat = main.get_fault_primitive("resources/fault_lists/ss_static")

# test = main.get_fault_primitive("test_fault_list")


def get_fp_string(fp_obj):
    return fp_obj[0]


def generate_combinations(fp_obj_list1, fp_obj_list2):
    self1_fc = it.combinations(fp_obj_list1, 2)
    full_combs = list(self1_fc)
    if fp_obj_list1 != fp_obj_list2:
        self2_fc = it.combinations(fp_obj_list2, 2)
        product = it.product(fp_obj_list1, fp_obj_list2)
        full_combs = full_combs + list(self2_fc) + list(product)
    return full_combs


def remove_unrealistic_tuples(full_combs):
    unrealistic_faults = []
    remove_result = full_combs.copy()

    for fobj_tup_index, fobj_tup_item in enumerate(full_combs):
        if (fobj_tup_item[0][1] is not fobj_tup_item[0][2]) or (fobj_tup_item[1][1] is not fobj_tup_item[1][2]):
            print("Illegal fault %s Found!\n" % fobj_tup_item[0])
            break
        else:
            seq_section = 2 * min(fobj_tup_item[0][1].SenOpsNum, fobj_tup_item[1][1].SenOpsNum)
            # no CFds contained. proof see paper
            if (fobj_tup_item[0][1].CFdsFlag == 1) or (fobj_tup_item[1][1].CFdsFlag == 1):
                continue
            elif (fobj_tup_item[0][1].Sen[-2] != 'r') or (fobj_tup_item[1][1].Sen[-2] != 'r'):
                continue
            elif (fobj_tup_item[0][1].Sen[-seq_section:] != fobj_tup_item[1][1].Sen[-seq_section:]) \
                    or (fobj_tup_item[0][1].vInit != fobj_tup_item[1][1].vInit):
                continue
            elif (fobj_tup_item[0][1].vFault == fobj_tup_item[1][1].vFault) and \
                    (fobj_tup_item[0][1].rdFlag == fobj_tup_item[1][1].rdFlag):
                continue

            unrealistic_faults.append((fobj_tup_item[0], fobj_tup_item[1]))
            print("%s" % str((fobj_tup_item[0][0], fobj_tup_item[1][0])))

    print("\nRemoved %d unrealistic faults." % len(unrealistic_faults))

    for fobj_tup in unrealistic_faults:
        remove_result.remove(fobj_tup)

    return remove_result


if __name__ == '__main__':
    fp_combs_obj_list = generate_combinations(s_stat, ss_stat)
    realistic_faults = remove_unrealistic_tuples(fp_combs_obj_list)

    file = open("resources/fault_lists/stat_stat_fault_list", 'w')
    for tup in realistic_faults:
        lf = get_fp_string(tup[0]) + '*' + get_fp_string(tup[1]) + '\n'
        file.write(lf)
    print("fault list is generated.")

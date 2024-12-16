# A script for generating all realistic 2-composite faults
# (includes up to 2-operation dynamic faults)

import fault_parser as ps
import itertools as it

model_name_2cF3 = '2cF_3'
model_name_2cF2aa = '2cF_2aa'

s_dyn = ps.get_fault_primitive("../resources/fault_lists/simple_dynamic", model_name_2cF3)
ss_dyn = ps.get_fault_primitive("../resources/fault_lists/ss_dynamic", model_name_2cF3)
s_stat = ps.get_fault_primitive("../resources/fault_lists/simple_static", model_name_2cF3)
ss_stat = ps.get_fault_primitive("../resources/fault_lists/ss_static", model_name_2cF3)
s_3dyn = ps.get_fault_primitive("../resources/fault_lists/simple_3dynamic", model_name_2cF3)
s_4dyn = ps.get_fault_primitive("../resources/fault_lists/simple_4dynamic", model_name_2cF3)
ss_3dyn = ps.get_fault_primitive("../resources/fault_lists/ss_3dynamic", model_name_2cF3)
s_max2 = ps.get_fault_primitive("../resources/fault_lists/single_fault_max2", model_name_2cF3)
s_max3 = ps.get_fault_primitive("../resources/fault_lists/single_fault_max3", model_name_2cF3)
asym_stat = ps.get_fault_primitive("../resources/fault_lists/asym_1", model_name_2cF3)
asym_dyn = ps.get_fault_primitive("../resources/fault_lists/asym_2", model_name_2cF3)
asym_3dyn = ps.get_fault_primitive("../resources/fault_lists/asym_3", model_name_2cF3)
asym_4dyn = ps.get_fault_primitive("../resources/fault_lists/asym_4", model_name_2cF3)
asym_max2 = ps.get_fault_primitive("../resources/fault_lists/asym_max2", model_name_2cF3)
asym_max3 = ps.get_fault_primitive("../resources/fault_lists/asym_max3", model_name_2cF3)

# test = ps.get_fault_primitive("test_fault_list")


def get_fp_string(fp_obj):
    return fp_obj[0]


def generate_combinations(fp_obj_list1, fp_obj_list2):
    self1_fc = it.combinations(fp_obj_list1, 2)
    full_combs = list(self1_fc)
    if fp_obj_list1 != fp_obj_list2:
        # self2_fc = it.combinations(fp_obj_list2, 2)
        self2_fc = []
        product = it.product(fp_obj_list1, fp_obj_list2)
        full_combs = full_combs + list(self2_fc) + list(product)
    return full_combs


def remove_unrealistic_tuples(full_combs, modelname):
    unrealistic_faults = []
    removed_result = full_combs.copy()

    for fobj_tup_index, fobj_tup_item in enumerate(full_combs):
        if (fobj_tup_item[0][1] is not fobj_tup_item[0][2]) or (fobj_tup_item[1][1] is not fobj_tup_item[1][2]):
            print("Illegal fault %s Found!\n" % fobj_tup_item[0])
            break
        else:
            seq_section = 2 * min(fobj_tup_item[0][1].SenOpsNum, fobj_tup_item[1][1].SenOpsNum)
            # Cannot sensitize simultaneously if CFds contained. proof see paper [TCAD'12]
            if (fobj_tup_item[0][1].CFdsFlag == 1) or (fobj_tup_item[1][1].CFdsFlag == 1):
                continue
            elif (fobj_tup_item[0][1].Sen[-2] != 'r') or (fobj_tup_item[1][1].Sen[-2] != 'r'):
                continue
            elif (fobj_tup_item[0][1].Sen[-seq_section:] != fobj_tup_item[1][1].Sen[-seq_section:]) \
                    or (fobj_tup_item[0][1].vInit != fobj_tup_item[1][1].vInit):
                continue
            # if the faulty value and the read-out value of the 2 FPs are identical, the simultaneous sensitization will not
            # cause uncertain state on v-cell, which means that the 2cF is realistic
            elif (fobj_tup_item[0][1].vFault == fobj_tup_item[1][1].vFault) and \
                    (fobj_tup_item[0][1].rdFlag == fobj_tup_item[1][1].rdFlag):
                continue
            # for 2cF2aa fault model, the <x1;y...rz/F/R>*<x2;y...rz/F/R> (x1 != x2) is realistic. See [TCAD'12]
            elif ((modelname == '2cF_2aa') and (fobj_tup_item[0][1].aInit != '-' != fobj_tup_item[1][1].aInit)
                  and (fobj_tup_item[0][1].aInit != fobj_tup_item[1][1].aInit)):
                continue

            unrealistic_faults.append((fobj_tup_item[0], fobj_tup_item[1]))
            # print("%s" % str((fobj_tup_item[0][0], fobj_tup_item[1][0])))

    print("\nRemoved %d unrealistic faults." % len(unrealistic_faults))

    for fobj_tup in unrealistic_faults:
        removed_result.remove(fobj_tup)

    return removed_result


if __name__ == '__main__':
    fp_combs_obj_list = generate_combinations(asym_stat, asym_stat)
    realistic_faults_2cF3 = remove_unrealistic_tuples(fp_combs_obj_list, model_name_2cF3)
    realistic_faults_2cF2aa = remove_unrealistic_tuples(fp_combs_obj_list, model_name_2cF2aa)

    with open("../resources/fault_lists/2cF_3/asym_1_2cF3", 'w') as file:
        for tup in realistic_faults_2cF3:
            lf = get_fp_string(tup[0]) + '*' + get_fp_string(tup[1]) + '\n'
            file.write(lf)

    with open("../resources/fault_lists/2cF_2aa/asym_1_2cF2aa", 'w') as file:
        for tup in realistic_faults_2cF2aa:
            lf = get_fp_string(tup[0]) + '*' + get_fp_string(tup[1]) + '\n'
            file.write(lf)

    print("fault list is generated.")

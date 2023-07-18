import evaluator as ev
import re

# ------------
#
# ------------

# march = ["up,w0",
#         "up,r0,r0,w0,r0,w1,w0,r0,w1,w1,r1",
#         "up,r1,r1,w1,r1,w0,w1,r1,w0,w0,r0",
#         "up,r0"
#        ]

# march = ["up,w0",
#         "up,r0,w0,r0,w1,w1",
#         "up,r1,w1,r1,w0,w0,r0",
#         "up,r0,w1,r1,r1,w0",
#         "up,r1,w0,r0,r0,w1,r1"
#        ]

# march = ["up,w0",
#         "up,r0,w0,r0,w1,w1,r1",
#         "up,r1,w0,r0,r0,w1",
#         "up,r1,w1,r1,w0,w0,r0",
#         "up,r0,w1,r1,r1,w0"
#        ]

# march = ["up,w0",
#         "up,r0,w1,r1,r1,w0,r0",
#         "up,r0,w0,r0,w1,w1",
#         "up,r1,w0,r0,r0,w1,r1",
#         "up,r1,w1,r1,w0,w0"
#        ]

# march = ["up,w0",
#         "up,r0,w0,r0,w1,w1",
#         "up,r1,w1,r1,w0,w0",
#         "up,r0,w1,w0,w1,r1",
#         "up,r1,w0,w1,w0,r0"
#        ]

# march = ["up,w1",
#         "down,r1,w0,r0,w0,r0",
#         "down,r0,w1,r1,w1,r1",
#         "up,r1,w0,r0,w0,r0",
#         "up,r0,w1,r1,w1,r1",
#         "up,r1"]

# march = ["up,w0",
#         "up,r0,w0,r0,w1,w1,r1,w0,r0",
#         "up,r0,w1",
#         "up,r1,w1,r1,w0,w0,r0,w1,r1",#non-CFds对结构还是有影响，去掉之后不能检测a-cell为1的CFwd
#         "up,r1"]

# March BDN
# march = ["up,w0",
#          "down,r0,w1,r1,w1,r1",
#          "down,r1,w0,r0,w0,r0",
#          "up,r0,w1,r1,w1,r1",
#          "up,r1,w0,r0,w0,r0",
#          "up,r0"]

# march = ["up,w0",
#          "up,r0,w0,w1,w1,r1,w0,r0",
#          "up,r0,w1",
#          "up,r1,w1,w0,w0,r0,w1,r1",
#          "up,r1"]

# march = ["up,w0",
#          "up,r0,w0,r0,w1,w0,r0,r0,w1,w1,r1",
#          "up,r1,w1,r1,w0,w1,r1,r1,w0,w0,r0",
#          "up,r0"]

march = ["any,w0",
         "any,r0,w0,r0,r0,r0,w1,r1,w1,r1,r1,w0,r0,r0,w0,w1,w0,w1,r1,w1,w0,w1,w0,r0,w1,w1,w1,w1,r1,r1,r1,w0,w0,w0,w0,r0",
         "any,r0",
         "any,w1",
         "any,r1,w1,r1,r1,r1,w0,r0,w0,r0,r0,w1,r1,r1,w1,w0,w1,w0,r0,w0,w1,w0,w1,r1,w0,w0,w0,w0,r0,r0,r0,w1,w1,w1,w1,r1",
         "any,r1"]

fault_list_file = 'linked_fault'


def get_fault_properties(fault_comps):
    fault_props = []
    for fp_index, fp in enumerate(fault_comps, 1):
        empty_dict = {'aCell': '', 'aInit': '',
                      'vCell': '', 'vInit': '',
                      'SenOpsNum': '', 'Sen': '',
                      'vFault': '',
                      'CFdsFlag': '', 'rdFlag': ''}

        vCell = 'v'
        # check if FP is CF
        if ';' in fp:
            if fp_index == 1:
                aCell = 'a1'
            elif fp_index == 2:
                aCell = 'a2'
            else:
                aCell = ''
                print("Illegal fault %s Found!\n" % fp)
                return fault_props

            fault_prop = fp.partition(';')
            a_prop = fault_prop[0]
            v_prop = fault_prop[2].split('/')
            aInit = a_prop[0]
            vInit = v_prop[0][0]
            vFault = v_prop[1]

            # check if FP is CFds
            if len(a_prop) > 1:
                CFdsFlag = 1
                Sen = a_prop[1:]
            else:
                CFdsFlag = 0
                Sen = v_prop[0][1:]

            # check rd/ir fault
            if (Sen[-2] == 'r') and (Sen[-1] != v_prop[-1]):
                if vFault == Sen[-1]:
                    # ir
                    rdFlag = 2
                else:
                    # rd
                    rdFlag = 1
            # check drd fault
            elif (Sen[-2] == 'r') and (Sen[-1] == v_prop[-1]):
                rdFlag = -1
            else:
                rdFlag = 0

        else:
            aCell = vCell
            fault_prop = fp.split('/')
            a_prop = fault_prop[0]
            v_prop = a_prop
            aInit = '-'
            vInit = v_prop[0]
            vFault = fault_prop[1]

            CFdsFlag = 0
            Sen = a_prop[1:]

            # check rd/ir fault
            if (Sen[-2] == 'r') and (Sen[-1] != fault_prop[-1]):
                if vFault == Sen[-1]:
                    # ir
                    rdFlag = 2
                else:
                    # rd
                    rdFlag = 1
            # check drd fault
            elif (Sen[-2] == 'r') and (Sen[-1] == fault_prop[-1]):
                rdFlag = -1
            else:
                rdFlag = 0

        SenOpsNum = Sen.count('r') + Sen.count('w')

        props_list = [aCell, aInit, vCell, vInit, SenOpsNum, Sen, vFault, CFdsFlag, rdFlag]
        fault_props_dict = dict(zip(empty_dict.keys(), props_list))
        fault_props.append(fault_props_dict)

    return fault_props


def get_fault_primitive(filename):
    # all_lines = ['<1;0r0/1/0>*<1;1w1/0/->']
    with open(filename, 'r') as fobj:
        all_lines = fobj.readlines()

    fobj_list = []

    # The content of lf is like: <S/F/R>*<S/F/R>
    for llf in all_lines:
        lf = llf.strip()
        if '*' in lf:
            # if contains "*", it is a 2-composite fault. extract FP1 and FP2 as ['S/F/R', 'S/F/R']
            fault_comps = re.findall(r"(?<=<).*?(?=>)", lf)
        else:
            fault_comps = [lf[1:-1]]

        fault_props = get_fault_properties(fault_comps)
        if len(fault_props) == 0:
            return fobj_list

        # Put FP1 and FP2 into object FP
        FP1 = ev.FP()
        FP1.__dict__.update(fault_props[0])
        if len(fault_props) > 1:
            FP2 = ev.FP()
            FP2.__dict__.update(fault_props[1])
        else:
            # if it is a simple fault, make FP1 and FP2 the same object
            FP2 = FP1

        fobj_list.append([lf, FP1, FP2])

    return fobj_list


def main(filename):
    fobj_list = get_fault_primitive(filename)
    if len(fobj_list) == 0:
        return ev.ERROR

    Undetected_fault = []

    for fobj in fobj_list:
        lf = fobj[0]
        FP1 = fobj[1]
        FP2 = fobj[2]

        print("\nDetecting fault %s" % lf)
        if ev.eval_2comp(FP1, FP2, march, ev.PROFOUND) == ev.UNDETECTED:
            print("!!! %s is NOT detected !!!" % lf)
            Undetected_fault.append(lf)
        else:
            pass
            print("%s is detected." % lf)

    if len(Undetected_fault) > 0:
        print("%d linked faults cannot be detected by this March sequence: " % (len(Undetected_fault)))
        [print(fault) for fault in Undetected_fault]
    else:
        print("\nCongratulations! All linked faults can be detected by this March sequence: ")
        [print(element) for element in march]


if __name__ == '__main__':
    main(fault_list_file)

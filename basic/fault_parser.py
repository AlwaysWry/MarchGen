# basic module of parsing the input fault list
import re

# unparsed fault list definitions
fault_list_file = '../resources/fault_lists/' + 'complete'
fault_model_name = '2cF_3'

NEST = True
NOT_NEST = False


# define fault primitive class
class SimpleFault:
    """simple fault primitive (SimpleFault) is a standard format of describing a functional memory
    fault. It contains involved cells, sensitize operations, fault value, etc."""

    props = {'aCell': '', 'aInit': '',
             'vCell': '', 'vInit': '',
             'SenOpsNum': '', 'Sen': '',
             'vFault': '',
             'CFdsFlag': '', 'rdFlag': '',
             'nestSenFlag': ''}

    fp_text = ''

    def __init__(self):
        """
        aCell: aggressor cell index
        aInit: sensitize initial state of a-cell
        vCell: victim cell index
        vInit: sensitize initial state of v-cell
        SenOpsNum: number of sensitize operations
        Sen: sensitize sequence
        vFault: fault state
        CFdsFlag: CFds flag
        rdFlag: read-verify operation
        nestSenFlag: nest sensitization available flag
        """

        self.__dict__.update(self.props)

    def dbg_Get_aInit(self):
        return self.__dict__['aInit']

    def dbg_Get_vInit(self):
        return self.__dict__['vInit']

    def dbg_Get_SenOpsNum(self):
        return self.__dict__['SenOpsNum']

    def dbg_Get_Sen(self):
        return self.__dict__['Sen']

    def dbg_Get_vFault(self):
        return self.__dict__['vFault']

    def dbg_Get_rdFlag(self):
        return self.__dict__['rdFlag']

    def dbg_Get_nestSenFlag(self):
        return self.__dict__['nestSenFlag']


def get_March_algorithm(filename):
    march = []

    with open(filename, 'r') as file:
        for ME in file:
            if ME.strip().startswith('#'):
                continue
            elif ME.strip():
                march.append(ME.strip())

    print("March test is successfully loaded.\n")
    return march


def arbit_nest_sensitization(Sen, CFdsFlag):
    operation_num = int(len(Sen) / 2)
    if not ((CFdsFlag == 0) and (Sen.startswith('w'))):
        return NOT_NEST
    elif operation_num % 2 == 0:
        # the number of Sen operations is even
        if Sen[:operation_num] == Sen[operation_num:]:
            return NOT_NEST
        else:
            return NEST
    elif operation_num > 1:
        if Sen[:operation_num + 2] == Sen[operation_num - 1:]:
            return NOT_NEST
        else:
            return NEST
    else:
        return NEST


def get_fault_properties(fault_comps, model):
    fault_props = []
    for fp_index, fp in enumerate(fault_comps, 1):
        empty_dict = {'aCell': '', 'aInit': '',
                      'vCell': '', 'vInit': '',
                      'SenOpsNum': '', 'Sen': '',
                      'vFault': '',
                      'CFdsFlag': '', 'rdFlag': '',
                      'nestSenFlag': ''}

        vCell = model[0]
        # check if SimpleFault is CF
        if ';' in fp:
            if fp_index == 1:
                aCell = model[1]
            elif fp_index == 2:
                aCell = model[2]
            else:
                aCell = ''
                print("Illegal fault %s found!\n" % fp)
                return

            fault_prop = fp.partition(';')
            a_prop = fault_prop[0]
            v_prop = fault_prop[2].split('/')

            try:
                aInit = a_prop[0]
                vInit = v_prop[0][0]
                vFault = v_prop[1]
            except IndexError:
                print("Illegal fault primitive. Check the fault list.\n")
                return

            # check if SimpleFault is CFds
            if len(a_prop) > 1:
                CFdsFlag = 1
                Sen = a_prop[1:]
            else:
                CFdsFlag = 0
                Sen = v_prop[0][1:]

            # check rd/ir feature
            if (Sen[-2] == 'r') and (Sen[-1] != v_prop[-1]):
                if vFault == Sen[-1]:
                    # ir
                    rdFlag = 2
                else:
                    # rd
                    rdFlag = 1
            # check drd feature
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

            # check rd/ir feature
            if (Sen[-2] == 'r') and (Sen[-1] != fault_prop[-1]):
                if vFault == Sen[-1]:
                    # ir
                    rdFlag = 2
                else:
                    # rd
                    rdFlag = 1
            # check drd feature
            elif (Sen[-2] == 'r') and (Sen[-1] == fault_prop[-1]):
                rdFlag = -1
            else:
                rdFlag = 0

        SenOpsNum = Sen.count('r') + Sen.count('w')

        # for non-CFds, check the nest sensitization conditions
        if arbit_nest_sensitization(Sen, CFdsFlag):
            if vInit != Sen[-1]:
                nestSenFlag = 'donor'
            else:
                nestSenFlag = 'receiver'
        else:
            nestSenFlag = 'invalid'

        props_list = [aCell, aInit, vCell, vInit, SenOpsNum, Sen, vFault, CFdsFlag, rdFlag, nestSenFlag]
        fault_props_dict = dict(zip(empty_dict.keys(), props_list))
        fault_props.append(fault_props_dict)

    return fault_props


def get_fault_primitive(filename, modelname):
    # all_lines = ['<1;0r0/1/0>*<1;1w1/0/->']
    try:
        fobj = open(filename, 'r')
    except OSError:
        print("Open fault list file failed. Make sure the file path and name is correct.\n")
        return

    all_lines = fobj.readlines()
    fobj.close()

    fobj_list = []

    # choose 2cF model
    model_dict = {'2cF_2aa': ['v', 'a1', 'a1'], '2cF_3': ['v', 'a1', 'a2']}

    # The content of lf is like: <S/F/R>*<S/F/R>
    for llf in all_lines:
        lf = llf.strip()
        if '*' in lf:
            # if contains "*", it is a 2-composite fault. extract FP1 and FP2 as ['S/F/R', 'S/F/R']
            fault_comps = re.findall(r"(?<=<).*?(?=>)", lf)
        else:
            fault_comps = [lf[1:-1]]

        fault_props = get_fault_properties(fault_comps, model_dict[modelname])
        if len(fault_props) == 0:
            return fobj_list

        # Put FP1 and FP2 into object SimpleFault
        FP1 = SimpleFault()
        FP1.__dict__.update(fault_props[0])
        FP1.fp_text = '<' + fault_comps[0] + '>'
        if len(fault_props) > 1:
            FP2 = SimpleFault()
            FP2.__dict__.update(fault_props[1])
            FP2.fp_text = '<' + fault_comps[1] + '>'
        else:
            # if it is a simple fault, make FP1 and FP2 the same object
            FP2 = FP1

        fobj_list.append([lf, FP1, FP2])

    return fobj_list

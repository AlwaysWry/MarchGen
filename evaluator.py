# a March test evaluation tool #

# define fault primitive class
class FP:
    """fault primitive (FP) is a standard format of describing a functional memory
    fault. It contains involved cells, sensitize operations, fault value, etc."""

    props = {'aCell': '', 'aInit': '',
             'vCell': '', 'vInit': '',
             'SenOpsNum': '', 'Sen': '',
             'vFault': '',
             'CFdsFlag': '', 'rdFlag': ''}

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

    def dbg_Get_rvFlag(self):
        return self.__dict__['rdFlag']


# return state definitions
INIT_ERROR = 0
ERROR = 0
UPDATED = 1
UPDATE_ERROR = 0
APPLIED = 1
DETECTED = 10
UNDETECTED = 100
PROFOUND = 1
FAST = 0

def init_cell_state(March):
    if len(March[0].split(',')) == 2:
        # cell_state is a dictionary, which contains current states of a1, a2, v respectively
        if '0' in March[0]:
            cell_state = {'a1': '0', 'a2': '0', 'v': '0'}
        elif '1' in March[0]:
            cell_state = {'a1': '1', 'a2': '1', 'v': '1'}
        else:
            return INIT_ERROR
        # print(cell_state)
        return cell_state
    else:
        return INIT_ERROR


def init_cell_order(mode):
    cell_order = []
    if mode == FAST:
        traverse_up_order = [['a1', 'a2', 'v']]
        traverse_down_order = [['v', 'a2', 'a1']]
    else:
        traverse_up_order = [['a1', 'a2', 'v'], ['a1', 'v', 'a2'], ['a2', 'a1', 'v']]
        traverse_down_order = [['v', 'a2', 'a1'], ['a2', 'v', 'a1'], ['v', 'a1', 'a2']]
    for tup in list(zip(traverse_up_order + traverse_down_order, traverse_down_order + traverse_up_order)):
        cell_order.append({'up': tup[0], 'down': tup[1]})
    return cell_order


def generate_operation_sequence(SenOpNum, op_num, ops):
    op_seq = ''
    for i in range(SenOpNum - 1, -1, -1):
        op_seq = op_seq + ops[op_num - i - 1]
    return op_seq


def update_fault_state(fp, cell_state, cell_history, visiting_cell, op_seq, op):
    if visiting_cell != 'v':
        if (fp.CFdsFlag == 1) and (cell_state[fp.vCell] == fp.vInit) and (op_seq == fp.Sen) and (
                cell_history[-fp.SenOpsNum][fp.aCell] == fp.aInit) and (visiting_cell == fp.aCell):
            cell_state[fp.vCell] = fp.vFault
            return UPDATED
    else:
        if 'w' in op:
            if (fp.CFdsFlag == 0) and (cell_state[fp.aCell] == fp.aInit) and (op_seq == fp.Sen) and (
                    cell_history[-fp.SenOpsNum][fp.vCell] == fp.vInit):
                cell_state[fp.vCell] = fp.vFault
                return UPDATED
        elif 'r' in op:
            if (fp.CFdsFlag == 0) and (fp.aInit == cell_state[fp.aCell]) and (op_seq == fp.Sen) and (
                    cell_history[-fp.SenOpsNum][fp.vCell] == fp.vInit):
                # check rdflag, if it is read-destruction fault, it is detected immediately
                if fp.rdFlag == 1:
                    return DETECTED
                else:
                    cell_state[fp.vCell] = fp.vFault
                    return UPDATED
            elif cell_state[fp.vCell] != op[1]:
                return DETECTED

    return UPDATE_ERROR


def apply_March_sequence(cell_state, traverse_cell_order, OPS, FP1, FP2):
    cell_history = []

    def set_cell_history():
        nonlocal cell_history
        if len(cell_history) == 0:
            for state in range(max(FP1.SenOpsNum, FP2.SenOpsNum)):
                cell_history.append(cell_state.copy())
        else:
            del cell_history[0]
            cell_history.append(cell_state.copy())
        return

    set_cell_history()

    for visiting_cell in traverse_cell_order:
        # operations of each element, index starts from 1
        ops = OPS[1:]
        update_flag = ['', '']
        for op_num, op in enumerate(ops, 1):
            # print("Operation: %s" % op)
            # a-cell case
            if visiting_cell != 'v':
                # write operation case
                if 'w' in op:
                    # Sensitization operation
                    if op_num >= FP1.SenOpsNum:
                        op_seq = generate_operation_sequence(FP1.SenOpsNum, op_num, ops)
                        update_fault_state(FP1, cell_state, cell_history, visiting_cell, op_seq, op)
                    #    if update_state == DETECTED:
                    #        return DETECTED

                    if FP2 is not FP1:
                        if op_num >= FP2.SenOpsNum:
                            op_seq = generate_operation_sequence(FP2.SenOpsNum, op_num, ops)
                            update_fault_state(FP2, cell_state, cell_history, visiting_cell, op_seq, op)
                        #    if update_state == DETECTED:
                        #        return DETECTED

                    # for write operation on a-cells, a-cell state will be changed
                    # no matter whether the fault is sensitized on v-cell
                    cell_state[visiting_cell] = op[1]

                    set_cell_history()

                # read operation case
                elif 'r' in op:
                    if cell_state[visiting_cell] != op[1]:
                        return DETECTED
                    else:
                        if op_num >= FP1.SenOpsNum:
                            op_seq = generate_operation_sequence(FP1.SenOpsNum, op_num, ops)
                            update_flag[0] = update_fault_state(FP1, cell_state, cell_history, visiting_cell, op_seq, op)
                            if update_flag[0] == DETECTED:
                                return DETECTED

                        if FP2 is not FP1:
                            if op_num >= FP2.SenOpsNum:
                                op_seq = generate_operation_sequence(FP2.SenOpsNum, op_num, ops)
                                update_flag[1] = update_fault_state(FP2, cell_state, cell_history, visiting_cell, op_seq, op)
                                if update_flag[1] == DETECTED:
                                    return DETECTED

                        set_cell_history()

            # v-cell case
            else:
                # write operation case
                if 'w' in op:
                    # Sensitization operation
                    if op_num >= FP1.SenOpsNum:
                        op_seq = generate_operation_sequence(FP1.SenOpsNum, op_num, ops)
                        update_flag[0] = update_fault_state(FP1, cell_state, cell_history, visiting_cell, op_seq, op)
                    #    if update_state == DETECTED:
                    #        return DETECTED
                    else:
                        update_flag[0] = UPDATE_ERROR

                    if FP2 is not FP1:
                        if op_num >= FP2.SenOpsNum:
                            op_seq = generate_operation_sequence(FP2.SenOpsNum, op_num, ops)
                            update_flag[1] = update_fault_state(FP2, cell_state, cell_history, visiting_cell, op_seq, op)
                        #    if update_state == DETECTED:
                        #        return DETECTED
                    else:
                        update_flag[1] = UPDATE_ERROR

                    # for write operation on v-cell, if the fault is sensitized,
                    # the final state of v-cell is the fault state, and the write operation will be suppressed
                    if UPDATED not in update_flag:
                        cell_state[visiting_cell] = op[1]

                    set_cell_history()

                # read operation case
                elif 'r' in op:
                    # Sensitization operation
                    if op_num >= FP1.SenOpsNum:
                        op_seq = generate_operation_sequence(FP1.SenOpsNum, op_num, ops)
                        if update_fault_state(FP1, cell_state, cell_history, visiting_cell, op_seq, op) == DETECTED:
                            return DETECTED
                    elif cell_state[visiting_cell] != op[1]:
                        return DETECTED

                    if FP2 is not FP1:
                        if op_num >= FP2.SenOpsNum:
                            op_seq = generate_operation_sequence(FP2.SenOpsNum, op_num, ops)
                            if update_fault_state(FP2, cell_state, cell_history, visiting_cell, op_seq, op) == DETECTED:
                                return DETECTED
                        elif cell_state[visiting_cell] != op[1]:
                            return DETECTED

                    set_cell_history()

    return APPLIED


def eval_2comp(FP1, FP2, March, mode):
    eval_flag = []
    cell_order = init_cell_order(mode)

    for order_index in range(len(cell_order)):
        cell_state = init_cell_state(March)
        # traverse march elements
        for m, element in enumerate(March[1:]):
            # set cell traverse order
            if 'down' in element:
                traverse_cell_order = cell_order[order_index]['down']
            else:
                traverse_cell_order = cell_order[order_index]['up']

            print("  evaluating element \"%s\" under %s" % (element, str(traverse_cell_order)))
            ops = element.split(',')
            if apply_March_sequence(cell_state, traverse_cell_order, ops, FP1, FP2) == DETECTED:
                eval_flag.append('success')
                print("    current fault is detected.")
                break
            else:
                print("    current fault is NOT detected.")

    if eval_flag.count('success') == len(cell_order):
        return DETECTED
    else:
        return UNDETECTED

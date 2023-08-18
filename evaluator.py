# a March test evaluation tool #
import copy


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

    def dbg_Get_rdFlag(self):
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


def update_fault_state(fp, cell_state, cell_state_temp, cell_history, visiting_cell, op_seq, op):
    if visiting_cell != 'v':
        # for dynamic faults, the init state of a-cell/v-cell should be kept in cell_history, or it will be
        # covered by new state. The init state also cannot be destroyed by fault state of other fault in the
        # middle of sensitization, like <0w1/0/->*<0w1w0/1/->, the FP2 can be sensitized by "O0, w1, w0",
        # because processing the fault means that it can be sensitized under correct conditions, i.e. "dynamic
        # faulty behavior can take place in the absence of static faults."
        # [Dynamic Faults in Random-Access-Memories: Concept, Fault Models and Tests, S.Hamdioui et al., 2002],
        # so we use the state in cell_history "cell_history[-fp.SenOpsNum][fp.aCell]" directly to check
        # the sensitization conditions.
        if (fp.CFdsFlag == 1) and (cell_state[fp.vCell] == fp.vInit) and (op_seq == fp.Sen) and (
                cell_history[-fp.SenOpsNum][fp.aCell] == fp.aInit) and (visiting_cell == fp.aCell):
            cell_state_temp[fp.vCell] = fp.vFault
            return UPDATED
    else:
        if 'w' in op:
            if (fp.CFdsFlag == 0) and ((cell_state[fp.aCell] == fp.aInit) or (fp.aInit == '-')) and (
                    op_seq == fp.Sen) and (
                    cell_history[-fp.SenOpsNum][fp.vCell] == fp.vInit):
                cell_state_temp[fp.vCell] = fp.vFault
                return UPDATED
        elif 'r' in op:
            if (fp.CFdsFlag == 0) and ((cell_state[fp.aCell] == fp.aInit) or (fp.aInit == '-')) and (
                    op_seq == fp.Sen) and (
                    cell_history[-fp.SenOpsNum][fp.vCell] == fp.vInit):
                # check rdflag, if it is read-destruction fault, it is detected immediately
                if (fp.rdFlag == 1) or (fp.rdFlag == 2):
                    return DETECTED
                else:
                    cell_state_temp[fp.vCell] = fp.vFault
                    return UPDATED
            elif cell_state[fp.vCell] != op[1]:
                return DETECTED

    return UPDATE_ERROR


def apply_March_element(cell_state, traverse_cell_order, OPS, FP1, FP2):

    # cell_history is a list contains multiple dictionaries.
    # it is used for saving the history states of the cells after each operation in several operation steps.
    # the length of cell_history depends on the longest sensitization sequence of current faults.
    # updating cell_history is by popping out the first dict and append the newest dict.
    cell_history = []

    def set_op_seq_history(fp, seq):
        nonlocal op_seq_history
        if fp is FP1:
            op_seq_history[0] = seq
        if fp is FP2:
            op_seq_history[1] = seq
        return

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
        # clear op_seq_history at the beginning of a March element apply
        op_seq_history = ['', '']

        for op_num, op in enumerate(ops, 1):
            # cell_state_temp is a copy of cell_state, it accepts the value updates during operations,
            # and update the final result to cell_state after an operation is applied.
            cell_state_temp = copy.deepcopy(cell_state)

            # a-cell case
            if visiting_cell != 'v':
                # write operation case
                if 'w' in op:
                    # Sensitization operation
                    if op_num >= FP1.SenOpsNum:
                        op_seq = generate_operation_sequence(FP1.SenOpsNum, op_num, ops)
                        update_fault_state(FP1, cell_state, cell_state_temp, cell_history, visiting_cell, op_seq, op)
                        set_op_seq_history(FP1, op_seq)

                    if FP2 is not FP1:
                        if op_num >= FP2.SenOpsNum:
                            op_seq = generate_operation_sequence(FP2.SenOpsNum, op_num, ops)
                            update_fault_state(FP2, cell_state, cell_state_temp, cell_history, visiting_cell, op_seq, op)
                            set_op_seq_history(FP2, op_seq)

                    # for write operation on a-cells, a-cell state will be changed
                    # no matter whether the fault is sensitized on v-cell
                    cell_state_temp[visiting_cell] = op[1]
                    # update the values in cell_state at the end of the operation
                    cell_state.update(cell_state_temp)
                    set_cell_history()

                # read operation case
                elif 'r' in op:
                    if cell_state[visiting_cell] != op[1]:
                        return DETECTED
                    else:
                        if op_num >= FP1.SenOpsNum:
                            op_seq = generate_operation_sequence(FP1.SenOpsNum, op_num, ops)
                            update_flag[0] = update_fault_state(FP1, cell_state, cell_state_temp, cell_history,
                                                                visiting_cell, op_seq, op)
                            set_op_seq_history(FP1, op_seq)
                            if update_flag[0] == DETECTED:
                                return DETECTED

                        if FP2 is not FP1:
                            if op_num >= FP2.SenOpsNum:
                                op_seq = generate_operation_sequence(FP2.SenOpsNum, op_num, ops)
                                update_flag[1] = update_fault_state(FP2, cell_state, cell_state_temp, cell_history,
                                                                    visiting_cell, op_seq, op)
                                set_op_seq_history(FP2, op_seq)
                                if update_flag[1] == DETECTED:
                                    return DETECTED

                        cell_state.update(cell_state_temp)
                        set_cell_history()

            # v-cell case
            else:
                # write operation case
                if 'w' in op:
                    # Sensitization operation
                    if op_num >= FP1.SenOpsNum:
                        op_seq = generate_operation_sequence(FP1.SenOpsNum, op_num, ops)
                        update_flag[0] = update_fault_state(FP1, cell_state, cell_state_temp, cell_history,
                                                            visiting_cell, op_seq, op)
                        set_op_seq_history(FP1, op_seq)
                    #    if update_state == DETECTED:
                    #        return DETECTED
                    else:
                        update_flag[0] = UPDATE_ERROR

                    if FP2 is not FP1:
                        if op_num >= FP2.SenOpsNum:
                            op_seq = generate_operation_sequence(FP2.SenOpsNum, op_num, ops)
                            update_flag[1] = update_fault_state(FP2, cell_state, cell_state_temp, cell_history,
                                                                visiting_cell, op_seq, op)
                            set_op_seq_history(FP2, op_seq)
                        #    if update_state == DETECTED:
                        #        return DETECTED
                    else:
                        update_flag[1] = UPDATE_ERROR

                    # for write operation on v-cell, if the fault is sensitized,
                    # the final state of v-cell is the fault state, and the write operation will be suppressed
                    if UPDATED not in update_flag:
                        cell_state_temp[visiting_cell] = op[1]

                    cell_state.update(cell_state_temp)
                    set_cell_history()

                # read operation case
                elif 'r' in op:
                    # Sensitization operation
                    if op_num >= FP1.SenOpsNum:
                        op_seq = generate_operation_sequence(FP1.SenOpsNum, op_num, ops)
                        update_flag[0] = update_fault_state(FP1, cell_state, cell_state_temp, cell_history,
                                                            visiting_cell, op_seq, op)

                        if update_flag[0] == DETECTED:
                            return DETECTED
                        # special case for drd
                        elif (FP1.rdFlag == -1) and (update_flag[0] == UPDATED) and (op_seq_history[0] == FP1.Sen) \
                                and (cell_state_temp[visiting_cell] != op[1]):
                            return DETECTED

                        set_op_seq_history(FP1, op_seq)
                    elif cell_state_temp[visiting_cell] != op[1]:
                        return DETECTED

                    if FP2 is not FP1:
                        if op_num >= FP2.SenOpsNum:
                            op_seq = generate_operation_sequence(FP2.SenOpsNum, op_num, ops)
                            update_flag[1] = update_fault_state(FP2, cell_state, cell_state_temp, cell_history,
                                                                visiting_cell, op_seq, op)

                            if update_flag[1] == DETECTED:
                                return DETECTED
                            # special case for drd
                            elif (FP1.rdFlag == -1) and (update_flag[1] == UPDATED) and (op_seq_history[1] == FP2.Sen) \
                                    and (cell_state_temp[visiting_cell] != op[1]):
                                return DETECTED

                            set_op_seq_history(FP2, op_seq)
                        elif cell_state_temp[visiting_cell] != op[1]:
                            return DETECTED

                    cell_state.update(cell_state_temp)
                    set_cell_history()

    return APPLIED


def eval_2comp(FP1, FP2, March, mode, logfile):
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

            logfile.write("  evaluating element \"%s\" under %s\n" % (element, str(traverse_cell_order)))
            ops = element.split(',')
            if apply_March_element(cell_state, traverse_cell_order, ops, FP1, FP2) == DETECTED:
                eval_flag.append('success')
                logfile.write("    current fault is detected.\n")
                break
            else:
                logfile.write("    current fault is NOT detected.\n")
                pass

    if eval_flag.count('success') == len(cell_order):
        return DETECTED
    else:
        return UNDETECTED

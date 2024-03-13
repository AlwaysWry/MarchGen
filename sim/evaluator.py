import copy

# return state definitions
INIT_ERROR = 0
INIT_OK = 1
ERROR = 0
UPDATED = 1
UPDATE_ERROR = 0
APPLIED = 1
DETECTED = 10
UNDETECTED = 100
PROFOUND = 1
FAST = 0


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


def init_cell_state(march, cell_state):
    cell_state.clear()

    if len(march) == 0:
        return INIT_ERROR

    if len(march[0].split(',')) == 2:
        # cell_state is a dictionary, which contains current states of a1, a2, v respectively
        if '0' in march[0]:
            cell_state['a1'] = '0'
            cell_state['a2'] = '0'
            cell_state['v'] = '0'
        elif '1' in march[0]:
            cell_state['a1'] = '1'
            cell_state['a2'] = '1'
            cell_state['v'] = '1'
        else:
            return INIT_ERROR
        # print(cell_state)
        return INIT_OK
    else:
        return INIT_ERROR


def init_cell_snapshot(cell_snapshot, cell_state, FP1, FP2):
    # cell_snapshot is a dictionary contains multiple lists.
    # it is used for saving the history states of each cell after each operation in several operation steps.
    # the length of lists in cell_snapshot depends on the longest sensitization sequence of current faults.
    # cell_snapshot is initialized as the cell states after the first write operation.
    cell_snapshot.clear()
    cell_snapshot['a1'] = []
    cell_snapshot['a2'] = []
    cell_snapshot['v'] = []

    for i in range(max(FP1.SenOpsNum, FP2.SenOpsNum)):
        cell_snapshot['a1'].append(cell_state['a1'])
        cell_snapshot['a2'].append(cell_state['a2'])
        cell_snapshot['v'].append(cell_state['v'])

    return


def update_cell_snapshot(cell_snapshot, cell_state, visiting_cell, update_flag, FP1, FP2):
    # cell_snapshot is updated when the state of the visiting cell is updated
    if (visiting_cell != 'v') and (UPDATED in update_flag):
        relevant_cells = [visiting_cell, 'v']
    else:
        relevant_cells = [visiting_cell]
    for cell in relevant_cells:
        cell_snapshot[cell].append(cell_state[cell])
        if len(cell_snapshot[cell]) > max(FP1.SenOpsNum, FP2.SenOpsNum):
            del cell_snapshot[cell][0]
    return


def init_operation_snapshot(op_snapshot):
    # op_snapshot is for recording all applied operations to each cell, and a counter for recording the total number
    # operations that each cell have been applied
    op_snapshot.clear()
    op_snapshot['a1'] = {'op_counter': 0, 'op_history': []}
    op_snapshot['a2'] = {'op_counter': 0, 'op_history': []}
    op_snapshot['v'] = {'op_counter': 0, 'op_history': []}

    return


def update_operation_history(op_snapshot, visiting_cell, FP1, FP2, op_seq_history):
    for fp_index, fp_item in enumerate([FP1, FP2]):
        op_seq = ''
        for op in op_snapshot[visiting_cell]['op_history'][-fp_item.SenOpsNum:]:
            op_seq = op_seq + op
        op_seq_history[fp_index] = op_seq

    return


def update_operation_snapshot(op_snapshot, visiting_cell, op, FP1, FP2):
    op_snapshot[visiting_cell]['op_history'].append(op)
    op_snapshot[visiting_cell]['op_counter'] = op_snapshot[visiting_cell]['op_counter'] + 1
    if op_snapshot[visiting_cell]['op_counter'] > max(FP1.SenOpsNum, FP2.SenOpsNum):
        del op_snapshot[visiting_cell]['op_history'][0]

    return


def get_relevant_seq(SenOpNum, snapshot, visiting_cell):
    op_seq = ''
    op_num = len(snapshot[visiting_cell]['op_history'])
    for i in range(SenOpNum - 1, -1, -1):
        op_seq = op_seq + snapshot[visiting_cell]['op_history'][op_num - i - 1]

    return op_seq


def update_fault_state(fp, cell_state, cell_state_temp, cell_snapshot, visiting_cell, op_seq, op):
    if visiting_cell != 'v':
        # for dynamic faults, the init state of a-cell/v-cell should be kept in cell_snapshot, or it will be
        # covered by new state. The init state also cannot be destroyed by fault state of other fault in the
        # middle of sensitization, like <0w1/0/->*<0w1w0/1/->, the FP2 can be sensitized by "O0, w1, w0",
        # because processing the fault means that it can be sensitized under correct conditions, i.e. "dynamic
        # faulty behavior can take place in the absence of static faults."
        # [Dynamic Faults in Random-Access-Memories: Concept, Fault Models and Tests, S.Hamdioui et al., 2002],
        # so we use the state in cell_snapshot "cell_snapshot[fp.aCell][-fp.SenOpsNum]" directly to check
        # the sensitization conditions.
        if (fp.CFdsFlag == 1) and (cell_state[fp.vCell] == fp.vInit) and (op_seq == fp.Sen) and (
                cell_snapshot[fp.aCell][-fp.SenOpsNum] == fp.aInit) and (visiting_cell == fp.aCell):
            cell_state_temp[fp.vCell] = fp.vFault
            return UPDATED
    else:
        if 'w' in op:
            if (fp.CFdsFlag == 0) and ((cell_state[fp.aCell] == fp.aInit) or (fp.aInit == '-')) and (
                    op_seq == fp.Sen) and (
                    cell_snapshot[fp.vCell][-fp.SenOpsNum] == fp.vInit):
                cell_state_temp[fp.vCell] = fp.vFault
                return UPDATED
        elif 'r' in op:
            if (fp.CFdsFlag == 0) and ((cell_state[fp.aCell] == fp.aInit) or (fp.aInit == '-')) and (
                    op_seq == fp.Sen) and (
                    cell_snapshot[fp.vCell][-fp.SenOpsNum] == fp.vInit):
                # check rdflag, if it is read-destruction fault, it is detected immediately
                if (fp.rdFlag == 1) or (fp.rdFlag == 2):
                    return DETECTED
                else:
                    cell_state_temp[fp.vCell] = fp.vFault
                    return UPDATED
            elif cell_state[fp.vCell] != op[1]:
                return DETECTED

    return UPDATE_ERROR


def apply_March_element(cell_state, cell_snapshot, op_snapshot, traverse_cell_order, OPS, FP1, FP2):

    for visiting_cell in traverse_cell_order:
        # operations of each element, index starts from 1
        ops = OPS[1:]
        update_flag = ['', '']
        # clear op_seq_history at the beginning of a March element apply
        op_seq_history = ['', '']

        for op_index, op in enumerate(ops):
            # cell_state_temp is a copy of cell_state, it accepts the value updates during operations,
            # and update the final result to cell_state after an operation is applied.
            cell_state_temp = copy.deepcopy(cell_state)
            # cell_state_temp = cell_state.copy()

            # update operation history and snapshot when a new operation is applied
            update_operation_history(op_snapshot, visiting_cell, FP1, FP2, op_seq_history)
            update_operation_snapshot(op_snapshot, visiting_cell, op, FP1, FP2)

            # a-cell case
            if visiting_cell != 'v':
                # write operation case
                if 'w' in op:
                    # Sensitization operation
                    if op_snapshot[visiting_cell]['op_counter'] >= FP1.SenOpsNum:
                        op_seq = get_relevant_seq(FP1.SenOpsNum, op_snapshot, visiting_cell)
                        update_flag[0] = update_fault_state(FP1, cell_state, cell_state_temp, cell_snapshot,
                                                            visiting_cell, op_seq, op)

                    if FP2 is not FP1:
                        if op_snapshot[visiting_cell]['op_counter'] >= FP2.SenOpsNum:
                            op_seq = get_relevant_seq(FP2.SenOpsNum, op_snapshot, visiting_cell)
                            update_flag[1] = update_fault_state(FP2, cell_state, cell_state_temp,
                                                                cell_snapshot, visiting_cell, op_seq, op)

                    # for write operation on a-cells, a-cell state will be changed
                    # no matter whether the fault is sensitized on v-cell
                    cell_state_temp[visiting_cell] = op[1]
                    # update the values in cell_state and cell_snapshot at the end of the operation
                    cell_state.update(cell_state_temp)
                    update_cell_snapshot(cell_snapshot, cell_state, visiting_cell, update_flag, FP1, FP2)

                # read operation case
                elif 'r' in op:
                    if cell_state[visiting_cell] != op[1]:
                        return DETECTED, op_index
                    else:
                        if op_snapshot[visiting_cell]['op_counter'] >= FP1.SenOpsNum:
                            op_seq = get_relevant_seq(FP1.SenOpsNum, op_snapshot, visiting_cell)
                            update_flag[0] = update_fault_state(FP1, cell_state, cell_state_temp, cell_snapshot,
                                                                visiting_cell, op_seq, op)
                            if update_flag[0] == DETECTED:
                                return DETECTED, op_index

                        if FP2 is not FP1:
                            if op_snapshot[visiting_cell]['op_counter'] >= FP2.SenOpsNum:
                                op_seq = get_relevant_seq(FP2.SenOpsNum, op_snapshot, visiting_cell)
                                update_flag[1] = update_fault_state(FP2, cell_state, cell_state_temp, cell_snapshot,
                                                                    visiting_cell, op_seq, op)
                                if update_flag[1] == DETECTED:
                                    return DETECTED, op_index

                        cell_state.update(cell_state_temp)
                        update_cell_snapshot(cell_snapshot, cell_state, visiting_cell, update_flag, FP1, FP2)

                else:
                    print("Illegal March operation found, program is terminated.\n")
                    print("Please check the March test file.\n")
                    return ERROR,

            # v-cell case
            else:
                # write operation case
                if 'w' in op:
                    # Sensitization operation
                    if op_snapshot[visiting_cell]['op_counter'] >= FP1.SenOpsNum:
                        op_seq = get_relevant_seq(FP1.SenOpsNum, op_snapshot, visiting_cell)
                        update_flag[0] = update_fault_state(FP1, cell_state, cell_state_temp, cell_snapshot,
                                                            visiting_cell, op_seq, op)
                    else:
                        update_flag[0] = UPDATE_ERROR

                    if FP2 is not FP1:
                        if op_snapshot[visiting_cell]['op_counter'] >= FP2.SenOpsNum:
                            op_seq = get_relevant_seq(FP2.SenOpsNum, op_snapshot, visiting_cell)
                            update_flag[1] = update_fault_state(FP2, cell_state, cell_state_temp, cell_snapshot,
                                                                visiting_cell, op_seq, op)
                    else:
                        update_flag[1] = UPDATE_ERROR

                    # for write operation on v-cell, if the fault is sensitized,
                    # the final state of v-cell is the fault state, and the write operation will be suppressed
                    if UPDATED not in update_flag:
                        cell_state_temp[visiting_cell] = op[1]

                    cell_state.update(cell_state_temp)
                    update_cell_snapshot(cell_snapshot, cell_state, visiting_cell, update_flag, FP1, FP2)

                # read operation case
                elif 'r' in op:
                    # Sensitization operation
                    if op_snapshot[visiting_cell]['op_counter'] >= FP1.SenOpsNum:
                        op_seq = get_relevant_seq(FP1.SenOpsNum, op_snapshot, visiting_cell)
                        update_flag[0] = update_fault_state(FP1, cell_state, cell_state_temp, cell_snapshot,
                                                            visiting_cell, op_seq, op)

                        if update_flag[0] == DETECTED:
                            return DETECTED, op_index
                        # special case for drd
                        elif (FP1.rdFlag == -1) and (update_flag[0] == UPDATED) and (op_seq_history[0] == FP1.Sen) \
                                and (cell_state_temp[visiting_cell] != op[1]):
                            return DETECTED, op_index

                    elif cell_state_temp[visiting_cell] != op[1]:
                        return DETECTED, op_index

                    if FP2 is not FP1:
                        if op_snapshot[visiting_cell]['op_counter'] >= FP2.SenOpsNum:
                            op_seq = get_relevant_seq(FP2.SenOpsNum, op_snapshot, visiting_cell)
                            update_flag[1] = update_fault_state(FP2, cell_state, cell_state_temp, cell_snapshot,
                                                                visiting_cell, op_seq, op)

                            if update_flag[1] == DETECTED:
                                return DETECTED, op_index
                            # special case for drd
                            elif (FP2.rdFlag == -1) and (update_flag[1] == UPDATED) and (op_seq_history[1] == FP2.Sen) \
                                    and (cell_state_temp[visiting_cell] != op[1]):
                                return DETECTED, op_index

                        elif cell_state_temp[visiting_cell] != op[1]:
                            return DETECTED, op_index

                    cell_state.update(cell_state_temp)
                    update_cell_snapshot(cell_snapshot, cell_state, visiting_cell, update_flag, FP1, FP2)

                else:
                    print("Illegal March operation found, program is terminated.\n")
                    print("Please check the March test file.\n")
                    return ERROR,

    return APPLIED,


def eval_2comp(FP1, FP2, march, mode):
    eval_flag = []
    eval_record = dict()
    cell_order = init_cell_order(mode)
    cell_state = {}
    cell_snapshot = {}
    op_snapshot = {}

    for order_index in range(len(cell_order)):

        if init_cell_state(march, cell_state) == INIT_ERROR:
            print("Illegal March test!\n")
            return ERROR

        init_cell_snapshot(cell_snapshot, cell_state, FP1, FP2)
        init_operation_snapshot(op_snapshot)

        order_key = str(cell_order[order_index]['up'])
        eval_record[order_key] = dict()

        # traverse march elements
        for m, element in enumerate(march[1:]):
            # set cell traverse order
            if 'down' in element:
                traverse_cell_order = cell_order[order_index]['down']
            else:
                traverse_cell_order = cell_order[order_index]['up']

            eval_record[order_key][element] = []
            ops = element.split(',')

            apply_result = apply_March_element(cell_state, cell_snapshot, op_snapshot, traverse_cell_order, ops, FP1, FP2)
            apply_flag = apply_result[0]

            if apply_flag == ERROR:
                return ERROR
            elif apply_flag == DETECTED:
                eval_flag.append('success')
                eval_record[order_key][element].extend(['success', apply_result[1]])
                break
            else:
                eval_record[order_key][element].extend(['fail', -1])
                pass

    if eval_flag.count('success') == len(cell_order):
        return DETECTED, eval_record
    else:
        return UNDETECTED, eval_record

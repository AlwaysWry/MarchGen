import os
import sys

sys.path.append('mwvc_solver/QUICK_VC')

import QUICK_VC_Py


def quickVC_solver(graph_file):
	arg_list = [graph_file, 'results/mwvclog', '10', '0', '0']
	QUICK_VC_Py.MWVC(arg_list)


if __name__ == '__main__':
	os.chdir("../")
	quickVC_solver('results/degenerated_2cF_graph.m2c')

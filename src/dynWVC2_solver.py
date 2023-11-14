import os
from ctypes import *


def dynWVC2_solver(graph_file):
	dyn_WVC2 = CDLL("mwvc_solver/DynWVC2/libDynWVC2.dll", winmode=0)

	graph = c_char_p(graph_file.encode())
	result = c_char_p(b'results/mwvclog')
	seed = c_char_p(b'37')
	cutoff_time = c_char_p(b'10')
	mode = c_char_p(b'0')

	dyn_WVC2.MWVC(graph, result, seed, cutoff_time, mode)
	return


if __name__ == '__main__':
	os.chdir("../")
	dynWVC2_solver('results/unlinked_2cF_graph.m2c')

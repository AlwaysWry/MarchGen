from ctypes import *


def dynWVC2_solver(graph_file):
	DynWVC2 = cdll.LoadLibrary("../mwvc_solver/DynWVC2/libDynWVC2.dll")

	graph = c_char_p(graph_file.encode())
	result = c_char_p(b'../results/mwvc.log')
	seed = c_char_p(b'37')
	cutoff_time = c_char_p(b'10')
	mode = c_char_p(b'0')

	DynWVC2.MWVC(graph, result, seed, cutoff_time, mode)


if __name__ == '__main__':
	dynWVC2_solver('../resources/unlinked_2cF.m2c')

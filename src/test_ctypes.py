from ctypes import *
import sys
sys.path.append('../mwvc_solver/QUICK_VC')
import QUICK_VC_Py

#lib = cdll.LoadLibrary("./DynWVC2.dll")
## lib.DynWVC2.argtypes = [c_int, c_wchar_p, c_wchar_p, c_wchar_p, c_wchar_p]
## lib.DynWVC2.restype = c_int
#arg_num = c_int(4)
#file = c_char_p(b'E:\\CLproj\\March2Comp\\resources\\unlinked_2cF.m2c')
#seed = c_char_p(b'37')
#cutoff_time = c_char_p(b'10')
#mode = c_char_p(b'0')
#
#lib.DynWVC2(arg_num, file, seed, cutoff_time, mode)
arg_list = ['../resources/unlinked_2cF.m2c', 'bio.txt', '10', '0', '0']
QUICK_VC_Py.MWVC(arg_list)

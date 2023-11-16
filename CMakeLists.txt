CMAKE_MINIMUM_REQUIRED(VERSION 3.24)
project(March2Comp)

if (CMAKE_SYSTEM_NAME MATCHES "Linux")

    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    set(CMAKE_CXX_FLAGS "-fopenmp")
    find_package(OpenMP)

    include_directories(mwvc_solver/QUICK_VC/include/mwvc)

    # settings for pybind11
    set(PYTHON_EXECUTABLE /home/zsr/.conda/envs/March2Comp/bin/python)
    # the directory of pybind11 needs to be changed before compiling
    set(pybind11_DIR "/home/zsr/.conda/envs/March2Comp/lib/python3.11/site-packages/pybind11/share/cmake/pybind11")
    find_package(pybind11 REQUIRED)

    set(LIBRARY_OUTPUT_PATH ../mwvc_solver/QUICK_VC)
    pybind11_add_module(QUICK_VC_Py mwvc_solver/QUICK_VC/QUICK_VC_Py.cpp)


elseif (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set(CMAKE_CXX_STANDARD 17)

    include_directories(mwvc_solver/DynWVC2)
    set(LIBRARY_OUTPUT_PATH mwvc_solver/DynWVC2)

    add_library(DynWVC2 SHARED mwvc_solver/DynWVC2/mwvc.cpp)

endif ()

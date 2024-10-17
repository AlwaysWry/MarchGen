# March2Comp: A Structure-oriented March Sequence Generator for Multi-operation Dynamic 2-composite Faults in Memory

  March2Comp is an automatic March test generator, aiming to generate March tests for 2-composite faults with arbitrary number of sensitization operations. It also integrates Sim2Comp, a fault simulator to check the fault coverage of a given March test under different 2-composite fault models.
    
  File directory structure:

  March2Comp
<br>    ├── common // basic modules, includes fault list parser, fault list generator, etc.
<br>    ├── mwvc_solver // solver of minimum weight vertex coverage (MWVC) problem
<br>    ├── resources // example fault lists and march tests
<br>    ├── sim // fault simulator
<br>    ├── src // source files of March test generator
<br>    ├── CMakeLists.txt
<br>    ├── March2Comp // start shell script
<br>    ├── March2Comp.py // entrance program of March2Comp

## Requirements
gcc/g++ 9.3.1 or higher
<br>  Python 3.11 or higher
<br>  pybind11
<br> CMake 3.24 or higher
## Users' Guide
Please change ```pybind11_DIR``` and ```PYTHON_EXECUTABLE``` in ```CMakeLists.txt``` to your path of pybind11 library and Python3 before compilation.
### Compile and run on Linux
A shell script is used to implement compilation and running of March2Comp. 

Commands for compilation:
<br>    ```March2Comp compile```

Commands for using March test generator:
<br>    ```March2Comp gen <fault list name> [fault model name=2cF_3, 2cF_2aa] [--nomp]```
*    Fault model name is optional, default value is '2cF_2aa'. Under 2cF_3 model, unrealistic faults (but realistic under 2cF_2aa model) will be removed from the list first.
*    "--nomp" is an option to close the multi-process simulation, which is open in default.

Commands for using fault simulator independently:
<br>    ```March2Comp sim <march test filename> <fault list filename> [--nomp]```

### Compile and run on Windows
You can use your own compiler and python interpreter to compile and run March2Comp. What follows is an example using CMake and MinGW in PowerShell or cmd (recommend running with administrator privilege):
<br>    ```cmake -G "MinGW Makefiles"```
<br>    ```mingw32-make```

Using March test generator:
<br>    ```python March2Comp.py <fault list file> [fault model name=2cF_3, 2cF_2aa] [--nomp]```

Using fault simulator independently:
<br>    ```python ./sim/sim2comp.py <march test filename> <fault list filename> [--nomp]```

### Instructions
1. Standard format of fault lists: a fault list consist of faults in standard fault primitive (FP) format: <br> ```<Sa;Sv/F/R>```
<br> ```<Sa;Sv/F/R>```
<br> ```<Sa;Sv/F/R>```
<br> ```...``` <br> where: <br> ```Sa```: sensitization operations on aggressor cell <br> ```Sv```: sensitization operations on victim cell <br> ```F```: faulty value of victim cell <br> ```R```: read-out value of victim cell

2. Standard format of March tests: a March test includes a series of March elements (ME). An ME consists of address order (AO) and operation sequence (OS). Here's an example:
<br> ```up,r0,w0,w1,w0,w1,r1``` <br> Possible AOs include ```up```, ```down``` and ```any```. <br> Possible operations in OS include ```r0```, ```w0```, ```r1``` and ```w1```. <br> There are some examples of March tests in ```resources/march_tests```.

3. Report of March test generator, fault simulator and test log are in ```generation_report.txt```, ```simulation_report.txt``` and ```testlog``` in ```results/``` respectively.

### Contact alwayswry@stu.pku.edu.cn if you have any problem.


 

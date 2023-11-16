# March2Comp: A Structure-oriented March Test Generator for Static and Multi-operation Dynamic 2-composite Faults in Memory

  It is a demo Python project for March2Comp.
    
  File dictionary structure:

  March2Comp
<br>    ├── basic // basic modules, includes fault list parser, fault list generator, etc.
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
### Compile and run on Linux
A shell script is used to implement compilation and running of March2Comp. 

Commands for compilation:
<br>    ```March2Comp compile```

Commands for using March test generator:
<br>    ```March2Comp gen <fault list name>```

Commands for using fault simulator independently:
<br>    ```March2Comp sim <march test filename> <fault list filename>```

### Compile and run on Windows
You can use your own compiler and python intepreter to compile and run March2Comp. What follows is an example using CMake and MinGW:

<br>    ```cmake -G "MinGW makefiles"```
<br>    ```mingw32-make```

Using March test generator:
<br>    ```python March2Comp.py <fault list file>```

Using fault simulator independently:
<br>    ```python ./sim/sim2comp.py <march test filename> <fault list filename>```


 

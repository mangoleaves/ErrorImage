# README

## How to build

1. Install Boost library (at least 1.78 version).
2. Use CMake to configure CMakeLists.txt.
3. Build all targets.

## How to run

1. The name of executable file is "exeGCLF-ImageTriSimp", put it with Qt5Core.dll and Qt5Gui.dll in the same directory.
2. "exeGCLF-ImageTriSimp" can be executed in three modes:
    * Given an initial triangulation, optimize the triangulation to be error-bounded (see example/error_bounded.bat).
    * Given an linear error-bounded triangulation, simplify it (see example/linear_simplify.bat).
    * Given an linear error-bounded triangulation, make it curved and simplify it (see example/curved_simplify.bat).
3. The color function and other parameters can be changed in "config.json".

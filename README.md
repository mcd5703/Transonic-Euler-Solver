---
title: Transonic Euler Solver in C++
author: Melik Demirel & Nicholas Sernberger
email: mcd5703@psu.edu, nms6249@psu.edu
date: 12/12/24
---

# Transonic Euler Solver in C++

This repository contains all files for the Transonic Euler Solver in C++. 

Created by Nick Sernberger and Melik Demirel.

Credits: Dr. Simon Miller and Dr. James Coder.

This project has primarily been tested on MacOS with `cmake`.

Capabilities are offered but not tested for:
- RHEL8 using hte ROAR cluster
- When we use `git` console in Windows or terminal in Linux.

## Dependencies

1. [Gnuplot](http://www.gnuplot.info)
   1. MacOS -- use `brew install gnuplot; brew install qt`
   2. Linux -- use `sudo apt-get install gnuplot`
   3. Windows -- install binaries at [url](http://www.gnuplot.info/download.html)
2. [GLEW](https://glew.sourceforge.net)
   1. MacOS -- can also use a system-level installation as `brew install glew`
   2. Linux -- system-level installation as `sudo apt-get install libglew-dev`
   3. Windows -- system-level installation as `https://sourceforge.net/projects/glew/files/glew/2.1.0/glew-2.1.0-win32.zip/download`
3. [GLFW](https://www.glfw.org)
4. [Matplot++](https://alandefreitas.github.io/matplotplusplus/)
5. [doxygen](https://doxygen.nl) for automated documentation html/latex
   1. [graphviz](https://graphviz.org) for inheritance graph drawing
   2. TeX Distribution for epstodpdf
      1. [MiKTex](https://miktex.org) for Window
      2. [MacTeX](https://tug.org/mactex/) for MacOS (MiKTeX hardly works with Sequoia)
   3. [ghostscript](https://www.ghostscript.com) for rendering doxygen

On MacOS, these will get you most of the way:
```zsh
brew install gnuplot
brew install qt
brew install glew
brew install graphviz
brew install ghostscript
brew install mactex
brew install doxygen
```

If using the [ROAR](rcportal.hpc.psu.edu/) system, be sure to load the necessary modules:

```console
module load cmake
module load gcc/9.1.0
module load gnuplot
```

And, if on MSYS UCRT64:

```console
pacman -S --needed mingw-w64-ucrt-x86_64-blas
pacman -S --needed mingw-w64-ucrt-x86_64-lapack
pacman -S --needed mingw-w64-ucrt-x86_64-vtk
pacman -S --needed mingw-w64-x86_64-libjpeg-turbo
pacman -S --needed mingw-w64-x86_64-libtiff
pacman -S --needed mingw-w64-x86_64-zlib
pacman -S --needed mingw-w64-x86_64-gnuplot
```

## Structure of the Code

After you build the whole program, you will get several directories.

```console
|--- build       <- where the build files and executables are stored
|--- include     <- header files
|--- src         <- main source files
|--- main        <- test code that has a `main` function -- one at a time for `add_executable`
|--- results     <- a folder to contain figures and stderr/stdout
|--- third_party <- to include third_party dependencies and build/install the these dependencies
|--- docs        <- documentation folder (open in html)
|--- img         <- documentation image icon
```

## How to Build

Just run the `build.sh`  `bash` file in the git console terminal in Windows or the default `terminal` in MacOS or Linux. Before running this script, make sure you have installed `git`, `cmake`, and related build toolchains. Run the build file that matches your system...

```bash
bash build.sh
```

In this shell script, we will download, compile and install the third_party dependencies (GLEW, glfw, Matplot++).

Build Script Testing:
* `build_linux.sh` -- tested on RCPortal ROAR
* `build_msvc.sh` -- tested on WIN10
* `build_clang.sh` -- tested on macOS Sequoia 15.1.1 with `homebrew`

#### Some notes on this:

The following directories should appear automatically when `bash build.sh` is run in the terminal of the project.
- `../build/`
- `../results/`
- `../third_party/`
- `../docs/`

ONLY run executables that are built inside the `../build/` folder after running `bash build.sh`. Do not run files in `../src` for they will not compile properly. You should not have to run any `.cpp` files if the `build.sh` file is run successfully and the executable file is successfully generated.

An alternative to `bash build.sh` is to run:

    chmod +x build.sh
    ./build.sh

## How to Run the Final Program

All the executables will be stored under the build folder.

Use the following command to run the executable.

```bash
echo "running the executable with some default parameters"
echo "./build/main -c config.inp > results.txt 2>&1"
echo "  the 2>&1 redirects the stderr to a 1 so we don't see the gnuplot problems"
mkdir -p ${project_root_dir}/results
./build/main -c config.inp > results/results.txt 2>&1
```

# Troubleshooting

Chances are that a `GLEW` is going to be a pain the ass to get working this way. In this method, we are compiling it all from scratch that, quite frankly, is not a great idea! This is because we can easily get a pre-compiled binary of it that aleviates all of our troubles as:

* Linux:

```bash
sudo apt-get install libglew-dev
```

* MacOS

```
brew reinstall glew
brew link glew
```

* Windows  
Install from the pre-compiled binaries provided by [GLEW](https://glew.sourceforge.net)
  * update CMakeLists.txt

```cmake
add_definitions(-DGLEW_STATIC) # static linking for glew only
set(GL_STATIC_LIBRARIES "${CMAKE_SOURCE_DIR}/lib/glew32sd.lib;${CMAKE_SOURCE_DIR}/lib/glfw3.lib")
...
target_link_libraries(main ${OPENGL_gl_LIBRARY} ${GL_STATIC_LIBRARIES})
```

## Extended System Requirements Guide for Troubleshooting on MacOS

If you are on MacOS, you **MUST** update your Mac and Xcode to the latest version, via Settings > Software Update and Apple Store > Xcode respectively. 
Only then will you be able to successfully follow the steps below.

Here's a step-by-step guide to obtain everything you need...

#### 1. VS Code
Download VS Code from the website https://code.visualstudio.com/.

You will also need (or want) the following extensions:
- C/C++ by Microsoft
- C/C++ Extension Pack by Microsoft
- C/C++ Themes by Microsoft
- CMake by twxs
- CMake Tools by Microsoft
- CodeLLDB by Vadim Chugunov
- indent-rainbow by oderwat
- Makefile Tools by Microsoft
- Pylance by Microsoft
- Python by Microsoft
- Python debugger by Microsoft

#### 2. Xcode Command Line Tools (xcode-tools): C++ Compiler
If you updated Xcode, Xcode may have deleted its Command Line Tools. 

Xcode Command Line Tools are required to run a C++ compiler and git on a mac, so you must reinstall Xcode's Commandline Tools. 

To do so, open terminal. Then, run the following: 

    xcode-select --install
    
You may also simply run:

    gcc 

Finally, check for installation with:

    gcc --version

Since git should have also reinstalled, you may check:

    git -v

#### 3. Homebrew: Package Management for Development Tool and Library Installation
To get Homebrew, go to the webpage (https://brew.sh/), and copy the command. Paste it in terminal and run:

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

Enter your password when asked.

Then, run:

    (echo; echo 'eval "$(/opt/homebrew/bin/brew shellenv)"') >> ~/.zprofile
    eval "$(/opt/homebrew/bin/brew shellenv)"

#### 4. Homebrew Packages
Now that you have Homebrew, you will need to install the following packages:

    brew install cmake
    brew install python3
    brew install numpy
    brew install gnuplot
    brew install python-matplotlib

#### 5. Summary
Make sure you have everything:
- VS Code (download from the site)
- Homebrew package manger 
- cmake 
    - make
- pip3
    - python3
        - numpy
        - matplotlib
- xcode-tools (xcode command line tools)
    - OpenGL
    - clang, clang++, gcc, g++ or some other c++ compiler
    - git
- gnuplot

You may check if you have all these installed with:

    #!/bin/bash
    cmake --version
    make --version
    python3 --version
    pip3 --version
    python3 -c "import numpy; print(numpy.__version__)"
    xcode-select --version
    gnuplot --version
    git --version
    brew --version
    clang --version || clang++ --version || gcc --version || g++ --version

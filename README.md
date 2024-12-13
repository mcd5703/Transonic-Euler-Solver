---
title: Homework 02
author: Simon W Miller
email: swm154@psu.edu
date: today
---
# Homework 02

This project has been tested on MacOS with `cmake` and in RHEL8 using hte ROAR cluster. And it should work when use `git` console in Windows or terminal in Linux.

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
```

## How to Build

Just run the `build.sh`  `bash` file in the git console terminal in Windows or the default `terminal` in MacOS or Linux. Before running this script, make sure you have installed `git`, `cmake`, and related build toolchains.

```bash
bash build.sh
```

In this shell script, we will download, compile and install the third_party dependencies (GLEW, glfw, Matplot++).

Build Script Testing:
* `build_linux.sh` -- tested on RCPortal ROAR
* `build_msvc.sh` -- tested on WIN10
* `build_clang.sh` -- tested on macOS Sequoia 15.1.1 with `homebrew`

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

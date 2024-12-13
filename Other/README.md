# AERSP424HW-2
This repository contains all files for AERSP 424 Homework #2. 
Created by Nick Sernberger and Melik Demirel.

## Notes

Please view Problem4.md, Problem5.md, Problem6.md in `..\src\` through online github to see our responses and notes about Problem4:

- [Problem4.md](https://github.com/nicksern/AERSP424HW-2/blob/main/src/Problem4.md)
- [Problem5.md](https://github.com/nicksern/AERSP424HW-2/blob/main/src/Problem5.md)
- [Problem6.md](https://github.com/nicksern/AERSP424HW-2/blob/main/src/Problem6.md)

Since `.md` files allow MathJax, we've provided a rundown of the [mathematical theory in problem 4](https://github.com/nicksern/AERSP424HW-2/blob/main/src/Problem4.md#the-mathematical-theory).

You **MUST** exit out of problem 4.2 to continue on to problem 4.1 (they are in reverse order, due to some potential bugs with 4.1).

## Instructions

### Step 1.

Simply run the following command in terminal once you are inside the project folder.

    bash build.sh

#### Some notes on this:

The following directories should appear automatically when `bash build.sh` is run in the terminal of the project.
- `../build/`
- `../results/`
- `../third_party/`

ONLY run executables that are built inside the `../build/` folder after running `bash build.sh`. Do not run files in `../src` for they will not compile properly. You should not have to run any `.cpp` files if the `build.sh` file is run successfully and the executable file is successfully generated.

An alternative to `bash build.sh` is to run:

    chmod +x build.sh
    ./build.sh

### Step 2.
Run the executable called `homework02` in `../build/`.

Text results will be shown in `..\results\results.txt`. For plots, you may need to run `homework02` in `../build/` again. If gnuplot is correctly installed, then the graphical results should automatically show.

## System Requirements Guide for Troubleshooting

### For MacOS
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

### For Windows
Unfortunately, we do not have extensive support for Windows. However, requirements listed in MacOS can help.

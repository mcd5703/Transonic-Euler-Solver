#!/bin/bash
current_directory="$(cd "$(dirname "$0")" && pwd)"
project_root_dir=$current_directory

# if on ICDS ROAR you will need to load modules to make this work
# module load anaconda
# module load cmake
# rem module load gcc/9.1.0

# if on Windows and using MSYS, might need a couple of tools
# https://code.visualstudio.com/docs/cpp/config-mingw
# pacman -S cmake
# pacman -S git
# pacman -S --needed base-devel mingw-w64-ucrt-x86_64-toolchain
# pacman -S mingw-w64-ucrt-x86_64-python-numpy
# pacman -S mingw-w64-ucrt-x86_64-blas
# pacman -S mingw-w64-ucrt-x86_64-openblas
# pacman -S mingw-w64-ucrt-x86_64-vtk

# if on MacOS, you will need to install cmake, python, numpy, xcode-tools (command line tools), and gnuplot

echo "Hello Friend! Let's try and build this application together :D"

mkdir -p ${project_root_dir}/third_party
cd ${project_root_dir}/third_party
#compile MATPLOTPLUSPLUS
if [ ! -d "matplotplusplus" ]; then
    git clone https://github.com/alandefreitas/matplotplusplus
    cd ${project_root_dir}/third_party/matplotplusplus
    rm -rf build install
    mkdir -p build && mkdir -p install && cd build
    cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/matplotplusplus/install -DCMAKE_BUILD_TYPE=Release ..
    cmake --build . -j 4
    cmake --install .
fi

cd ${project_root_dir}/third_party
if [ ! -d "glfw" ]; then
    git clone https://github.com/glfw/glfw.git glfw
    cd ${project_root_dir}/third_party/glfw
    rm -rf build install
    mkdir -p build && mkdir -p install && cd build
    cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/glfw/install -DCMAKE_BUILD_TYPE=Release ..
    cmake --build . -j 4
    cmake --install .
fi

cd ${project_root_dir}/third_party
if [ ! -d "glm" ]; then
    git clone https://github.com/g-truc/glm
    cd ${project_root_dir}/third_party/glm
    rm -rf build install
    mkdir -p build && mkdir -p install && cd build
    cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/glm/install -DGLM_BUILD_TESTS:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=OFF ..
    cmake --build . -j 4
    cmake --install .
fi

# compile the rest of the application
cd ${project_root_dir}

# remove the build directory that has the current code in it
echo "deleting the BUILD directory"
rm -rf ${project_root_dir}/build
rm -rf ${project_root_dir}/install

echo "make a new BUILD directory to start the compiling process"
mkdir -p ${project_root_dir}/build
cd ${project_root_dir}/build

echo "cmake engage!"
cmake ..

echo "convert this to an executable application -- let's go!!"
cmake --build . -j 4
cmake --install .
cpack
cd ${project_root_dir}
echo "declare success -- hooray!"

echo "running the executable with some default parameters"
echo "./build/homework02 -c config.inp > results.txt 2>&1"
echo "  the 2>&1 redirects the stderr to a 1 so we don't see the gnuplot problems"
mkdir -p ${project_root_dir}/results
./build/homework02 -c config.inp > results/results.txt 2>&1

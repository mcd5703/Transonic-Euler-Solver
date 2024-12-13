#!/bin/sh
current_directory="$(cd "$(dirname "$0")" && pwd)"
project_root_dir=$current_directory
# export CC=gcc-14 CXX=g++-14
# export SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX15.1.sdk/

# if on MacOS, you will need to install cmake, python, numpy, xcode-tools, gnuplot, glew

echo "Hello Friend! Let's try and build this application together :D"

mkdir -p ${project_root_dir}/third_party
cd ${project_root_dir}/third_party
#compile MATPLOTPLUSPLUS
if [ ! -d "matplotplusplus" ]; then
    git clone https://github.com/alandefreitas/matplotplusplus
    cd ${project_root_dir}/third_party/matplotplusplus
    rm -rf build install
    mkdir -p build && mkdir -p install && cd build
    cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/matplotplusplus/install -DCMAKE_BUILD_TYPE=Release -DMATPLOTPP_BUILD_TESTS:BOOL=OFF -DMATPLOTPP_BUILD_EXAMPLES:BOOL=OFF -DMATPLOTPP_BUILD_WITH_SANITIZERS:BOOL=OFF ..
    cmake --build . --config Release -j 4
    cmake --install . --config Release
fi

cd ${project_root_dir}/third_party
if [ ! -d "glew" ]; then
    git clone https://github.com/Perlmint/glew-cmake.git glew
    cd ${project_root_dir}/third_party/glew
    #cd build
    #cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/glew/install -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_UTILS:BOOL=ON -DCMAKE_BUILD_TYPE=Debug ./cmake
    rm -rf build install
    mkdir -p build && mkdir -p install && cd build
    cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/glew/install -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_LIBDIR=lib -DBUILD_SHARED_LIBS:BOOL=ON ..
    cmake --build . -j 4 --config Release
    cmake --build . --config Release --target glew_s
    cmake --build . --config Release --target INSTALL
    cmake --install . --config Release
    cp -r lib/* ../install/lib
    cp -r lib/Release/* ../install/lib
fi

cd ${project_root_dir}/third_party
if [ ! -d "glfw" ]; then
    git clone https://github.com/glfw/glfw.git glfw
    cd ${project_root_dir}/third_party/glfw
    rm -rf build install
    mkdir -p build && mkdir -p install && cd build
    # cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/glfw/install -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_LIBDIR=lib -DGLFW_BUILD_EXAMPLES:BOOL=OFF -DGLFW_BUILD_WAYLAND:BOOL=OFF -DCMAKE_USER_MAKE_RULES_OVERRIDE=${project_root_dir}/ClangOverrides.txt ..
    cmake -DCMAKE_INSTALL_PREFIX=${project_root_dir}/third_party/glfw/install -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_LIBDIR=lib -DGLFW_BUILD_EXAMPLES:BOOL=OFF -DGLFW_BUILD_WAYLAND:BOOL=OFF ..
    cmake --build . -j 4 --config Release
    cmake --build . --config Release --target INSTALL
    cmake --install . --config Release
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
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ..

echo "convert this to an executable application -- let's go!!"
cmake --build . -j 4 --config Release 
cmake --install . --config Release
cd ${project_root_dir}
echo "declare success -- hooray!"

echo "running the executable with some default parameters"
echo "./build/main -c config.inp > results.txt 2>&1"
echo "  the 2>&1 redirects the stderr to a 1 so we don't see the gnuplot problems"
mkdir -p ${project_root_dir}/results
#./build/main -c config.inp > results/results.txt 2>&1
./build/main

# build the documentation
doxygen DOXYFILE

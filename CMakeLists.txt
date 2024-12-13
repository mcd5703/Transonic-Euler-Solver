cmake_minimum_required(VERSION 3.22) # Specify CMake version.
project(homework02
    DESCRIPTION "homework 02 solution set"
    VERSION 0.0.1) # A project name. This can be anything.
set(CMAKE_CXX_STANDARD 17) # Specify C++ version.
if(MSVC)
    message(NOTICE "Using the MSCV options")

    # Force to always compile with W4
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
        string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(NOTICE "Using the clang options")
    SET(GCC_COVERAGE_COMPILE_FLAGS "-fsanitize=address")
    SET(GCC_COVERAGE_LINK_FLAGS "-lgcov")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic")

    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    message(NOTICE "Using the gcc/c++ options")

    SET(GCC_COVERAGE_COMPILE_FLAGS "-fsanitize=address")
    SET(GCC_COVERAGE_LINK_FLAGS "-lgcov")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic")

    # SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
    # SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
    # SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
endif()
message(STATUS "CMakeLists.txt is at ${PROJECT_SOURCE_DIR}")

find_package(OpenGL REQUIRED)

find_package(GLEW)
if(MSVC)
    set(GLEW_FOUND FALSE)
endif()

if(GLEW_FOUND)
    include_directories(${GLEW_INCLUDE_DIRS})
    message(WARNING "Using the SYSTEM version of GLEW")
else()
    message(WARNING "Using the CUSTOM build version of GLEW")
    set(GLEW_DIR "${PROJECT_SOURCE_DIR}/third_party/glew/install")
    include_directories("${GLEW_DIR}/include")
    add_definitions(-DGLEW_STATIC) #to ensure static library building
    add_compile_definitions(GLEW_STATIC)
    if(MSVC)
        set(GLEW_STATIC_LIBRARIES "${GLEW_DIR}/lib/libglew32d.lib")
    else()
        set(GLEW_STATIC_LIBRARIES "${GLEW_DIR}/lib/libGLEW.a")
    endif()
endif()


list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/third_party/glfw/install")
find_package(glfw3 REQUIRED)

if(glfw3_FOUND)
    message(WARNING "Using the FOUND version of glfw3")
    set(GLFW_DIR "${PROJECT_SOURCE_DIR}/third_party/glfw/install")
    include_directories("${GLFW_DIR}/include")
    if(MSVC)
        set(GLFW_STATIC_LIBRARIES "${GLFW_DIR}/lib/glfw3.lib")
    else()
        set(GLFW_STATIC_LIBRARIES "${GLFW_DIR}/lib/libglfw3.a")
    endif()
else()
    message(FATAL_ERROR "glfw3 could not be found ... why??")
endif()

list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/third_party/matplotplusplus/install")
find_package(Matplot++ REQUIRED)
include(GNUInstallDirs)


include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

file(GLOB all_SRCS
    "${PROJECT_SOURCE_DIR}/include/*.h"
    "${PROJECT_SOURCE_DIR}/include/*.hpp"
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/*.c"
)
# Create a binary file. 
# The first argument is the name of the binary file.
# It can be anything. In this case, it is set to be the project name.
# The other arguments are source files.
# Header files are not needed in add_executable(), assuming they are properly included.
if(MSVC)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()
add_executable(main main/main.cpp ${all_SRCS})

target_link_libraries(main OpenGL::GL Matplot++::matplot)
target_link_libraries(main ${GLFW_STATIC_LIBRARIES})
if(MSVC)
    #    target_link_libraries(main ${GLFW_STATIC_LIBRARIES})
else()
    target_link_libraries(main m)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_link_libraries(main glfw)
    else()
        target_link_libraries(main ${GLFW_DIR}/lib/libglfw3.a)
    endif()
    target_link_libraries(main ${CMAKE_DL_LIBS})
    # target_link_libraries(main -fsanitize=undefined)
endif()



if(GLEW_FOUND)
    target_link_libraries(main GLEW::glew) # fix for MacOS / linux
else()
    target_link_libraries(main ${GLEW_STATIC_LIBRARIES})
endif()

add_custom_command(TARGET main
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:main> ${PROJECT_SOURCE_DIR}/build)

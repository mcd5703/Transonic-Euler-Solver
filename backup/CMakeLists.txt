cmake_minimum_required(VERSION 3.22)
project(homework02
    DESCRIPTION "homework 02 solution set"
    VERSION 0.0.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wall -Wpedantic -Wextra")
message(STATUS "CMakeLists.txt is at ${PROJECT_SOURCE_DIR}")
include(CTest)
enable_testing()

# Dependencies
find_package(OpenGL REQUIRED)

list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/third_party/glfw/install")
find_package(glfw3 REQUIRED)

list(APPEND CMAKE_PREFIX_PATH "${PROJECT_SOURCE_DIR}/third_party/glm/install")
find_package(glm REQUIRED)

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

add_executable(${CMAKE_PROJECT_NAME} ${all_SRCS})
target_link_libraries(${CMAKE_PROJECT_NAME} m)
target_link_libraries(${CMAKE_PROJECT_NAME} OpenGL::GL glfw Matplot++::matplot glm::glm)

# Copy timedata.csv to the build directory after building the executable
add_custom_command(TARGET ${CMAKE_PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/include/timedata.csv
    $<TARGET_FILE_DIR:${CMAKE_PROJECT_NAME}>)

# Enable packaging
set(CMAKE_INSTALL_PREFIX .)
install(FILES ${all_SRCS} DESTINATION .
    PERMISSIONS OWNER_READ OWNER_WRITE WORLD_READ WORLD_EXECUTE COMPONENT deb
)

set(CPACK_PROJECT_NAME ${CMAKE_PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${CMAKE_PROJECT_VERSION})
set(CPACK_GENERATOR "ZIP;TGZ")
set(CPACK_SOURCE_IGNORE_FILES
    /.git
    /dist
    /results
    /.*build.*
    /\\\\.DS_Store
)
set(CPACK_ARCHIVE_COMPONENT_INSTALL 1)
include(CPack)

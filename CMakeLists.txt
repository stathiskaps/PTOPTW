# CMakeList.txt : CMake project for AMTOPTW_algorithm, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")
set(CMAKE_BUILD_TYPE "DEBUG")

project ("AMTOPTW_algorithm")
set(INCLUDE_DIRECTORY "${PROJECT_SOURCE_DIR}/include")

set(HEADER_FILES "${INCLUDE_DIRECTORY}/Controller.h" "${INCLUDE_DIRECTORY}/Definitions.h"  "${INCLUDE_DIRECTORY}/Divider.h" "${INCLUDE_DIRECTORY}/ILS.h" "${INCLUDE_DIRECTORY}/List.h"  "${INCLUDE_DIRECTORY}/Processor.h" "${INCLUDE_DIRECTORY}/Solution.h" "${INCLUDE_DIRECTORY}/TaskManager.h" "${INCLUDE_DIRECTORY}/Unifier.h")

include_directories(${INCLUDE_DIRECTORY})

# Add source to this project's executable.
add_executable (AMTOPTW_algorithm "main.cpp" "Controller.cpp" "Definitions.cpp"  "Divider.cpp"  "Processor.cpp" "TaskManager.cpp" "Unifier.cpp" ${HEADER_FILES} "ILS.cpp"  "include/OP.h" "OP.cpp")

#target_include_directories(AMTOPTW_algorithm PRIVATE ${INCLUDE_DIRECTORY})

# TODO: Add tests and install targets if needed.

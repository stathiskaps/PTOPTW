# CMakeList.txt : CMake project for AMTOPTW_algorithm, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -pthread")
set(CMAKE_BUILD_TYPE "DEBUG")

project ("AMTOPTW_algorithm")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

set(HEADER_FILES "${INCLUDE_DIR}/Definitions.h" "${INCLUDE_DIR}/ILS.h" "${INCLUDE_DIR}/List.h" "${INCLUDE_DIR}/Solution.h" "${INCLUDE_DIR}/OP.h" "${INCLUDE_DIR}/Custom.h" "${INCLUDE_DIR}/Graphics.h")
set(SOURCE "main.cpp" "Definitions.cpp" "ILS.cpp" "Solution.cpp" "OP.cpp" "Graphics.cpp")

include_directories(${INCLUDE_DIR})

# Add source to this project's executable.
add_executable (AMTOPTW ${HEADER_FILES} ${SOURCE})
target_link_libraries(AMTOPTW GL GLU glut)

#target_include_directories(AMTOPTW_algorithm PRIVATE ${INCLUDE_DIR})

# TODO: Add tests and install targets if needed.

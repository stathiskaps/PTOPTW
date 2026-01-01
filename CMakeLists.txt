cmake_minimum_required(VERSION 3.16)
project(PTOPTW LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(ptoptw
  src/main.cpp
  src/ILS.cpp
  src/OP.cpp
  src/Solution.cpp
  src/Definitions.cpp
)

target_include_directories(ptoptw PRIVATE
  include
  third_party
)

# optional: warnings
if (MSVC)
  target_compile_options(ptoptw PRIVATE /W4)
else()
  target_compile_options(ptoptw PRIVATE -Wall -Wextra -Wpedantic)
endif()

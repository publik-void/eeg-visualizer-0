cmake_minimum_required(VERSION 3.13)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# CMake FAQ, maybe this will come in handy:
# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (${CMAKE_GENERATOR} STREQUAL "Ninja Multi-Config")
  set(CMAKE_DEFAULT_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif()

project(eeg-visualizer-0)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0" " -ggdb3")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang") # "Clang" or "AppleClang"
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0" " -glldb")
endif()

add_compile_options(-W -Wall -Wextra)

# For "Ninja" or "Ninja Multi-Config", colorized output has to be forced
if (${CMAKE_GENERATOR} MATCHES "Ninja")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    string(APPEND CMAKE_CXX_FLAGS " -fdiagnostics-color=always")
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang") # "Clang" or "AppleClang"
    string(APPEND CMAKE_CXX_FLAGS " -fcolor-diagnostics")
  endif()
endif()

add_library(lsl-imported INTERFACE IMPORTED)
target_link_libraries(lsl-imported INTERFACE lsl)
target_include_directories(lsl-imported SYSTEM INTERFACE /usr/local/include)
target_link_directories(lsl-imported INTERFACE /usr/local/lib)

add_library(lsl-customized STATIC src/lsl-customized.cpp)
target_link_libraries(lsl-customized PUBLIC lsl-imported)

add_library(signal-handling STATIC src/signal-handling.cpp)

add_library(args STATIC src/args.cpp)

add_library(dsp STATIC src/dsp.cpp)

add_library(server-dummy STATIC src/server-dummy.cpp)
target_link_libraries(server-dummy PUBLIC lsl-customized)

add_library(client-visual STATIC src/client-visual.cpp)
target_link_libraries(client-visual PUBLIC lsl-customized sfml-system
  sfml-window sfml-graphics fftw3f)

add_executable(main src/main.cpp)
target_link_libraries(main PUBLIC signal-handling args dsp server-dummy
  client-visual)


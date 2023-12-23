cmake_minimum_required(VERSION 3.19)
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

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)

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

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
find_package(LSL REQUIRED)
find_package(SFML  COMPONENTS window graphics REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3f IMPORTED_TARGET)

add_library(lsl-imported INTERFACE IMPORTED)
target_link_libraries(lsl-imported INTERFACE LSL::lsl)

add_library(lsl-customized STATIC src/lsl-customized.cpp)
target_link_libraries(lsl-customized PUBLIC lsl-imported)

add_library(signal-handling STATIC src/signal-handling.cpp)

add_library(args STATIC src/args.cpp)

add_library(number-printing STATIC src/number-printing.cpp)

add_library(dsp INTERFACE src/dsp.hpp)

add_library(server-dummy STATIC src/server-dummy.cpp)
target_link_libraries(server-dummy PUBLIC lsl-customized)

add_library(plotting STATIC src/plotting.cpp)
target_link_libraries(plotting PUBLIC sfml-graphics)

add_library(client-visualizer STATIC src/client-visualizer.cpp)
target_link_libraries(client-visualizer PUBLIC lsl-customized plotting
  number-printing sfml-window sfml-graphics PkgConfig::FFTW)

add_executable(eeg-visualizer src/main.cpp)
target_link_libraries(eeg-visualizer PUBLIC signal-handling args dsp
  server-dummy client-visualizer Threads::Threads)


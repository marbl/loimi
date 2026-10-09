# SPDX-License-Identifier: MIT-0
cmake_minimum_required(VERSION 3.18)
project(hpc_consensus VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(VERKKO_ROOT "" CACHE PATH "Verkko checkout with build/lib/libverkko.a")
# cnspack serialization is native-ABI dependent: compile against the same
# headers/build used for the packages rather than vendoring a guessed reader.
set(VUTIL "${VERKKO_ROOT}/src/canu/src/utility/src")
if(NOT EXISTS "${VERKKO_ROOT}/build/lib/libverkko.a")
  message(FATAL_ERROR "Set VERKKO_ROOT to a built Verkko checkout")
endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(SYSTEM_DEPS REQUIRED IMPORTED_TARGET openssl libcurl liblzma)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
# The imported OpenMP target supplies compiler flags as well as the runtime;
# locating a libomp/libgomp filename alone is not portable across Linux/macOS.
add_library(consensus src/core.cpp)
target_include_directories(consensus PUBLIC src)
# Quote-only lookup avoids utility/math.H shadowing the system math.h on macOS.
target_compile_options(consensus PUBLIC "-iquote${VUTIL}" "-iquote${VERKKO_ROOT}/src/canu/src/stores")
target_link_libraries(consensus PUBLIC "${VERKKO_ROOT}/build/lib/libverkko.a" ${SYSTEM_DEPS_LINK_LIBRARIES} ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES} Threads::Threads OpenMP::OpenMP_CXX)
target_compile_options(consensus PRIVATE -Wall -Wextra -Wpedantic)
add_executable(hpc-consensus src/main.cpp)
target_link_libraries(hpc-consensus PRIVATE consensus)
install(TARGETS hpc-consensus RUNTIME DESTINATION bin)
enable_testing()
# Core tests need no input datasets. Integration generates small native packages
# with Verkko tools, and exercises CRAM only when samtools is available on PATH.
add_executable(core-tests tests/core.cpp)
target_link_libraries(core-tests PRIVATE consensus)
add_executable(extract-regions tests/extract_regions.cpp)
target_link_libraries(extract-regions PRIVATE consensus)
add_executable(apply-run-model tests/apply_run_model.cpp)
target_link_libraries(apply-run-model PRIVATE consensus)
add_test(NAME core COMMAND core-tests)
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(NAME native-package COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration.py" "$<TARGET_FILE:hpc-consensus>" "${VERKKO_ROOT}")
add_test(NAME calibration COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/calibration.py")

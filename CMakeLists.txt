cmake_minimum_required(VERSION 3.13)
project(astap_cpp_platesolver CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif ()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Optional SYCL offload. Works with either toolchain:
#   AdaptiveCpp:  cmake -DASTAP_SYCL=ON -DAdaptiveCpp_DIR=...
#   oneAPI/DPC++: cmake -DASTAP_SYCL=ON -DCMAKE_CXX_COMPILER=icpx
option(ASTAP_SYCL "Build the SYCL offload backend" OFF)
set(ASTAP_SYCL_TARGETS "" CACHE STRING
        "SYCL targets, e.g. nvptx64-nvidia-cuda or spir64. Empty lets the toolchain decide.")

set(ASTAP_SYCL_SOURCES "")
if (ASTAP_SYCL)
    find_package(AdaptiveCpp CONFIG QUIET)
    if (AdaptiveCpp_FOUND)
        set(ASTAP_SYCL_FLAVOUR "AdaptiveCpp")
    else ()
        include(CheckCXXSourceCompiles)
        set(CMAKE_REQUIRED_FLAGS "-fsycl")
        set(CMAKE_REQUIRED_LINK_OPTIONS "-fsycl")
        check_cxx_source_compiles("#include <sycl/sycl.hpp>\nint main(){ sycl::queue q; return 0; }"
                ASTAP_HAVE_FSYCL)
        unset(CMAKE_REQUIRED_FLAGS)
        unset(CMAKE_REQUIRED_LINK_OPTIONS)
        if (NOT ASTAP_HAVE_FSYCL)
            message(FATAL_ERROR
                    "ASTAP_SYCL=ON but no SYCL toolchain found. Use -DCMAKE_CXX_COMPILER=icpx "
                    "for oneAPI, or point AdaptiveCpp_DIR at an AdaptiveCpp install.")
        endif ()
        set(ASTAP_SYCL_FLAVOUR "DPC++")
        # icpx defaults to -fp-model=fast, which reassociates and uses reduced
        # precision math. The whole port depends on IEEE-754 semantics, so pin
        # them for host and device alike, or results stop matching the reference.
        add_compile_options(-fp-model=precise -ffp-contract=off)
    endif ()
    set(ASTAP_SYCL_SOURCES src/quad_batch_sycl.cpp)
    message(STATUS "SYCL offload enabled (${ASTAP_SYCL_FLAVOUR})")
endif ()

add_library(astap_solver
        src/astro_math.cpp
        src/parallel.cpp
  src/quad_batch.cpp
        src/calc_trans_cubic.cpp
        src/fits.cpp
        src/matching.cpp
        src/quads.cpp
        src/star_database.cpp
        src/star_detection.cpp
        src/solver.cpp
        ${ASTAP_SYCL_SOURCES}
)
target_include_directories(astap_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(astap_solver PUBLIC Threads::Threads)

if (ASTAP_SYCL)
    target_compile_definitions(astap_solver PUBLIC ASTAP_WITH_SYCL)
    if (AdaptiveCpp_FOUND)
        # Without this the device backend fuses a*b+c into an FMA, which changes
        # the last bit and breaks the bit-identical guarantee. NVPTX contracts by
        # default, so this is not optional.
        target_compile_options(astap_solver PRIVATE -ffp-contract=off)
        add_sycl_to_target(TARGET astap_solver SOURCES ${ASTAP_SYCL_SOURCES})
    else ()
        # -ffp-contract=off keeps the device from fusing a*b+c into an FMA, which
        # would change the last bit and break the bit-identical guarantee.
        set(_sycl_flags -fsycl -ffp-contract=off)
        if (ASTAP_SYCL_TARGETS)
            list(APPEND _sycl_flags -fsycl-targets=${ASTAP_SYCL_TARGETS})
        endif ()
        target_compile_options(astap_solver PRIVATE ${_sycl_flags})
        target_link_options(astap_solver PUBLIC ${_sycl_flags})
    endif ()
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(astap_solver PRIVATE -Wall -Wextra)
    # The port relies on IEEE-754 semantics and on the default rounding mode being
    # round-half-to-even, the same as FPC's Round(). Do not enable -ffast-math.
endif ()

add_executable(astap_solve src/main.cpp)
target_link_libraries(astap_solve PRIVATE astap_solver)

enable_testing()
add_executable(solver_tests tests/solver_tests.cpp)
target_link_libraries(solver_tests PRIVATE astap_solver)
add_test(NAME solver_tests COMMAND solver_tests)

add_executable(quad_batch_tests tests/quad_batch_tests.cpp)
target_link_libraries(quad_batch_tests PRIVATE astap_solver)
add_test(NAME quad_batch_tests COMMAND quad_batch_tests)

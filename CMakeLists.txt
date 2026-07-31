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

add_library(astap_solver
        src/astro_math.cpp
        src/parallel.cpp
        src/calc_trans_cubic.cpp
        src/fits.cpp
        src/matching.cpp
        src/quads.cpp
        src/star_database.cpp
        src/star_detection.cpp
        src/solver.cpp
)
target_include_directories(astap_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(astap_solver PUBLIC Threads::Threads)

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

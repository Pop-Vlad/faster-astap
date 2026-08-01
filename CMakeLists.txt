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
        src/quad_batch.cpp
        src/quad_index.cpp
        src/index_solver.cpp
        src/calc_trans_cubic.cpp
        src/matching.cpp
        src/quads.cpp
        src/star_database.cpp
        src/star_detection.cpp
        src/solver.cpp
)
target_include_directories(astap_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(astap_solver PUBLIC Threads::Threads)

# --- image module -------------------------------------------------------------
# Everything that turns a file on disk into an ImageArray, the equivalent of
# load_image() in the original. FITS (including Rice compressed .fz), Netpbm and
# BMP need nothing beyond the compiler. The other formats astap_cli reads are
# decoded by the usual system libraries when they are installed; without them
# the build reports those formats as unsupported and everything else still
# works. Set any of these to OFF to leave a format out on purpose.
option(ASTAP_ZLIB "Read GZIP compressed FITS tiles when zlib is available" ON)
option(ASTAP_PNG "Read PNG images when libpng is available" ON)
option(ASTAP_JPEG "Read JPEG images when libjpeg is available" ON)
option(ASTAP_TIFF "Read TIFF images when libtiff is available" ON)
option(ASTAP_LIBRAW "Read raw camera images when LibRaw is available" ON)

add_library(astap_image
        src/image/image_io.cpp
        src/image/fits.cpp
        src/image/rice.cpp
        src/image/image_pnm.cpp
        src/image/image_bmp.cpp
        src/image/image_png.cpp
        src/image/image_jpeg.cpp
        src/image/image_tiff.cpp
        src/image/image_raw.cpp
)
target_include_directories(astap_image PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# For the shared types, astro_math and the thread pool. Nothing in the solver
# depends on this module, so the two form no cycle.
target_link_libraries(astap_image PUBLIC astap_solver)

set(ASTAP_IMAGE_FORMATS "FITS;FITS.fz (Rice);Netpbm;BMP")

if (ASTAP_ZLIB)
    # Needed for the tiles fpack stores as GZIP because Rice did not shrink
    # them, and for a wholly GZIP compressed .fz.
    find_package(ZLIB QUIET)
    if (ZLIB_FOUND)
        target_compile_definitions(astap_image PRIVATE ASTAP_WITH_ZLIB)
        target_link_libraries(astap_image PRIVATE ZLIB::ZLIB)
        list(APPEND ASTAP_IMAGE_FORMATS "FITS.fz (GZIP)")
    endif ()
endif ()

if (ASTAP_PNG)
    find_package(PNG QUIET)
    if (PNG_FOUND)
        target_compile_definitions(astap_image PRIVATE ASTAP_WITH_PNG)
        target_link_libraries(astap_image PRIVATE PNG::PNG)
        list(APPEND ASTAP_IMAGE_FORMATS "PNG")
    endif ()
endif ()

if (ASTAP_JPEG)
    find_package(JPEG QUIET)
    if (JPEG_FOUND)
        target_compile_definitions(astap_image PRIVATE ASTAP_WITH_JPEG)
        target_link_libraries(astap_image PRIVATE JPEG::JPEG)
        list(APPEND ASTAP_IMAGE_FORMATS "JPEG")
    endif ()
endif ()

if (ASTAP_TIFF)
    find_package(TIFF QUIET)
    if (TIFF_FOUND)
        target_compile_definitions(astap_image PRIVATE ASTAP_WITH_TIFF)
        target_link_libraries(astap_image PRIVATE TIFF::TIFF)
        list(APPEND ASTAP_IMAGE_FORMATS "TIFF")
    endif ()
endif ()

if (ASTAP_LIBRAW)
    # LibRaw ships a pkg-config file rather than a CMake config.
    find_package(PkgConfig QUIET)
    if (PkgConfig_FOUND)
        pkg_check_modules(LIBRAW QUIET IMPORTED_TARGET libraw)
    endif ()
    if (LIBRAW_FOUND)
        target_compile_definitions(astap_image PRIVATE ASTAP_WITH_LIBRAW)
        target_link_libraries(astap_image PRIVATE PkgConfig::LIBRAW)
        list(APPEND ASTAP_IMAGE_FORMATS "raw camera files")
    endif ()
endif ()

# supported_image_extensions() reports the same set at run time.
string(REPLACE ";" ", " ASTAP_IMAGE_FORMATS_TEXT "${ASTAP_IMAGE_FORMATS}")
message(STATUS "Image formats: ${ASTAP_IMAGE_FORMATS_TEXT}")

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(astap_solver PRIVATE -Wall -Wextra)
    target_compile_options(astap_image PRIVATE -Wall -Wextra)
    # The port relies on IEEE-754 semantics and on the default rounding mode being
    # round-half-to-even, the same as FPC's Round(). Do not enable -ffast-math.
endif ()

add_executable(astap_solve src/main.cpp)
target_link_libraries(astap_solve PRIVATE astap_image)

# The index method, same options and same .ini/.wcs output as astap_solve.
add_executable(astap_index_solve src/index_main.cpp)
target_link_libraries(astap_index_solve PRIVATE astap_image)

enable_testing()
add_executable(solver_tests tests/solver_tests.cpp)
target_link_libraries(solver_tests PRIVATE astap_solver)
add_test(NAME solver_tests COMMAND solver_tests)

add_executable(image_io_tests tests/image_io_tests.cpp)
target_link_libraries(image_io_tests PRIVATE astap_image)
target_compile_definitions(image_io_tests PRIVATE
        ASTAP_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data"
        ASTAP_TEST_TMP_DIR="${CMAKE_CURRENT_BINARY_DIR}")
add_test(NAME image_io_tests COMMAND image_io_tests)

add_executable(quad_index_bench tests/quad_index_bench.cpp)
target_link_libraries(quad_index_bench PRIVATE astap_image)

add_executable(corpus_harness tools/corpus_harness.cpp)
target_link_libraries(corpus_harness PRIVATE astap_image)

add_executable(quad_batch_tests tests/quad_batch_tests.cpp)
target_link_libraries(quad_batch_tests PRIVATE astap_solver)
add_test(NAME quad_batch_tests COMMAND quad_batch_tests)

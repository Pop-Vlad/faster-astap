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

# A MinGW build otherwise needs libstdc++, libgcc and libwinpthread next to the
# .exe or on PATH, and an .exe that will not start unless MSYS2 is on PATH is of
# no use to the imaging application that launches it. So this is on by default
# on Windows, where those libraries are not something the system provides
# anyway; off elsewhere, where a package shares the system copies as it should.
if (WIN32)
    set(ASTAP_STATIC_RUNTIME_DEFAULT ON)
else ()
    set(ASTAP_STATIC_RUNTIME_DEFAULT OFF)
endif ()
option(ASTAP_STATIC_RUNTIME "Link the runtime and image libraries into the executables (Windows)"
        ${ASTAP_STATIC_RUNTIME_DEFAULT})

# The image libraries are linked from their static archives in the same build,
# so the two executables are the whole distribution. -static is what makes the
# linker prefer libfoo.a over the libfoo.dll.a import library for each -lfoo,
# which is why the link line below has to come from pkg-config as -l flags
# rather than as the absolute paths the CMake find modules return.
set(ASTAP_STATIC_IMAGE_LIBS OFF)
if (ASTAP_STATIC_RUNTIME AND WIN32)
    if (MINGW)
        add_link_options(-static-libgcc -static-libstdc++ -static)
        set(ASTAP_STATIC_IMAGE_LIBS ON)
    elseif (MSVC AND NOT CMAKE_VERSION VERSION_LESS 3.15)
        cmake_policy(SET CMP0091 NEW)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif ()
endif ()

add_library(astap_solver
        src/astro_math.cpp
        src/mapped_file.cpp
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
        # The per image pipeline of both solvers. It needs no image decoder —
        # the caller supplies pixels — so it belongs to this module, and the
        # file front ends that read and write around it live in astap_image.
        src/solve_service.cpp
        src/spiral_service.cpp
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
        # The file front ends: read an image, hand it to the pipeline in
        # astap_solver, write the .ini and .wcs. Here because this is the module
        # that has the decoder, and the pipeline itself deliberately does not.
        src/solve_service_files.cpp
        src/spiral_service_files.cpp
)
target_include_directories(astap_image PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# For the shared types, astro_math and the thread pool. Nothing in the solver
# depends on this module, so the two form no cycle.
target_link_libraries(astap_image PUBLIC astap_solver)

set(ASTAP_IMAGE_FORMATS "FITS;FITS.fz (Rice);Netpbm;BMP")

find_package(PkgConfig QUIET)

# Adds one optional image library, and records the formats it brings.
#
# A static build asks pkg-config for the link line, because only that lists the
# transitive dependencies: libtiff on its own needs jbig, webp, sharpyuv, zstd,
# lzma, Lerc and deflate, and LibRaw needs lcms2 and OpenMP. That list is read
# from the libraries' own metadata rather than written out here, so it keeps up
# with them. Otherwise the ordinary CMake find module is used and the library
# stays a DLL.
#
# A macro, not a function, so that appending to ASTAP_IMAGE_FORMATS is visible
# to the caller.
macro(astap_image_library pkg_module find_name imported_target define)
    set(_astap_done FALSE)
    if (ASTAP_STATIC_IMAGE_LIBS AND PkgConfig_FOUND)
        pkg_check_modules(_astap_pc QUIET ${pkg_module})
        if (_astap_pc_FOUND)
            target_include_directories(astap_image PRIVATE ${_astap_pc_INCLUDE_DIRS})
            # The static cflags matter as much as the libraries: LibRaw declares
            # its symbols __declspec(dllimport) unless -DLIBRAW_NODLL says the
            # link is static, and its .pc file is what knows that.
            target_compile_options(astap_image PRIVATE ${_astap_pc_STATIC_CFLAGS_OTHER})
            target_link_libraries(astap_image PRIVATE ${_astap_pc_STATIC_LDFLAGS})
            set(_astap_done TRUE)
        endif ()
    endif ()
    if (NOT _astap_done AND NOT "${find_name}" STREQUAL "")
        find_package(${find_name} QUIET)
        if (${find_name}_FOUND)
            target_link_libraries(astap_image PRIVATE ${imported_target})
            set(_astap_done TRUE)
        endif ()
    endif ()
    if (_astap_done)
        target_compile_definitions(astap_image PRIVATE ${define})
        list(APPEND ASTAP_IMAGE_FORMATS ${ARGN})
    endif ()
endmacro()

if (ASTAP_ZLIB)
    # Needed for the tiles fpack stores as GZIP because Rice did not shrink
    # them, and for a wholly GZIP compressed .fz.
    astap_image_library(zlib ZLIB ZLIB::ZLIB ASTAP_WITH_ZLIB "FITS.fz (GZIP)")
endif ()

if (ASTAP_PNG)
    astap_image_library(libpng16 PNG PNG::PNG ASTAP_WITH_PNG "PNG")
endif ()

if (ASTAP_JPEG)
    astap_image_library(libjpeg JPEG JPEG::JPEG ASTAP_WITH_JPEG "JPEG")
endif ()

if (ASTAP_TIFF)
    astap_image_library(libtiff-4 TIFF TIFF::TIFF ASTAP_WITH_TIFF "TIFF")
endif ()

if (ASTAP_LIBRAW)
    # LibRaw ships a pkg-config file rather than a CMake config, so there is no
    # find module to fall back on. Windows rarely has pkg-config outside MSYS2;
    # vcpkg installs a CMake config under the unofficial- name instead.
    set(_astap_done FALSE)
    if (PkgConfig_FOUND)
        if (ASTAP_STATIC_IMAGE_LIBS)
            pkg_check_modules(LIBRAW QUIET libraw)
            if (LIBRAW_FOUND)
                target_include_directories(astap_image PRIVATE ${LIBRAW_INCLUDE_DIRS})
                # -DLIBRAW_NODLL, without which the headers ask for dllimport
                # symbols that a static archive does not provide.
                target_compile_options(astap_image PRIVATE ${LIBRAW_STATIC_CFLAGS_OTHER})
                target_link_libraries(astap_image PRIVATE ${LIBRAW_STATIC_LDFLAGS})
                # htonl and friends live in ws2_32 on Windows. LibRaw's .pc file
                # leaves them out because linking the DLL resolves them on its
                # own; a static archive does not.
                if (WIN32)
                    target_link_libraries(astap_image PRIVATE ws2_32)
                endif ()
                set(_astap_done TRUE)
            endif ()
        else ()
            pkg_check_modules(LIBRAW QUIET IMPORTED_TARGET libraw)
            if (LIBRAW_FOUND)
                target_link_libraries(astap_image PRIVATE PkgConfig::LIBRAW)
                set(_astap_done TRUE)
            endif ()
        endif ()
    endif ()
    if (NOT _astap_done)
        find_package(unofficial-libraw CONFIG QUIET)
        if (unofficial-libraw_FOUND)
            target_link_libraries(astap_image PRIVATE unofficial::libraw::raw)
            set(_astap_done TRUE)
        endif ()
    endif ()
    if (_astap_done)
        target_compile_definitions(astap_image PRIVATE ASTAP_WITH_LIBRAW)
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
elseif (MSVC)
    # /utf-8 because the sources are UTF-8 and MSVC otherwise reads them in the
    # system codepage. /fp:precise is the same IEEE-754 requirement as above:
    # do not switch this to /fp:fast.
    target_compile_options(astap_solver PRIVATE /W3 /utf-8 /permissive- /fp:precise)
    target_compile_options(astap_image PRIVATE /W3 /utf-8 /permissive- /fp:precise)
    # getenv, snprintf and fopen are used deliberately and portably.
    target_compile_definitions(astap_solver PUBLIC _CRT_SECURE_NO_WARNINGS)
endif ()

if (WIN32)
    # Windows headers otherwise define min/max as macros, which breaks
    # std::min/std::max at every call site that includes them.
    target_compile_definitions(astap_solver PUBLIC NOMINMAX WIN32_LEAN_AND_MEAN)
endif ()

add_executable(astap_solve src/main.cpp)
target_link_libraries(astap_solve PRIVATE astap_image)

# The index method, same options and same .ini/.wcs output as astap_solve.
add_executable(astap_index_solve src/index_main.cpp)
target_link_libraries(astap_index_solve PRIVATE astap_image)

# --- install, and the Windows bundle --------------------------------------
# `cmake --install build --prefix dist` gives a directory holding the two
# solvers. On Windows they go in the root, so the directory can be zipped and
# the DLLs collected below sit next to the .exe that loads them; elsewhere they
# go in bin/, as a prefix is expected to be laid out.
if (WIN32)
    set(ASTAP_RUNTIME_DEST .)
else ()
    set(ASTAP_RUNTIME_DEST bin)
endif ()
install(TARGETS astap_solve astap_index_solve RUNTIME DESTINATION ${ASTAP_RUNTIME_DEST})

# --- the N.I.N.A. integration ------------------------------------------------
# A module of its own, in nina-plugin/: a resident server that keeps the index in
# memory between images, and the C# plugin that drives it. It consumes the
# libraries above and adds nothing to them, so turning it off changes nothing
# about the two solvers.
option(ASTAP_NINA_INTEGRATION "Build the N.I.N.A. integration in nina-plugin/" ON)
if (ASTAP_NINA_INTEGRATION AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/nina-plugin/CMakeLists.txt)
    add_subdirectory(nina-plugin)
endif ()

# With the image libraries switched on, a Windows build links libpng, libjpeg,
# libtiff, LibRaw and zlib as DLLs, and each of those pulls in more of its own
# (libtiff alone wants jbig, webp, zstd, lzma, Lerc and deflate). Linking them
# statically instead would mean hard coding that list, which changes with every
# library release, so the install step resolves the closure and copies it next
# to the .exe. The result runs on a machine with no MSYS2 on it.
#
# Only the third-party DLLs are collected: everything shipped with Windows is
# excluded, and so is anything the resolver could not find, which is reported
# rather than silently dropped.
if (WIN32 AND NOT CMAKE_VERSION VERSION_LESS 3.16)
    # The DLLs sit beside the compiler, in the environment's bin directory.
    get_filename_component(ASTAP_TOOLCHAIN_BIN "${CMAKE_CXX_COMPILER}" DIRECTORY)
    install(CODE "set(ASTAP_DLL_DIRS [[${ASTAP_TOOLCHAIN_BIN}]])")
    install(CODE [[
        # The resolver hands back paths with mixed separators
        # ("C:\WINDOWS\system32/kernel32.dll"). Normalising them is what the
        # newer policies do, and without these CMake warns once per DLL.
        if (POLICY CMP0207)
            cmake_policy(SET CMP0207 NEW)
        endif ()
        if (POLICY CMP0177)
            cmake_policy(SET CMP0177 NEW)
        endif ()
        file(GET_RUNTIME_DEPENDENCIES
                EXECUTABLES "${CMAKE_INSTALL_PREFIX}/astap_solve.exe"
                            "${CMAKE_INSTALL_PREFIX}/astap_index_solve.exe"
                RESOLVED_DEPENDENCIES_VAR astap_resolved
                UNRESOLVED_DEPENDENCIES_VAR astap_unresolved
                DIRECTORIES ${ASTAP_DLL_DIRS})
        # Bundle exactly what came out of the toolchain directory. Filtering on
        # where a DLL lives rather than on its name is what keeps kernel32 and
        # friends out: Windows reports its own path in a mix of cases
        # (C:/WINDOWS/System32), which no single regex matches reliably.
        # Compared as a plain prefix rather than as a regex, so a toolchain
        # installed under a path holding '.' or '+' still matches.
        string(TOLOWER "${ASTAP_DLL_DIRS}/" astap_bin_lower)
        foreach (astap_dll ${astap_resolved})
            string(TOLOWER "${astap_dll}" astap_dll_lower)
            string(FIND "${astap_dll_lower}" "${astap_bin_lower}" astap_at)
            if (astap_at EQUAL 0)
                file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}" TYPE SHARED_LIBRARY
                        FOLLOW_SYMLINK_CHAIN FILES "${astap_dll}")
            endif ()
        endforeach ()
        # api-ms-win-* are the Universal CRT, part of Windows 10 and later, so
        # they are expected here and are not a problem.
        foreach (astap_missing ${astap_unresolved})
            if (NOT astap_missing MATCHES "^api-ms-win")
                message(WARNING "not bundled, and not part of Windows: ${astap_missing}")
            endif ()
        endforeach ()
    ]])
endif ()

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

# Solving an array has to give the same answer as solving the file it came from.
# Reports itself skipped, and passes, without a star database or a corpus image.
add_executable(service_api_tests tests/service_api_tests.cpp)
target_link_libraries(service_api_tests PRIVATE astap_image)
target_compile_definitions(service_api_tests PRIVATE
        ASTAP_TEST_CORPUS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/corpus"
        ASTAP_TEST_TMP_DIR="${CMAKE_CURRENT_BINARY_DIR}")
add_test(NAME service_api_tests COMMAND service_api_tests)

# Reports itself skipped, and passes, when the machine has no star database to
# build an index from.
add_executable(index_file_tests tests/index_file_tests.cpp)
target_link_libraries(index_file_tests PRIVATE astap_solver)
target_compile_definitions(index_file_tests PRIVATE
        ASTAP_TEST_TMP_DIR="${CMAKE_CURRENT_BINARY_DIR}")
add_test(NAME index_file_tests COMMAND index_file_tests)

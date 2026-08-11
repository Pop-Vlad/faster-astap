# The N.I.N.A. integration, built as a module of its own.
#
# Nothing in ../src or ../include refers to anything here. This directory
# consumes the solver libraries the parent builds and adds one executable to
# them; removing the directory, or turning ASTAP_NINA_INTEGRATION off, leaves
# astap_solve and astap_index_solve exactly as they were.
#
# The C# plugin that goes with this executable is in plugin/, and is built by
# the .NET SDK rather than by CMake - see README.md.

# The version resource is not decoration: N.I.N.A. reads the executable's
# FileVersion before it will pass its default options. See solver/version.rc.
if (WIN32)
    enable_language(RC)
    add_executable(astap_nina_solve
            solver/nina_solve_main.cpp
            solver/version.rc
    )
else ()
    add_executable(astap_nina_solve
            solver/nina_solve_main.cpp
    )
endif ()

target_link_libraries(astap_nina_solve PRIVATE astap_image)

install(TARGETS astap_nina_solve RUNTIME DESTINATION ${ASTAP_RUNTIME_DEST})

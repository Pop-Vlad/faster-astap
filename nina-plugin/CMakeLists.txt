# The N.I.N.A. integration, built as a module of its own.
#
# Nothing in ../src or ../include refers to anything here. This directory
# consumes the solver libraries the parent builds and adds one executable to
# them; removing the directory, or turning ASTAP_NINA_INTEGRATION off, leaves
# astap_solve and astap_index_solve exactly as they were.
#
# The C# plugin that goes with this executable is in src/ as well, and is built
# by the .NET SDK rather than by CMake — see README.md.

# The version resource is not decoration: N.I.N.A. reads the executable's
# FileVersion before it will pass its default options. See server/version.rc.
if (WIN32)
    enable_language(RC)
    add_executable(astap_index_server
            server/index_server_main.cpp
            server/ipc.cpp
            server/version.rc
    )
else ()
    add_executable(astap_index_server
            server/index_server_main.cpp
            server/ipc.cpp
    )
endif ()

target_include_directories(astap_index_server PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/server)
target_link_libraries(astap_index_server PRIVATE astap_image)

install(TARGETS astap_index_server RUNTIME DESTINATION ${ASTAP_RUNTIME_DEST})

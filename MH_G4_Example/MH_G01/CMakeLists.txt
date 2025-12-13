cmake_minimum_required(VERSION 3.16)
project(MH_G01 LANGUAGES CXX)

# --------------------------------------------------
# Find Geant4 (GDML only)
# --------------------------------------------------
find_package(Geant4 REQUIRED gdml)
include(${Geant4_USE_FILE})

# --------------------------------------------------
# Include directories
# --------------------------------------------------
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIRS}
)

# --------------------------------------------------
# Common sources (everything in src/)
# --------------------------------------------------
file(GLOB COMMON_SOURCES
    "${PROJECT_SOURCE_DIR}/src/*.cc"
)

# --------------------------------------------------
# Shared library
# --------------------------------------------------
add_library(MH_G01_lib ${COMMON_SOURCES})
target_link_libraries(MH_G01_lib ${Geant4_LIBRARIES})

# --------------------------------------------------
# Executables (these contain main())
# --------------------------------------------------
add_executable(load_gdml load_gdml.cc)
add_executable(load_gdml_color load_gdml_color.cc)

target_link_libraries(load_gdml
    MH_G01_lib
    ${Geant4_LIBRARIES}
)

target_link_libraries(load_gdml_color
    MH_G01_lib
    ${Geant4_LIBRARIES}
)

install(TARGETS load_gdml load_gdml_color
        RUNTIME DESTINATION bin)


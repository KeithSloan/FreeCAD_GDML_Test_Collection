# Set the minimum CMake version
cmake_minimum_required(VERSION 3.16...3.27)

# Project name
project(MyGeant4App)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_VIS "Build example with Geant4 Vis drivers" ON)
if(WITH_GEANT4_VIS)
  find_package(Geant4 REQUIRED gdml vis_all)
else()
  find_package(Geant4 REQUIRED gdml)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
# Find the Geant4 package using CMake


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)



# Set C++ standard to C++11 or C++14 (depending on your Geant4 version)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the sources (replace main.cpp with your actual source file names)
#set(SOURCES src/parseCheck.cpp src/YourOtherFiles.cpp)
set(SOURCES parseCheck.cpp)

# Define the executable
add_executable(MyGeant4App ${SOURCES})

# Link the Geant4 libraries to your project
target_link_libraries(MyGeant4App ${Geant4_LIBRARIES})

# Include the Geant4 include directories
target_include_directories(MyGeant4App PUBLIC ${Geant4_INCLUDE_DIRS})

# Optionally enable all warnings
target_compile_options(MyGeant4App PRIVATE -Wall -Wextra -pedantic)

# Include Geant4's mandatory CMake macros and variables
include(${Geant4_USE_FILE})

# Optionally install the application
install(TARGETS MyGeant4App DESTINATION bin)


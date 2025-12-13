# Install script for directory: /Users/ksloan/geant4-v11.3.2

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/ksloan/geant4-v11.3.2-full-install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/source/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Development" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/Geant4/geant4make" TYPE FILE MESSAGE_LAZY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/InstallTreeFiles/geant4make.sh")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Development" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/Geant4/geant4make" TYPE FILE MESSAGE_LAZY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/InstallTreeFiles/geant4make.csh")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Development" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/Geant4/geant4make" TYPE DIRECTORY MESSAGE_LAZY FILES "/Users/ksloan/geant4-v11.3.2/config" FILES_MATCHING REGEX "/[^/]*\\.gmk$" REGEX "/cvs$" EXCLUDE REGEX "/\\.svn$" EXCLUDE REGEX "/scripts\\/$" EXCLUDE)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory "$ENV{DESTDIR}/Users/ksloan/geant4-v11.3.2-full-install/lib/Geant4-11.3.2")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink .. Darwin-clang WORKING_DIRECTORY "$ENV{DESTDIR}/Users/ksloan/geant4-v11.3.2-full-install/lib/Geant4-11.3.2")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Runtime" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE MESSAGE_LAZY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/InstallTreeFiles/geant4.sh")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Runtime" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE MESSAGE_LAZY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/InstallTreeFiles/geant4.csh")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Development" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE FILE MESSAGE_LAZY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/InstallTreeFiles/geant4-config")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Development" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Geant4/Geant4LibraryDepends.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Geant4/Geant4LibraryDepends.cmake"
         "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/CMakeFiles/Export/2561fa484b1aba2acff7bf3a9f332413/Geant4LibraryDepends.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Geant4/Geant4LibraryDepends-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Geant4/Geant4LibraryDepends.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Geant4" TYPE FILE MESSAGE_LAZY FILES "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/CMakeFiles/Export/2561fa484b1aba2acff7bf3a9f332413/Geant4LibraryDepends.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Geant4" TYPE FILE MESSAGE_LAZY FILES "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/CMakeFiles/Export/2561fa484b1aba2acff7bf3a9f332413/Geant4LibraryDepends-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Development" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Geant4" TYPE FILE MESSAGE_LAZY FILES
    "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/InstallTreeFiles/Geant4Config.cmake"
    "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/Geant4ConfigVersion.cmake"
    "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/G4EXPATShim.cmake"
    "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/G4FreetypeShim.cmake"
    "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/G4HDF5Shim.cmake"
    "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/G4MotifShim.cmake"
    "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/G4X11Shim.cmake"
    "/Users/ksloan/geant4-v11.3.2/cmake/Templates/UseGeant4.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Development" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Geant4" TYPE DIRECTORY MESSAGE_LAZY FILES "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/Modules")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Examples" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/Geant4" TYPE DIRECTORY MESSAGE_LAZY FILES "/Users/ksloan/geant4-v11.3.2/examples" REGEX "/cvs$" EXCLUDE REGEX "/\\.svn$" EXCLUDE)
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/ksloan/github/FreeCAD_GDML_Test_Collection/MH_G4_Example/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()

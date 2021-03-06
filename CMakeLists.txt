####################################################################################################
# CMakeLists file for Trax
# Generated by Gianluca Elia
# 2020-04-11
####################################################################################################

####################################################################################################
# basic project config
cmake_minimum_required(VERSION 3.5)
set(project_name "Trax")
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 11)

####################################################################################################
# load modules
include(SuperColliderServerPlugin RESULT_VARIABLE server_plugin_found)
if(NOT server_plugin_found)
    message(FATAL_ERROR "Could not find server plugin functions module")
endif()

include(SuperColliderCompilerConfig RESULT_VARIABLE compiler_config_found)
if(NOT compiler_config_found)
    message(FATAL_ERROR "Could not find compiler config module")
endif()

# Windows - puts redistributable DLLs in install directory
include(InstallRequiredSystemLibraries)

sc_check_sc_path("${SC_PATH}")
message(STATUS "Found SuperCollider: ${SC_PATH}")
set(SC_PATH "${SC_PATH}" CACHE PATH
    "Path to SuperCollider source. Relative paths are treated as relative to this script" FORCE)

include("${SC_PATH}/SCVersion.txt")
message(STATUS "Building plugins for SuperCollider version: ${SC_VERSION}")

# set project here to avoid SCVersion.txt clobbering our version info
project(${project_name})
sc_do_initial_compiler_config() # do after setting project so compiler ID is available

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR NOT CMAKE_INSTALL_PREFIX)
    message(WARNING "No install prefix provided, defaulting to $BUILD_DIR/install")
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install prefix" FORCE)
endif()

message(STATUS "Install directory set to: ${CMAKE_INSTALL_PREFIX}")

####################################################################################################
# options
option(SUPERNOVA "Build plugins for supernova" ON)
option(SCSYNTH "Build plugins for scsynth" ON)
option(NATIVE "Optimize for native architecture" OFF)
option(STRICT "Use strict warning flags" OFF)
option(NOVA_SIMD "Build plugins with nova-simd support." ON)

####################################################################################################
# include libraries

if (NOVA_SIMD)
	add_definitions(-DNOVA_SIMD)
	include_directories(${SC_PATH}/external_libraries/nova-simd)
endif()

####################################################################################################
# Begin target Trax

set(Trax_cpp_files
    plugins/Trax/Trax.cpp
    plugins/Trax/Trax.hpp
)
set(Trax_sc_files
    plugins/Trax/Trax.sc
)
set(Trax_schelp_files
    plugins/Trax/Trax.schelp
)

sc_add_server_plugin(
    "Trax/Trax" # desination directory
    "Trax" # target name
    "${Trax_cpp_files}"
    "${Trax_sc_files}"
    "${Trax_schelp_files}"
)

# End target Trax
####################################################################################################

####################################################################################################
# Begin target HPS

set(HPS_cpp_files
    plugins/HPS/HPS.hpp
    plugins/HPS/HPS.cpp
)
set(HPS_sc_files
    plugins/HPS/HPS.sc
)
set(HPS_schelp_files
    plugins/HPS/HPS.schelp
)

sc_add_server_plugin(
    "Trax/HPS" # desination directory
    "HPS" # target name
    "${HPS_cpp_files}"
    "${HPS_sc_files}"
    "${HPS_schelp_files}"
)

# End target HPS
####################################################################################################

####################################################################################################
# Begin target PV_ShiftConv

set(PV_ShiftConv_cpp_files
    plugins/PV_ShiftConv/PV_ShiftConv.hpp
    plugins/PV_ShiftConv/PV_ShiftConv.cpp
)
set(PV_ShiftConv_sc_files
    plugins/PV_ShiftConv/PV_ShiftConv.sc
)
set(PV_ShiftConv_schelp_files
    plugins/PV_ShiftConv/PV_ShiftConv.schelp
)

sc_add_server_plugin(
    "Trax/PV_ShiftConv" # desination directory
    "PV_ShiftConv" # target name
    "${PV_ShiftConv_cpp_files}"
    "${PV_ShiftConv_sc_files}"
    "${PV_ShiftConv_schelp_files}"
)

# End target PV_ShiftConv
####################################################################################################

####################################################################################################
# Begin target PV_TPV

set(PV_TPV_cpp_files
    plugins/PV_TPV/PV_TPV.cpp
    plugins/PV_TPV/PV_TPV.hpp
)
set(PV_TPV_sc_files
    plugins/PV_TPV/PV_TPV.sc
)
set(PV_TPV_schelp_files
    plugins/PV_TPV/PV_TPV.schelp
)

sc_add_server_plugin(
    "Trax/PV_TPV" # desination directory
    "PV_TPV" # target name
    "${PV_TPV_cpp_files}"
    "${PV_TPV_sc_files}"
    "${PV_TPV_schelp_files}"
)

# End target PV_TPV
####################################################################################################

####################################################################################################
# Begin target PeaX

set(PeaX_cpp_files
    plugins/PeaX/PeaX.cpp
    plugins/PeaX/PeaX.hpp
)
set(PeaX_sc_files
    plugins/PeaX/PeaX.sc
)
set(PeaX_schelp_files
    plugins/PeaX/PeaX.schelp
)

sc_add_server_plugin(
    "Trax/PeaX" # desination directory
    "PeaX" # target name
    "${PeaX_cpp_files}"
    "${PeaX_sc_files}"
    "${PeaX_schelp_files}"
)

# End target PeaX
####################################################################################################

####################################################################################################
# Begin target PeaXShift

set(PeaXShift_cpp_files
    plugins/PeaXShift/PeaXShift.hpp
    plugins/PeaXShift/PeaXShift.cpp
)
set(PeaXShift_sc_files
    plugins/PeaXShift/PeaXShift.sc
)
set(PeaXShift_schelp_files
    plugins/PeaXShift/PeaXShift.schelp
)

sc_add_server_plugin(
    "Trax/PeaXShift" # desination directory
    "PeaXShift" # target name
    "${PeaXShift_cpp_files}"
    "${PeaXShift_sc_files}"
    "${PeaXShift_schelp_files}"
)

# End target PeaXShift
####################################################################################################

####################################################################################################
# Begin target PeaXShift2

set(PeaXShift2_cpp_files
    plugins/PeaXShift2/PeaXShift2.hpp
    plugins/PeaXShift2/PeaXShift2.cpp
)
set(PeaXShift2_sc_files
    plugins/PeaXShift2/PeaXShift2.sc
)
set(PeaXShift2_schelp_files
    plugins/PeaXShift2/PeaXShift2.schelp
)

sc_add_server_plugin(
    "Trax/PeaXShift2" # desination directory
    "PeaXShift2" # target name
    "${PeaXShift2_cpp_files}"
    "${PeaXShift2_sc_files}"
    "${PeaXShift2_schelp_files}"
)

# End target PeaXShift2
####################################################################################################

####################################################################################################
# END PLUGIN TARGET DEFINITION
####################################################################################################

message(STATUS "Generating plugin targets done")

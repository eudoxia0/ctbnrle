cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)
project(GLPK C)

get_filename_component(GLPK_ROOT_DIR glpk-4.47 ABSOLUTE CACHE INTERNAL "GLPK root directory")

set(GLPK_INCLUDE_DIRS ${GLPK_ROOT_DIR}/src CACHE INTERNAL "GLPK include directories")

set(GLPK_LIBRARIES ${GLPK_ROOT_DIR}/src/.libs/libglpk.a CACHE INTERNAL "Do not edit this, but do run configure and make from within the GLPK directory before attempting to build CTBN.")

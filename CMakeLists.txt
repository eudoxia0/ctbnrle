cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

project(CTBNRLE)

enable_testing()

add_subdirectory(glpk)

set(CTBN_INCLUDE_DIRS
    ${GLPK_INCLUDE_DIRS}
    ${EIGEN_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/hdr
)

add_subdirectory(src)
add_subdirectory(tst)
add_subdirectory(demo)


# (C) Copyright 2019 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_FORTRAN_STANDARD 08)
set(CMAKE_FORTRAN_STANDARD_REQUIRED ON)
set(CMAKE_FORTRAN_EXTENSIONS OFF)

if( NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_definitions( -DNDEBUG )
endif( )

#######################################################################################
# Fortran
#######################################################################################

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
  include( compiler_flags_GNU_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM" )
  include( compiler_flags_IntelLLVM_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Intel" )
  include( compiler_flags_Intel_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "XL" )
  include( compiler_flags_XL_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
  include( compiler_flags_Cray_Fortran )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
  include( compiler_flags_NVHPC_Fortran )
else()
  message( STATUS "Fortran compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")
endif()

#######################################################################################
# C
#######################################################################################

# todo

#######################################################################################
# C++
#######################################################################################

if( CMAKE_CXX_COMPILER_ID MATCHES "GNU" )
  include( compiler_flags_GNU_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
  include( compiler_flags_Intel_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "XL" )
  include( compiler_flags_XL_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Cray" )
  include( compiler_flags_Cray_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  include( compiler_flags_Clang_CXX )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "NVHPC" )
  include( compiler_flags_NVHPC_CXX )
else()
  message( STATUS "C++ compiler with ID ${CMAKE_CXX_COMPILER_ID} will be used with CMake default options")
endif()

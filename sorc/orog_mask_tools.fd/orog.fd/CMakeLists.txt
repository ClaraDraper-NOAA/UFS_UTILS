set(lib_src netcdf_io.F90)
set(exe_src mtnlm7_oclsm.f)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -convert big_endian -assume byterecl")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fconvert=big-endian -fno-range-check")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()
endif()

set(exe_name orog)

add_library(orog_lib STATIC ${lib_src})
add_executable(${exe_name} mtnlm7_oclsm.f)

set(mod_dir "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(orog_lib PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
target_include_directories(orog_lib INTERFACE ${mod_dir})

target_link_libraries(
  orog_lib
  PUBLIC
  bacio::bacio_4
  w3nco::w3nco_d
  ip::ip_d
  sp::sp_d
  NetCDF::NetCDF_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(orog_lib PUBLIC OpenMP::OpenMP_Fortran)
endif()

target_link_libraries(orog PRIVATE orog_lib)

install(TARGETS ${exe_name} RUNTIME DESTINATION ${exec_dir})

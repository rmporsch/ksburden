project(ksburden)
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
include(ExternalProject)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
include_directories(inc)
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )

# ==========
# Check compiler version
# ==========
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
		message(FATAL_ERROR "GCC version must be at least 4.8!")
	endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

	if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.2)
		message(FATAL_ERROR "Clang version must be at least 3.2!")
	endif()

else()
	message(WARNING "You are using an unsupported compiler!
	Compilation has only been tested with Clang and GCC.")
endif()


# ================
# External projects
# ================

# gflags
ExternalProject_Add(gflagsDown
	GIT_REPOSITORY https://github.com/gflags/gflags.git
       	PREFIX ${LIBRARY_OUTPUT_PATH}/gflags
	BUILD_IN_SOURCE 1
	CONFIGURE_COMMAND export CXXFLAGS="-fPIC" && cmake .
	BUILD_COMMAND make
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	)

set(gflags_DIR ${LIBRARY_OUTPUT_PATH}/gflags/src/gflagsDown)
add_library(gflags STATIC IMPORTED)
add_dependencies(gflags gflagsDown)
set_property(TARGET gflags PROPERTY IMPORTED_LOCATION ${gflags_DIR}/lib/libgflags.a)
include_directories(${gflags_DIR}/include/)

#============
# LIBSTATGEN
#============
ExternalProject_Add(libStatGen
	GIT_REPOSITORY https://github.com/statgen/libStatGen.git
	PREFIX ${LIBRARY_OUTPUT_PATH}/libStatGen
	CONFIGURE_COMMAND ""
	BUILD_IN_SOURCE 1
	BUILD_COMMAND make CFLAGS='-fPIC'
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	)

add_library(libStatGenLib STATIC IMPORTED)
set_property(TARGET libStatGenLib PROPERTY IMPORTED_LOCATION ${LIBRARY_OUTPUT_PATH}/libStatGen/src/libStatGen/libStatGen.a)
add_dependencies(libStatGenLib libStatGen)
include_directories(${LIBRARY_OUTPUT_PATH}/libStatGen/src/libStatGen/include/)

#============
# easylogging
#============
ExternalProject_Add(easylogging
	PREFIX ${LIBRARY_OUTPUT_PATH}/easylogging
	DOWNLOAD_COMMAND wget https://github.com/easylogging/easyloggingpp/releases/download/v9.80/easyloggingpp_v9.80.tar.gz
	CONFIGURE_COMMAND tar -zxvf ../easyloggingpp_v9.80.tar.gz 
	BUILD_IN_SOURCE 1
	BUILD_COMMAND ""
	INSTALL_COMMAND ""
	UPDATE_COMMAND ""
	)

include_directories(${LIBRARY_OUTPUT_PATH}/easylogging/src/easylogging)


#============
# OpenMP
#============ 
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#========
# Boost
#========
set(BOOST_DIR ${LIBRARY_OUTPUT_PATH}/boost)
include_directories(BOOST_DIR)

#============
# Armadillo
#============
ExternalProject_Add(ArmadilloDownload
	PREFIX ${LIBRARY_OUTPUT_PATH}/armadillo-7.100.3
	DOWNLOAD_COMMAND wget http://sourceforge.net/projects/arma/files/armadillo-8.100.0.tar.xz
	CONFIGURE_COMMAND tar -xf ../armadillo-7.100.3.tar.xz -C ../ArmadilloDownload  --strip-components=1
	BUILD_IN_SOURCE 1
	BUILD_COMMAND cmake .
	INSTALL_COMMAND make
	UPDATE_COMMAND ""
	)

set(Armadillo_DIR ${LIBRARY_OUTPUT_PATH}/armadillo-7.100.3/src/ArmadilloDownload)
set(ARMADILLO_LIBRARY ${Armadillo_DIR}/libarmadillo.so)
set(ARMADILLO_INCLUDE_DIR ${Armadillo_DIR}/include/)
include_directories(${ARMADILLO_INCLUDE_DIR})
#=============================
#=============================

#=============================
# Library
#=============================
# loadVariantFile
add_library(variantFile SHARED "src/loadVarinatFile.cpp") 
add_dependencies(variantFile easylogging) 
target_compile_features(variantFile PRIVATE cxx_range_for)
# models
add_library(models SHARED "src/models.cpp") 
add_dependencies(models ArmadilloDownload easylogging) 
target_compile_features(models PRIVATE cxx_range_for)

# load_vcf
add_library(loadvcf SHARED "src/load_vcf.cpp")
add_dependencies(loadvcf ArmadilloDownload libStatGen easylogging variantFile)
target_link_libraries(loadvcf libStatGenLib)
target_compile_features(loadvcf PRIVATE cxx_range_for)

# load_plink
add_library(loadplink SHARED "src/load_plink.cpp")
add_dependencies(loadplink ArmadilloDownload libStatGen easylogging variantFile)
target_compile_features(loadplink PRIVATE cxx_range_for)

# liability_model
add_library(liabilitymodel SHARED "src/liability_model.cpp")
add_dependencies(liabilitymodel ArmadilloDownload loadvcf loadplink)
target_link_libraries(liabilitymodel loadvcf)
target_compile_features(liabilitymodel PRIVATE cxx_range_for)

# simulation_model
add_library(simulation SHARED "src/simulation.cpp")
add_dependencies(simulation ArmadilloDownload loadvcf models liabilitymodel)
target_link_libraries(simulation loadvcf models liabilitymodel)
target_compile_features(simulation PRIVATE cxx_range_for)

# analysis
add_library(analysis SHARED "src/data_analysis.cpp")
add_dependencies(analysis ArmadilloDownload loadvcf)
target_link_libraries(analysis loadvcf)
target_compile_features(analysis PRIVATE cxx_range_for)


# sim_flags
add_library(sim_flags SHARED "src/simulation_flags.cpp")
add_dependencies(sim_flags gflags)
target_link_libraries(sim_flags gflags pthread)
target_compile_features(sim_flags PRIVATE cxx_range_for)
# see https://github.com/halide/Halide/issues/175 why I added 'pthread'

# ksburden_flags
add_library(ksburden_flags SHARED "src/flags_ksburden.cpp")
add_dependencies(ksburden_flags gflags)
target_link_libraries(ksburden_flags gflags pthread)
target_compile_features(ksburden_flags PRIVATE cxx_range_for)
# see https://github.com/halide/Halide/issues/175 why I added 'pthread'

install(TARGETS models loadvcf liabilitymodel sim_flags analysis DESTINATION lib)

#=========
# programs
#=========

# simulation program
add_executable(sim ${CMAKE_SOURCE_DIR}/src/sim.cpp)
target_link_libraries(sim
	simulation
	models
	liabilitymodel
	loadvcf
  loadplink
  variantFile
	${ARMADILLO_LIBRARIES}
	z
	sim_flags
	-O2 blas lapack -DEBUG
	)
target_compile_features(sim PRIVATE cxx_range_for)

add_executable(sim2 ${CMAKE_SOURCE_DIR}/src/sim_2.cpp)
target_link_libraries(sim2
	simulation
	models
	liabilitymodel
	loadvcf
  loadplink
  variantFile
	${ARMADILLO_LIBRARIES}
	z
	sim_flags
	-O2 blas lapack -DEBUG
	)
target_compile_features(sim2 PRIVATE cxx_range_for)

# analysis program
add_executable(ksburden ${CMAKE_SOURCE_DIR}/src/ksburden.cpp)
target_link_libraries(ksburden
	models
	analysis
	loadvcf
  variantFile
	${ARMADILLO_LIBRARIES}
	z
	ksburden_flags
	-O2 blas lapack
	)
target_compile_features(ksburden PRIVATE cxx_range_for)

#tests
add_executable(tests ${CMAKE_SOURCE_DIR}/src/tests/run_tests.cpp)

# gflags
ExternalProject_Add(googleTest_git
	GIT_REPOSITORY https://github.com/google/googletest.git
       	PREFIX ${LIBRARY_OUTPUT_PATH}/googletest
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_SOURCE_DIR}
	UPDATE_COMMAND ""
	)

ExternalProject_Get_Property(googleTest_git install_dir)
include_directories(${install_dir}/include)
set(gtest_SOURCE_DIR ${LIBRARY_OUTPUT_PATH}/googletest/src/googleTest_git/googletest)
include_directories(
  ${gtest_SOURCE_DIR}/include
  ${gtest_SOURCE_DIR})
set(gtest_main ${LIBRARY_OUTPUT_PATH}/libgtest_main.a)
set(gtest ${LIBRARY_OUTPUT_PATH}/libgtest.a)

add_dependencies(tests googleTest_git)
target_link_libraries(tests
	simulation
	models
	analysis
	liabilitymodel
	loadvcf
  loadplink
  variantFile
	${ARMADILLO_LIBRARIES}
	z
	${gtest}
	-O2 blas lapack -DEBUG
	)
target_compile_features(tests PRIVATE cxx_range_for)

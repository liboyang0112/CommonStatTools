# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
set (CMAKE_CXX_STANDARD 14)
else()
set (CMAKE_CXX_STANDARD 17)
endif()

project(ExoStats)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Graf Hist Gpad Matrix MathCore RIO Net Tree RooFitCore RooFit RooStats Minuit HistFactory )

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
if(DEFINED ENV{ROOT_USE_FILE})
  include(${ROOT_USE_FILE})
endif()

include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/cmake/modules)
include(RootNewMacros OPTIONAL) #only for root < 6.20.02
#ROOT_GENERATE_DICTIONARY(G__ExoStats Event.h LINKDEF EventLinkDef.h)
#ROOT_GENERATE_DICTIONARY(G__ExoStats Minimization.h AsimovDataMaking.h AsymptoticsCLsRunner.h HistFactoryInspector.h RooExpandedFitResult.h LINKDEF LinkDef.h)
ROOT_GENERATE_DICTIONARY(G__ExoStats RooExpandedFitResult.h LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(ExoStats SHARED Minimization.C AsimovDataMaking.C AsymptoticsCLsRunner.C HistFactoryInspector.C RooExpandedFitResult.C FitCrossCheckForLimits.C G__ExoStats.cxx runSig.C runAsymptoticsCLs.C)
target_link_libraries(ExoStats ${ROOT_LIBRARIES})

if(DEFINED ENV{AtlasProject})
  add_dependencies(ExoStats StatisticsTools)
endif()

#---Create  a main program using the library
#add_executable(runAsymptoticsCLs runAsymptoticsCLs.C)
#target_link_libraries(runAsymptoticsCLs ExoStats)

# install the library
install(TARGETS ExoStats
        LIBRARY DESTINATION lib)

# install pcm and rootmap into lib
if(DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY AND NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY STREQUAL "")
  set(pcm_output_dir ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
else()
  set(pcm_output_dir ${CMAKE_CURRENT_BINARY_DIR})
endif()
install(FILES ${pcm_output_dir}/libExoStats_rdict.pcm DESTINATION lib)
install(FILES ${pcm_output_dir}/libExoStats.rootmap DESTINATION lib)

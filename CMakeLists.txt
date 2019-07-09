# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
set (CMAKE_CXX_STANDARD 14)

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
include(RootNewMacros) 
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


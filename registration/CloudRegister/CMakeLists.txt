cmake_minimum_required(VERSION 2.8)
project(cloud_regist)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall  -O3 -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall   -O3 -march=native")

ADD_DEFINITIONS(-DUBUNTU_SWITCH)

# Check C++11 or C++0x support
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX14)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
   add_definitions(-DCOMPILEDWITHC14)
   message(STATUS "Using flag -std=c++14.")
elseif(COMPILER_SUPPORTS_CXX11)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   add_definitions(-DCOMPILEDWITHC11)
   message(STATUS "Using flag -std=c++11.")
elseif(COMPILER_SUPPORTS_CXX0X)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
   add_definitions(-DCOMPILEDWITHC0X)
   message(STATUS "Using flag -std=c++0x.")
else()
   message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

#find_package(OpenCV 3.0 QUIET)
#if(NOT OpenCV_FOUND)
   #find_package(OpenCV 2.4.3 QUIET)
   #if(NOT OpenCV_FOUND)
      #message(FATAL_ERROR "OpenCV > 2.4.3 not found.")
   #endif()
#endif()

find_package(Eigen3 3.1.0 REQUIRED)
#find_package(Pangolin REQUIRED)
find_package(PCL REQUIRED)
find_package(glog REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})

link_directories(${PROJECT_SOURCE_DIR}/G2O/g2o_ubuntu/)

add_definitions(${PCL_DEFINITIONS})

include_directories(
${PROJECT_SOURCE_DIR}
${PROJECT_SOURCE_DIR}/G2O/include/g2o/
${PROJECT_SOURCE_DIR}/CloudRegister
${PROJECT_SOURCE_DIR}/CloudExample
#${EIGEN3_INCLUDE_DIR}
#${Pangolin_INCLUDE_DIRS}
)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_library(${PROJECT_NAME} SHARED
CloudRegister/CADModel.cpp  
CloudRegister/CloudRegister.cpp   
CloudRegister/EdgeConstraint.cpp  
CloudRegister/funHelper.cpp  
CloudRegister/GeometryUtils.cpp  
CloudRegister/TransformOptimize.cpp
CloudRegister/CloudSegment.cpp
CloudRegister/SimpleViewer.cpp
CloudRegister/CalcNetHeight.cpp
CloudRegister/CalcHoleMeasure.cpp
CloudRegister/CalcBayAndDepthMeasure.cpp
CloudRegister/CalcMeasureHelper.cpp
CloudRegister/CalcWallVerticality.cpp
CloudRegister/CalcRootFlatness.cpp
CloudRegister/Threshold.cpp
CloudRegister/CalcWallFlatness.cpp
CloudRegister/CalcCorner.cpp
CloudRegister/CalcSquareness.cpp
)

target_link_libraries(${PROJECT_NAME}
${PCL_LIBRARIES}
g2o_core
g2o_solver_eigen
g2o_types_data
g2o_stuff
g2o_types_sba
g2o_types_sim3
#${OpenCV_LIBS}
#${EIGEN3_LIBS}
#${Pangolin_LIBRARIES}
)

# Build examples

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/CloudExample)

add_executable(cloud_examp 
CloudExample/CloudExample.cpp) 
target_link_libraries(cloud_examp ${PROJECT_NAME} glog::glog) #${PCL_LIBRARIES})


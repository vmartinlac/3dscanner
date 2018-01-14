cmake_minimum_required(VERSION 3.0)

project(3dscanner CXX)

include(GNUInstallDirs)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#find_package(MPI COMPONENTS CXX REQUIRED)
#find_package(Qt5 COMPONENTS Core Widgets Gui PrintSupport Test Multimedia REQUIRED)
#find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
#find_package(Boost 1.60.0 REQUIRED COMPONENTS math_tr1 random program_options)
#find_package(VTK REQUIRED)

#include(${VTK_USE_FILE})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/qcustomplot")

#enable_testing()

add_executable(
   3dscanner
   src/RotatingLaserPlaneSolver.cpp
   src/main.cpp)
target_link_libraries(3dscanner Eigen3::Eigen)


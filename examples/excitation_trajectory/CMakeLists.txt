#
cmake_minimum_required(VERSION 3.0)

project(excitation_trajectory)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_executable(excitation_trajectory
  src/excitation_trajectory.cpp
)
target_include_directories(excitation_trajectory
  PRIVATE 
  /usr/include
)
target_link_libraries(excitation_trajectory
  librobotdynamics.so 
  libnlopt.so 
  libm.so
)





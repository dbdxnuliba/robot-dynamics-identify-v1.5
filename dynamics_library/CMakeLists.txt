#
cmake_minimum_required(VERSION 3.0)

project(dynamics_library)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_library(robotdynamics SHARED 
  src/robot_model.cpp
  src/robot_identify.cpp
)
target_include_directories(robotdynamics
  PRIVATE 
  include
)

install(DIRECTORY include/
  DESTINATION /usr/include
)

install(TARGETS robotdynamics
  LIBRARY DESTINATION /usr/lib
)






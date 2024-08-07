cmake_minimum_required(VERSION 3.0.2)
project(switchlib)

include_directories(
  include
  ${EIGEN_INCLUDE_DIR}
  "/usr/include/eigen3"
  ${G2O_INCLUDE_DIR}
)

add_library(switchlib SHARED 
src/edge_se2Switchable.cpp
src/edge_se3Switchable.cpp
src/edge_switchPrior.cpp
src/jacobian.cpp
src/vertex_switchLinear.cpp)

target_link_libraries(switchlib
${G2O_TYPES_SLAM2D}
${G2O_TYPES_SLAM3D}
${G2O_TYPES_DATA}
${G2O_CORE_LIBRARY}
${G2O_STUFF_LIBRARY}
${G2O_SOLVER_EIGEN}
)
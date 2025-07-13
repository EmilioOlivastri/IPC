cmake_minimum_required(VERSION 3.5.2)
project(IPC)

set(CMAKE_CXX_FLAGS "-std=c++14 -O3")
set(G2O_ROOT /home/slam-emix/Workspace/lib/g2o)
set(CMAKE_BUILD_TYPE Release)
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(Boost REQUIRED COMPONENTS system)
find_package(Glog REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(G2O REQUIRED)
find_package(yaml-cpp REQUIRED)

include_directories(
  include
  ${EIGEN_INCLUDE_DIR}
  "/usr/include/eigen3"
  ${G2O_INCLUDE_DIR}
  ${YAML_CPP_INCLUDE_DIR}
  ${GLOG_INCLUDE_DIR}/glog
)

add_executable(ipc_tester_2D
examples/ipc_tester_2D.cpp
src/utils.cpp
src/consensus_utils.cpp
src/consensus.cpp
src/simulation.cpp
)

target_link_libraries(ipc_tester_2D
  ${G2O_TYPES_SLAM2D}
  ${G2O_TYPES_SLAM3D}
  ${G2O_TYPES_DATA}
  ${G2O_CORE_LIBRARY}
  ${G2O_STUFF_LIBRARY}
  ${G2O_SOLVER_EIGEN}
  ${YAML_CPP_LIBRARIES}
  ${GLOG_LIBRARY}
 )


add_executable(ipc_tester_3D
examples/ipc_tester_3D.cpp
src/utils.cpp
src/consensus_utils.cpp
src/consensus.cpp
src/simulation.cpp
)

target_link_libraries(ipc_tester_3D
  ${G2O_TYPES_SLAM2D}
  ${G2O_TYPES_SLAM3D}
  ${G2O_TYPES_DATA}
  ${G2O_CORE_LIBRARY}
  ${G2O_STUFF_LIBRARY}
  ${G2O_SOLVER_EIGEN}
  ${YAML_CPP_LIBRARIES}
  ${GLOG_LIBRARY}
 )


add_executable(graph_fixer
examples/graph_fixer.cpp
src/utils.cpp
src/consensus_utils.cpp
src/consensus.cpp
)

target_link_libraries(graph_fixer
  ${G2O_TYPES_SLAM2D}
  ${G2O_TYPES_SLAM3D}
  ${G2O_TYPES_DATA}
  ${G2O_CORE_LIBRARY}
  ${G2O_STUFF_LIBRARY}
  ${G2O_SOLVER_EIGEN}
  ${YAML_CPP_LIBRARIES}
  ${GLOG_LIBRARY}
 )

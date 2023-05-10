cmake_minimum_required(VERSION 3.3)

project(LEARNING)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()


find_package(dlib REQUIRED)


include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scip_learning
   src/cppmain.cpp
   src/nodesel_random.cpp src/branch_unrealistic.cpp src/branch_unrealistic.h src/Utils.cpp src/Utils.h src/DatasetWriter.cpp src/DatasetWriter.h src/FeaturesCalculator.cpp src/FeaturesCalculator.h src/dialog_generateDataset.cpp src/dialog_generateDataset.h src/RegressionModel.cpp src/RegressionModel.h src/branch_unrealisticTrained.cpp src/branch_unrealisticTrained.h src/dialog_loadModel.cpp src/dialog_loadModel.h src/dialog_trainmodel.cpp src/dialog_trainmodel.h src/mpi/Worker.cpp src/mpi/Worker.h src/EventhdlrUpdateFeatures.cpp src/EventhdlrUpdateFeatures.h src/BranchingHistory.cpp src/BranchingHistory.h)

target_link_libraries(scip_learning ${SCIP_LIBRARIES})
target_link_libraries(scip_learning  dlib::dlib)
target_link_libraries(scip_learning ${MPI_LIBRARIES})

find_package(ZLIB)
if(ZLIB_FOUND)
    target_link_libraries(scip_learning ${ZLIB_LIBRARIES})
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CPP_FLAGS "${CMAKE_CPP_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

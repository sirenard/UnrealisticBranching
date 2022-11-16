cmake_minimum_required(VERSION 3.3)

project(LEARNING)



if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

#set(DLIB_PATH /home/simon/Téléchargements/dlib-19.24)
find_package(dlib REQUIRED)
#include_directories(${DLIB_PATH})
#target_link_libraries(scip_learning ${DLIB_PATH} dlib)
#add_subdirectory(${DLIB_PATH} dlib_build)


include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scip_learning
   src/cppmain.cpp
   src/nodesel_random.cpp src/branch_unrealistic.cpp src/branch_unrealistic.h src/Utils.cpp src/Utils.h src/DatasetWriter.cpp src/DatasetWriter.h src/FeaturesCalculator.cpp src/FeaturesCalculator.h src/dialog_generateDataset.cpp src/dialog_generateDataset.h src/RegressionModel.cpp src/RegressionModel.h src/branch_unrealisticTrained.cpp src/branch_unrealisticTrained.h src/dialog_loadModel.cpp src/dialog_loadModel.h src/dialog_trainmodel.cpp src/dialog_trainmodel.h src/mpi/Worker.cpp src/mpi/Worker.h)

target_link_libraries(scip_learning ${SCIP_LIBRARIES})
target_link_libraries(scip_learning  dlib::dlib)

if(ZLIB)
    find_package(ZLIB)
    if(ZLIB_FOUND)
        target_link_libraries(scip_learning ${ZLIB_LIBRARIES})
    endif()
endif()

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_CPP_COMPILER mpicxx)
SET(CMAKE_CXX_COMPILER mpicxx)

if( TARGET examples )
    add_dependencies( examples scip_learning )
endif()


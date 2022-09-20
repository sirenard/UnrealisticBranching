cmake_minimum_required(VERSION 3.3)

project(LEARNING)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(/home/simon/Téléchargements/OpenNN}) # Voir si ça marche
include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scip_learning
   src/cppmain.cpp
   src/nodesel_random.cpp src/branch_unrealistic.cpp src/branch_unrealistic.h src/Utils.cpp src/Utils.h src/DatasetWriter.cpp src/DatasetWriter.h src/FeaturesCalculator.cpp src/FeaturesCalculator.h)

target_link_libraries(scip_learning ${SCIP_LIBRARIES})

if(ZLIB)
    find_package(ZLIB)
    if(ZLIB_FOUND)
        target_link_libraries(scip_learning ${ZLIB_LIBRARIES})
    endif()
endif()

if( TARGET examples )
    add_dependencies( examples scip_learning )
endif()


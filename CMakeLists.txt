cmake_minimum_required(VERSION 3.14..3.16)
project(placer)

set(CMAKE_CXX_STANDARD 17)
set(SOURCES
        src/circuit/Circuit.cpp
        src/circuit/Parser.cpp
        src/circuit/Evaluator.cpp

        src/dataStructures/Die.cpp
        src/dataStructures/Instance.cpp
        src/dataStructures/Net.cpp
        src/dataStructures/Pin.cpp

        src/algorithm/place.cpp
        src/algorithm/FFT_calculator.cpp
        src/algorithm/math/matrixSolver.cpp

        src/tools/visualizer.cpp

        src/etc/example.cpp
        src/etc/dbExample.cpp
        )

set(INCLUDE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/include/algorithm
        ${PROJECT_SOURCE_DIR}/include/circuit
        ${PROJECT_SOURCE_DIR}/include/dataStructures
        ${PROJECT_SOURCE_DIR}/include/tools
        ${PROJECT_SOURCE_DIR}/submodule/CImg
        )

# Find OpenMP
find_package(OpenMP)

add_executable(placer src/main_placer.cpp ${SOURCES})
add_executable(qPlacer src/main_qPlacer.cpp ${SOURCES})
add_executable(evaluator src/main_evaluator.cpp ${SOURCES})

target_include_directories(placer PRIVATE ${INCLUDE_DIRS})
target_include_directories(qPlacer PRIVATE ${INCLUDE_DIRS})
target_include_directories(evaluator PRIVATE ${INCLUDE_DIRS})

add_subdirectory(submodule/OpenDB)
# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(placer PRIVATE OpenMP::OpenMP_CXX opendb X11 pthread)
    target_link_libraries(qPlacer PRIVATE OpenMP::OpenMP_CXX opendb X11 pthread)
    target_link_libraries(evaluator PRIVATE OpenMP::OpenMP_CXX opendb X11 pthread)
else()
    target_link_libraries(placer PRIVATE opendb X11 pthread)
    target_link_libraries(qPlacer PRIVATE opendb X11 pthread)
    target_link_libraries(evaluator PRIVATE opendb X11 pthread)
endif()
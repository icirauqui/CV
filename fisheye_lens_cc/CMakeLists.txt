cmake_minimum_required(VERSION 3.5)

project(lense LANGUAGES CXX)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall -O3 -march=native ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native -std=c++17")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w -O3 -march=native -std=c++17") 
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -w -O3 -march=native") 

# Find package OpenCV, needed because new opencv4 is in use
find_package(OpenCV REQUIRED)

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${OpenCV_INCLUDE_DIRS}
)

file(CREATE_LINK "${PROJECT_SOURCE_DIR}/images"
        "${CMAKE_CURRENT_BINARY_DIR}/images" SYMBOLIC)

#file(CREATE_LINK "${PROJECT_SOURCE_DIR}/images"
#        "${CMAKE_CURRENT_BINARY_DIR}/images" SYMBOLIC)

add_executable(lense
  src/fe_lens/fe_lens.cpp
  src/fe_lens/fe_lens.hpp
  src/main.cpp
)

target_link_libraries(lense
    ${OpenCV_LIBS}
)


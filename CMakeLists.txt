cmake_minimum_required( VERSION 2.8 )
project (mystim)

set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_BUILD_TYPE "Release" )
set( CMAKE_CXX_FLAGS "-std=c++11 -march=native -O3" )

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )
set( LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib )

# 寻找OpenCV库
find_package( OpenCV 3.1 REQUIRED )
# 添加头文件
include_directories( ${OpenCV_INCLUDE_DIRS} )
# Eigen
include_directories( "/usr/include/eigen3" )

set( THIRD_PARTY_LIBS 
    ${OpenCV_LIBS}
)

include_directories( ${PROJECT_SOURCE_DIR}/include )
add_subdirectory( src )
add_subdirectory( test )



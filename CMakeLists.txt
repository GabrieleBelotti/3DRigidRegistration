cmake_minimum_required(VERSION 3.14)

project(ITK3DRegClass)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

#add_executable(ImageRegistration3DClass UserMainRegistration.cpp _gbe_3DRegistration.cpp _gbe_3DRegistration.h)
add_executable(ImageRegistration3DClass UserMainRegistration.cpp)
target_include_directories(ImageRegistration3DClass PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ImageRegistration3DClass ${ITK_LIBRARIES})
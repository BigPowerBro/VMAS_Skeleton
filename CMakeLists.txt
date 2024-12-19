# 设置 CMake 的最低版本要求
cmake_minimum_required(VERSION 3.10)

# 设置项目名称
project(VMAS_skeleton)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 11)

# 设置 OpenMesh 和 Eigen 的路径
set(OPENMESH_DIR "package/OpenMesh")
set(EIGEN3_INCLUDE_DIR "package/Eigen3")

# 添加 include 目录
include_directories(${OPENMESH_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)

# 添加 lib 目录
link_directories(${CMAKE_SOURCE_DIR}/lib)

# 添加源代码目录
file(GLOB SOURCES "src/*.cpp")
file(GLOB INCLUDES "include/*.h")

# 查找 OpenMesh 库
find_library(OPENMESH_CORE_LIBRARY OpenMeshCored PATHS ${CMAKE_SOURCE_DIR}/lib)
find_library(OPENMESH_TOOLS_LIBRARY OpenMeshToolsd PATHS ${CMAKE_SOURCE_DIR}/lib)

# 创建可执行文件
add_executable(VMAS_skeleton ${SOURCES} ${INCLUDES})

target_compile_definitions(VMAS_skeleton PRIVATE _USE_MATH_DEFINES)

# 链接 OpenMesh 库
target_link_libraries(VMAS_skeleton
    ${OPENMESH_CORE_LIBRARY}
    ${OPENMESH_TOOLS_LIBRARY}
)
SET(target final)

add_executable(${target} main.cpp sim.cpp sync.cpp glad.c)
TARGET_LINK_LIBRARIES(${target})

set_property(TARGET ${target} PROPERTY CXX_STANDARD 17)
set_property(TARGET ${target} PROPERTY
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

FIND_PACKAGE(glfw3 REQUIRED)
TARGET_LINK_LIBRARIES(${target} glfw3)

FIND_PACKAGE(OpenGL REQUIRED)
TARGET_LINK_LIBRARIES(${target} OpenGL::GL)

FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)
TARGET_LINK_LIBRARIES(${target} Eigen3::Eigen)

TARGET_INCLUDE_DIRECTORIES(${target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../include")

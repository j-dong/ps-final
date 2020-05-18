SET(target final)

add_executable(${target} main.cpp sim.cpp sync.cpp glad.c imgui.cpp imgui_demo.cpp imgui_draw.cpp imgui_widgets.cpp imgui_impl_glfw.cpp imgui_impl_opengl3.cpp)
TARGET_LINK_LIBRARIES(${target})

set_property(TARGET ${target} PROPERTY CXX_STANDARD 17)
set_property(TARGET ${target} PROPERTY
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

FIND_PACKAGE(glfw3 REQUIRED)
TARGET_LINK_LIBRARIES(${target} glfw)

FIND_PACKAGE(OpenGL REQUIRED)
TARGET_LINK_LIBRARIES(${target} OpenGL::GL)

FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)
TARGET_LINK_LIBRARIES(${target} Eigen3::Eigen)

FIND_PACKAGE(Threads REQUIRED)
TARGET_LINK_LIBRARIES(${target} Threads::Threads)

FIND_PACKAGE(OpenMP)
TARGET_LINK_LIBRARIES(${target} OpenMP::OpenMP_CXX)

TARGET_INCLUDE_DIRECTORIES(${target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../include")

if(CMAKE_DL_LIBS)
  TARGET_LINK_LIBRARIES(${target} ${CMAKE_DL_LIBS})
endif()

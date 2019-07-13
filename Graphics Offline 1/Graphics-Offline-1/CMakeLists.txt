
cmake_minimum_required(VERSION 3.10)
project(Graphics_Offline_1)

set(CMAKE_CXX_STANDARD 11)

include_directories(OpenGL)
include_directories(OpenGL/include) # OpenGL/include has to contain the required OpenGL's .h files
include_directories(OpenGL/lib) # OpenGL/lib has to contain the required OpenGL's .lib files
# glut32.dll must be present in "project-directory/OpenGL/dll/"
include_directories(codes)

add_custom_target(glutdlllib
        COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/OpenGL/dll/glut32.dll ${CMAKE_BINARY_DIR}
        )

# required .lib files to be copied into compiler's proper directory
set(OpenGlLibs glaux glu32 glui32 glut32 opengl32)

add_executable(demo codes/demoMain.cpp)

target_link_libraries(demo ${OpenGlLibs})

add_dependencies(demo glutdlllib)
#######################################################
add_executable(cubeToSphere codes/Sphere-to-Cube.cpp)

target_link_libraries(cubeToSphere ${OpenGlLibs})

add_dependencies(cubeToSphere glutdlllib)

#######################################################
add_executable(wheel codes/Wheel.cpp)

target_link_libraries(wheel ${OpenGlLibs})

add_dependencies(wheel glutdlllib)

#######################################################
add_executable(testing codes/test.cpp)

target_link_libraries(testing ${OpenGlLibs})

add_dependencies(testing glutdlllib)
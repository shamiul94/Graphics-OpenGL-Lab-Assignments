cmake_minimum_required(VERSION 3.10)
project(OpenGLLinuxTest)

set (CMAKE_CXX_STANDARD 14)

set(OpenGlLinkers -lglut -lGLU -lGL)

add_executable(spheretocube codes/sphere-to-cube.cpp)
target_link_libraries(spheretocube ${OpenGlLinkers})

#add_executable(demo-ray-tracer ray-tracer-offline-codes/demo.cpp)
#target_link_libraries(demo-ray-tracer ${OpenGlLinkers})


add_executable(1505038-ray-tracer my_codes/1505038.cpp)
target_link_libraries(1505038-ray-tracer ${OpenGlLinkers})

add_executable(testing-ray ray-tracer-offline-codes/test.cpp)
target_link_libraries(testing-ray ${OpenGlLinkers})

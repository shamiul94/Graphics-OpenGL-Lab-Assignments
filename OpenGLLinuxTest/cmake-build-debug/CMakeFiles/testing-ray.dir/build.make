# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = "/media/shamiul93/Software and Media/Necessary Softwares/Linux Softwares/clion-2018.1.5/bin/cmake/bin/cmake"

# The command to remove a file.
RM = "/media/shamiul93/Software and Media/Necessary Softwares/Linux Softwares/clion-2018.1.5/bin/cmake/bin/cmake" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/shamiul93/CLionProjects/OpenGLLinuxTest

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shamiul93/CLionProjects/OpenGLLinuxTest/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/testing-ray.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/testing-ray.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testing-ray.dir/flags.make

CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o: CMakeFiles/testing-ray.dir/flags.make
CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o: ../ray-tracer-offline-codes/test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shamiul93/CLionProjects/OpenGLLinuxTest/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o -c /home/shamiul93/CLionProjects/OpenGLLinuxTest/ray-tracer-offline-codes/test.cpp

CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shamiul93/CLionProjects/OpenGLLinuxTest/ray-tracer-offline-codes/test.cpp > CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.i

CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shamiul93/CLionProjects/OpenGLLinuxTest/ray-tracer-offline-codes/test.cpp -o CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.s

CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o.requires:

.PHONY : CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o.requires

CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o.provides: CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o.requires
	$(MAKE) -f CMakeFiles/testing-ray.dir/build.make CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o.provides.build
.PHONY : CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o.provides

CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o.provides.build: CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o


# Object files for target testing-ray
testing__ray_OBJECTS = \
"CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o"

# External object files for target testing-ray
testing__ray_EXTERNAL_OBJECTS =

testing-ray: CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o
testing-ray: CMakeFiles/testing-ray.dir/build.make
testing-ray: CMakeFiles/testing-ray.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/shamiul93/CLionProjects/OpenGLLinuxTest/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable testing-ray"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testing-ray.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testing-ray.dir/build: testing-ray

.PHONY : CMakeFiles/testing-ray.dir/build

CMakeFiles/testing-ray.dir/requires: CMakeFiles/testing-ray.dir/ray-tracer-offline-codes/test.cpp.o.requires

.PHONY : CMakeFiles/testing-ray.dir/requires

CMakeFiles/testing-ray.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testing-ray.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testing-ray.dir/clean

CMakeFiles/testing-ray.dir/depend:
	cd /home/shamiul93/CLionProjects/OpenGLLinuxTest/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shamiul93/CLionProjects/OpenGLLinuxTest /home/shamiul93/CLionProjects/OpenGLLinuxTest /home/shamiul93/CLionProjects/OpenGLLinuxTest/cmake-build-debug /home/shamiul93/CLionProjects/OpenGLLinuxTest/cmake-build-debug /home/shamiul93/CLionProjects/OpenGLLinuxTest/cmake-build-debug/CMakeFiles/testing-ray.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testing-ray.dir/depend


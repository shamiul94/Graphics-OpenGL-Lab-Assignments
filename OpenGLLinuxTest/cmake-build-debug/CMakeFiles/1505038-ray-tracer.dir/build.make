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
CMAKE_SOURCE_DIR = /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/1505038-ray-tracer.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/1505038-ray-tracer.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/1505038-ray-tracer.dir/flags.make

CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o: CMakeFiles/1505038-ray-tracer.dir/flags.make
CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o: ../my_codes/1505038.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o -c /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/my_codes/1505038.cpp

CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/my_codes/1505038.cpp > CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.i

CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/my_codes/1505038.cpp -o CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.s

CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o.requires:

.PHONY : CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o.requires

CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o.provides: CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o.requires
	$(MAKE) -f CMakeFiles/1505038-ray-tracer.dir/build.make CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o.provides.build
.PHONY : CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o.provides

CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o.provides.build: CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o


# Object files for target 1505038-ray-tracer
1505038__ray__tracer_OBJECTS = \
"CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o"

# External object files for target 1505038-ray-tracer
1505038__ray__tracer_EXTERNAL_OBJECTS =

1505038-ray-tracer: CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o
1505038-ray-tracer: CMakeFiles/1505038-ray-tracer.dir/build.make
1505038-ray-tracer: CMakeFiles/1505038-ray-tracer.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable 1505038-ray-tracer"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/1505038-ray-tracer.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/1505038-ray-tracer.dir/build: 1505038-ray-tracer

.PHONY : CMakeFiles/1505038-ray-tracer.dir/build

CMakeFiles/1505038-ray-tracer.dir/requires: CMakeFiles/1505038-ray-tracer.dir/my_codes/1505038.cpp.o.requires

.PHONY : CMakeFiles/1505038-ray-tracer.dir/requires

CMakeFiles/1505038-ray-tracer.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/1505038-ray-tracer.dir/cmake_clean.cmake
.PHONY : CMakeFiles/1505038-ray-tracer.dir/clean

CMakeFiles/1505038-ray-tracer.dir/depend:
	cd /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/cmake-build-debug /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/cmake-build-debug /home/shamiul93/Documents/Graphics-all-offline-git/Graphics-OpenGL-Lab-Assignments/OpenGLLinuxTest/cmake-build-debug/CMakeFiles/1505038-ray-tracer.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/1505038-ray-tracer.dir/depend

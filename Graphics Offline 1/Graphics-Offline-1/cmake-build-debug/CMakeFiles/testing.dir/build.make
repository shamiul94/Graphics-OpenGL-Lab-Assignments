# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.10

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2018.1\bin\cmake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2018.1\bin\cmake\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/testing.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/testing.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testing.dir/flags.make

CMakeFiles/testing.dir/codes/test.cpp.obj: CMakeFiles/testing.dir/flags.make
CMakeFiles/testing.dir/codes/test.cpp.obj: CMakeFiles/testing.dir/includes_CXX.rsp
CMakeFiles/testing.dir/codes/test.cpp.obj: ../codes/test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testing.dir/codes/test.cpp.obj"
	C:\PROGRA~2\CODEBL~1\MinGW\bin\G__~1.EXE  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\testing.dir\codes\test.cpp.obj -c "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1\codes\test.cpp"

CMakeFiles/testing.dir/codes/test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testing.dir/codes/test.cpp.i"
	C:\PROGRA~2\CODEBL~1\MinGW\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1\codes\test.cpp" > CMakeFiles\testing.dir\codes\test.cpp.i

CMakeFiles/testing.dir/codes/test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testing.dir/codes/test.cpp.s"
	C:\PROGRA~2\CODEBL~1\MinGW\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1\codes\test.cpp" -o CMakeFiles\testing.dir\codes\test.cpp.s

CMakeFiles/testing.dir/codes/test.cpp.obj.requires:

.PHONY : CMakeFiles/testing.dir/codes/test.cpp.obj.requires

CMakeFiles/testing.dir/codes/test.cpp.obj.provides: CMakeFiles/testing.dir/codes/test.cpp.obj.requires
	$(MAKE) -f CMakeFiles\testing.dir\build.make CMakeFiles/testing.dir/codes/test.cpp.obj.provides.build
.PHONY : CMakeFiles/testing.dir/codes/test.cpp.obj.provides

CMakeFiles/testing.dir/codes/test.cpp.obj.provides.build: CMakeFiles/testing.dir/codes/test.cpp.obj


# Object files for target testing
testing_OBJECTS = \
"CMakeFiles/testing.dir/codes/test.cpp.obj"

# External object files for target testing
testing_EXTERNAL_OBJECTS =

testing.exe: CMakeFiles/testing.dir/codes/test.cpp.obj
testing.exe: CMakeFiles/testing.dir/build.make
testing.exe: CMakeFiles/testing.dir/linklibs.rsp
testing.exe: CMakeFiles/testing.dir/objects1.rsp
testing.exe: CMakeFiles/testing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable testing.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\testing.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testing.dir/build: testing.exe

.PHONY : CMakeFiles/testing.dir/build

CMakeFiles/testing.dir/requires: CMakeFiles/testing.dir/codes/test.cpp.obj.requires

.PHONY : CMakeFiles/testing.dir/requires

CMakeFiles/testing.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\testing.dir\cmake_clean.cmake
.PHONY : CMakeFiles/testing.dir/clean

CMakeFiles/testing.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1" "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1" "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1\cmake-build-debug" "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1\cmake-build-debug" "C:\Users\Heisenberg\Desktop\4-1\graphics sessional\OpenGL transformation, camera rotation, wheel movement offline\Graphics-Offline-1\cmake-build-debug\CMakeFiles\testing.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/testing.dir/depend


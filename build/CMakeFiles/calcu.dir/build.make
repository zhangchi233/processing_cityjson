# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.25.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.25.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/asdfasd/Downloads/calcu

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/asdfasd/Downloads/calcu/build

# Include any dependencies generated for this target.
include CMakeFiles/calcu.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/calcu.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/calcu.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/calcu.dir/flags.make

CMakeFiles/calcu.dir/src/geomtools.cpp.o: CMakeFiles/calcu.dir/flags.make
CMakeFiles/calcu.dir/src/geomtools.cpp.o: /Users/asdfasd/Downloads/calcu/src/geomtools.cpp
CMakeFiles/calcu.dir/src/geomtools.cpp.o: CMakeFiles/calcu.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/asdfasd/Downloads/calcu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/calcu.dir/src/geomtools.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/calcu.dir/src/geomtools.cpp.o -MF CMakeFiles/calcu.dir/src/geomtools.cpp.o.d -o CMakeFiles/calcu.dir/src/geomtools.cpp.o -c /Users/asdfasd/Downloads/calcu/src/geomtools.cpp

CMakeFiles/calcu.dir/src/geomtools.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/calcu.dir/src/geomtools.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/asdfasd/Downloads/calcu/src/geomtools.cpp > CMakeFiles/calcu.dir/src/geomtools.cpp.i

CMakeFiles/calcu.dir/src/geomtools.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/calcu.dir/src/geomtools.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/asdfasd/Downloads/calcu/src/geomtools.cpp -o CMakeFiles/calcu.dir/src/geomtools.cpp.s

CMakeFiles/calcu.dir/src/main.cpp.o: CMakeFiles/calcu.dir/flags.make
CMakeFiles/calcu.dir/src/main.cpp.o: /Users/asdfasd/Downloads/calcu/src/main.cpp
CMakeFiles/calcu.dir/src/main.cpp.o: CMakeFiles/calcu.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/asdfasd/Downloads/calcu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/calcu.dir/src/main.cpp.o"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/calcu.dir/src/main.cpp.o -MF CMakeFiles/calcu.dir/src/main.cpp.o.d -o CMakeFiles/calcu.dir/src/main.cpp.o -c /Users/asdfasd/Downloads/calcu/src/main.cpp

CMakeFiles/calcu.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/calcu.dir/src/main.cpp.i"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/asdfasd/Downloads/calcu/src/main.cpp > CMakeFiles/calcu.dir/src/main.cpp.i

CMakeFiles/calcu.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/calcu.dir/src/main.cpp.s"
	g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/asdfasd/Downloads/calcu/src/main.cpp -o CMakeFiles/calcu.dir/src/main.cpp.s

# Object files for target calcu
calcu_OBJECTS = \
"CMakeFiles/calcu.dir/src/geomtools.cpp.o" \
"CMakeFiles/calcu.dir/src/main.cpp.o"

# External object files for target calcu
calcu_EXTERNAL_OBJECTS =

calcu: CMakeFiles/calcu.dir/src/geomtools.cpp.o
calcu: CMakeFiles/calcu.dir/src/main.cpp.o
calcu: CMakeFiles/calcu.dir/build.make
calcu: /usr/local/lib/libgmpxx.dylib
calcu: /usr/local/lib/libmpfr.dylib
calcu: /usr/local/lib/libgmp.dylib
calcu: CMakeFiles/calcu.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/asdfasd/Downloads/calcu/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable calcu"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/calcu.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/calcu.dir/build: calcu
.PHONY : CMakeFiles/calcu.dir/build

CMakeFiles/calcu.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/calcu.dir/cmake_clean.cmake
.PHONY : CMakeFiles/calcu.dir/clean

CMakeFiles/calcu.dir/depend:
	cd /Users/asdfasd/Downloads/calcu/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/asdfasd/Downloads/calcu /Users/asdfasd/Downloads/calcu /Users/asdfasd/Downloads/calcu/build /Users/asdfasd/Downloads/calcu/build /Users/asdfasd/Downloads/calcu/build/CMakeFiles/calcu.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/calcu.dir/depend

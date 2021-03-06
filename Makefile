# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /usr/local/lib/python3.8/dist-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /usr/local/lib/python3.8/dist-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/dhairya/Documents/Complex Plotter"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/dhairya/Documents/Complex Plotter"

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/lib/python3.8/dist-packages/cmake/data/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/local/lib/python3.8/dist-packages/cmake/data/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start "/home/dhairya/Documents/Complex Plotter/CMakeFiles" "/home/dhairya/Documents/Complex Plotter/CMakeFiles/progress.marks"
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start "/home/dhairya/Documents/Complex Plotter/CMakeFiles" 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named complex_plotter

# Build rule for target.
complex_plotter: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 complex_plotter
.PHONY : complex_plotter

# fast build rule for target.
complex_plotter/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/build
.PHONY : complex_plotter/fast

src/complexPlotter.o: src/complexPlotter.cpp.o

.PHONY : src/complexPlotter.o

# target to build an object file
src/complexPlotter.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/src/complexPlotter.cpp.o
.PHONY : src/complexPlotter.cpp.o

src/complexPlotter.i: src/complexPlotter.cpp.i

.PHONY : src/complexPlotter.i

# target to preprocess a source file
src/complexPlotter.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/src/complexPlotter.cpp.i
.PHONY : src/complexPlotter.cpp.i

src/complexPlotter.s: src/complexPlotter.cpp.s

.PHONY : src/complexPlotter.s

# target to generate assembly for a file
src/complexPlotter.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/src/complexPlotter.cpp.s
.PHONY : src/complexPlotter.cpp.s

src/interpretFunction.o: src/interpretFunction.cpp.o

.PHONY : src/interpretFunction.o

# target to build an object file
src/interpretFunction.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/src/interpretFunction.cpp.o
.PHONY : src/interpretFunction.cpp.o

src/interpretFunction.i: src/interpretFunction.cpp.i

.PHONY : src/interpretFunction.i

# target to preprocess a source file
src/interpretFunction.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/src/interpretFunction.cpp.i
.PHONY : src/interpretFunction.cpp.i

src/interpretFunction.s: src/interpretFunction.cpp.s

.PHONY : src/interpretFunction.s

# target to generate assembly for a file
src/interpretFunction.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/src/interpretFunction.cpp.s
.PHONY : src/interpretFunction.cpp.s

src/specialFunctions.o: src/specialFunctions.cpp.o

.PHONY : src/specialFunctions.o

# target to build an object file
src/specialFunctions.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/src/specialFunctions.cpp.o
.PHONY : src/specialFunctions.cpp.o

src/specialFunctions.i: src/specialFunctions.cpp.i

.PHONY : src/specialFunctions.i

# target to preprocess a source file
src/specialFunctions.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/src/specialFunctions.cpp.i
.PHONY : src/specialFunctions.cpp.i

src/specialFunctions.s: src/specialFunctions.cpp.s

.PHONY : src/specialFunctions.s

# target to generate assembly for a file
src/specialFunctions.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/complex_plotter.dir/build.make CMakeFiles/complex_plotter.dir/src/specialFunctions.cpp.s
.PHONY : src/specialFunctions.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... complex_plotter"
	@echo "... src/complexPlotter.o"
	@echo "... src/complexPlotter.i"
	@echo "... src/complexPlotter.s"
	@echo "... src/interpretFunction.o"
	@echo "... src/interpretFunction.i"
	@echo "... src/interpretFunction.s"
	@echo "... src/specialFunctions.o"
	@echo "... src/specialFunctions.i"
	@echo "... src/specialFunctions.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system


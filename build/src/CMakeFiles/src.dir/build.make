# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /scrfs/apps/cmake/cmake-3.13.3/bin/cmake

# The command to remove a file.
RM = /scrfs/apps/cmake/cmake-3.13.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/adenduku/IsingNucleation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/adenduku/IsingNucleation/build

# Include any dependencies generated for this target.
include src/CMakeFiles/src.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/src.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/src.dir/flags.make

src/CMakeFiles/src.dir/simul.c.o: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/simul.c.o: ../src/simul.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adenduku/IsingNucleation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/src.dir/simul.c.o"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/src.dir/simul.c.o   -c /home/adenduku/IsingNucleation/src/simul.c

src/CMakeFiles/src.dir/simul.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/simul.c.i"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/adenduku/IsingNucleation/src/simul.c > CMakeFiles/src.dir/simul.c.i

src/CMakeFiles/src.dir/simul.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/simul.c.s"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/adenduku/IsingNucleation/src/simul.c -o CMakeFiles/src.dir/simul.c.s

src/CMakeFiles/src.dir/monte_carlo.c.o: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/monte_carlo.c.o: ../src/monte_carlo.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adenduku/IsingNucleation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object src/CMakeFiles/src.dir/monte_carlo.c.o"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/src.dir/monte_carlo.c.o   -c /home/adenduku/IsingNucleation/src/monte_carlo.c

src/CMakeFiles/src.dir/monte_carlo.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/monte_carlo.c.i"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/adenduku/IsingNucleation/src/monte_carlo.c > CMakeFiles/src.dir/monte_carlo.c.i

src/CMakeFiles/src.dir/monte_carlo.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/monte_carlo.c.s"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/adenduku/IsingNucleation/src/monte_carlo.c -o CMakeFiles/src.dir/monte_carlo.c.s

src/CMakeFiles/src.dir/ising_lattice.c.o: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/ising_lattice.c.o: ../src/ising_lattice.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adenduku/IsingNucleation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object src/CMakeFiles/src.dir/ising_lattice.c.o"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/src.dir/ising_lattice.c.o   -c /home/adenduku/IsingNucleation/src/ising_lattice.c

src/CMakeFiles/src.dir/ising_lattice.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/ising_lattice.c.i"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/adenduku/IsingNucleation/src/ising_lattice.c > CMakeFiles/src.dir/ising_lattice.c.i

src/CMakeFiles/src.dir/ising_lattice.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/ising_lattice.c.s"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/adenduku/IsingNucleation/src/ising_lattice.c -o CMakeFiles/src.dir/ising_lattice.c.s

src/CMakeFiles/src.dir/util.c.o: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/util.c.o: ../src/util.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adenduku/IsingNucleation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object src/CMakeFiles/src.dir/util.c.o"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/src.dir/util.c.o   -c /home/adenduku/IsingNucleation/src/util.c

src/CMakeFiles/src.dir/util.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/util.c.i"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/adenduku/IsingNucleation/src/util.c > CMakeFiles/src.dir/util.c.i

src/CMakeFiles/src.dir/util.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/util.c.s"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/adenduku/IsingNucleation/src/util.c -o CMakeFiles/src.dir/util.c.s

src/CMakeFiles/src.dir/media.c.o: src/CMakeFiles/src.dir/flags.make
src/CMakeFiles/src.dir/media.c.o: ../src/media.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/adenduku/IsingNucleation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object src/CMakeFiles/src.dir/media.c.o"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/src.dir/media.c.o   -c /home/adenduku/IsingNucleation/src/media.c

src/CMakeFiles/src.dir/media.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/src.dir/media.c.i"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/adenduku/IsingNucleation/src/media.c > CMakeFiles/src.dir/media.c.i

src/CMakeFiles/src.dir/media.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/src.dir/media.c.s"
	cd /home/adenduku/IsingNucleation/build/src && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/adenduku/IsingNucleation/src/media.c -o CMakeFiles/src.dir/media.c.s

# Object files for target src
src_OBJECTS = \
"CMakeFiles/src.dir/simul.c.o" \
"CMakeFiles/src.dir/monte_carlo.c.o" \
"CMakeFiles/src.dir/ising_lattice.c.o" \
"CMakeFiles/src.dir/util.c.o" \
"CMakeFiles/src.dir/media.c.o"

# External object files for target src
src_EXTERNAL_OBJECTS =

src/libsrc.a: src/CMakeFiles/src.dir/simul.c.o
src/libsrc.a: src/CMakeFiles/src.dir/monte_carlo.c.o
src/libsrc.a: src/CMakeFiles/src.dir/ising_lattice.c.o
src/libsrc.a: src/CMakeFiles/src.dir/util.c.o
src/libsrc.a: src/CMakeFiles/src.dir/media.c.o
src/libsrc.a: src/CMakeFiles/src.dir/build.make
src/libsrc.a: src/CMakeFiles/src.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/adenduku/IsingNucleation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking C static library libsrc.a"
	cd /home/adenduku/IsingNucleation/build/src && $(CMAKE_COMMAND) -P CMakeFiles/src.dir/cmake_clean_target.cmake
	cd /home/adenduku/IsingNucleation/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/src.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/src.dir/build: src/libsrc.a

.PHONY : src/CMakeFiles/src.dir/build

src/CMakeFiles/src.dir/clean:
	cd /home/adenduku/IsingNucleation/build/src && $(CMAKE_COMMAND) -P CMakeFiles/src.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/src.dir/clean

src/CMakeFiles/src.dir/depend:
	cd /home/adenduku/IsingNucleation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/adenduku/IsingNucleation /home/adenduku/IsingNucleation/src /home/adenduku/IsingNucleation/build /home/adenduku/IsingNucleation/build/src /home/adenduku/IsingNucleation/build/src/CMakeFiles/src.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/src.dir/depend

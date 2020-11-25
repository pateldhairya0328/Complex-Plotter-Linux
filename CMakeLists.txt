# cmake needs this line
cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 17)

# Define project name
project(complex_plotter_project)
include_directories(include)
file(GLOB SOURCES "src/*.cpp")

# Find OpenCV, you may need to set OpenCV_DIR variable
# to the absolute path to the directory containing OpenCVConfig.cmake file
# via the command line or GUI
find_package(OpenCV REQUIRED)

# If the package has been found, several variables will
# be set, you can find the full list with descriptions
# in the OpenCVConfig.cmake file.
# Print some message showing some of them
message(STATUS "OpenCV library status:")
message(STATUS "    config: ${OpenCV_DIR}")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")

# Declare the executable target built from your sources
add_executable(complex_plotter ${SOURCES})

add_library(complex_bessel SHARED IMPORTED)
set_property(TARGET complex_bessel PROPERTY IMPORTED_LOCATION "/usr/lib/libcomplex_bessel.so")

# Link your application with OpenCV libraries
target_link_libraries(complex_plotter PRIVATE ${OpenCV_LIBS} complex_bessel)

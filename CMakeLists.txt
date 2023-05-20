cmake_minimum_required(VERSION 3.26)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES 
 /usr/include/c++/11
 /usr/include/x86_64-linux-gnu/c++/11
 /usr/include/c++/11/backward
 /usr/lib/gcc/x86_64-linux-gnu/11/include
 /usr/local/include
 /usr/include/x86_64-linux-gnu
 /usr/include
)
set(CMAKE_CXX_FLAGS -Wl,--copy-dt-needed-entries)

project(pathtracer
	VERSION 1.0
	LANGUAGES CXX)

#find_package(Boost REQUIRED components MATH)

file(GLOB_RECURSE deps_SRC CONFIGURE_DEPENDS "dep/*.cpp" "dep/*.hpp" "dep/*.h")
#message(STATUS "files: ${deps_SRC}")
add_library(Dependencies ${deps_SRC})
target_include_directories(Dependencies PUBLIC dep) 

file(GLOB framework_src CONFIGURE_DEPENDS 
	"*.cpp" "*.h"
       	"geometry/*.cpp" "geometry/*.h"
	"acceleration/*.cpp" "acceleration/*.h"
	"bxdf/*.cpp" "bxdf/*.h"
	"integrator/*.cpp" "integrator/*.h"
	"medium/*.cpp" "medium/*.h"
)
add_library(Framework ${framework_src})
target_link_libraries(Framework PUBLIC Dependencies)
target_link_libraries(Framework PUBLIC openvdb)


add_executable(Pathtracer pathtracer.cpp)

#target_link_libraries(Pathtracer PRIVATE glm::glm)
#target_link_libraries(Pathtracer PRIVATE deps_SRC)
target_link_libraries(Pathtracer PRIVATE Dependencies)
target_link_libraries(Pathtracer PRIVATE Framework)
#target_link_libraries(Pathtracer PRIVATE Lib)

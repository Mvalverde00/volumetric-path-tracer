cmake_minimum_required(VERSION 3.26)

project(pathtracer
	VERSION 1.0
	LANGUAGES CXX)

#find_package(glm REQUIRED)

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
)
add_library(Framework ${framework_src})
target_link_libraries(Framework PUBLIC Dependencies)


add_executable(Pathtracer pathtracer.cpp)

#target_link_libraries(Pathtracer PRIVATE glm::glm)
#target_link_libraries(Pathtracer PRIVATE deps_SRC)
target_link_libraries(Pathtracer PRIVATE Dependencies)
target_link_libraries(Pathtracer PRIVATE Framework)
#target_link_libraries(Pathtracer PRIVATE Lib)

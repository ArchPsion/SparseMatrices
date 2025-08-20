cmake_minimum_required(VERSION 3.16)

project(foo VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra -Warith-conversion -pedantic -Wpedantic -g -ggdb")

find_package(Qt6 REQUIRED COMPONENTS Widgets)

qt_standard_project_setup()

qt_add_executable(	foo
			
			HexSparseMatrix.hpp
			QSparseMatrixWindow.hpp
			
			Main.cpp
)

target_link_libraries(foo PRIVATE Qt6::Widgets)

set_target_properties(		foo
				PROPERTIES
				WIN32_EXECUTABLE ON
    				MACOSX_BUNDLE ON
)

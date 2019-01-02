cmake_minimum_required(VERSION 3.1)

project(FunctionPlot)

include(FindPkgConfig)

pkg_search_module(SDL2 REQUIRED sdl2)
find_package(Magnum REQUIRED GL MeshTools Primitives Shaders Sdl2Application)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_executable(function-plot FunctionPlot.cpp)
target_include_directories(function-plot PRIVATE ${SDL2_INCLUDE_DIRS})
target_link_libraries(function-plot PRIVATE
    Magnum::Application
    Magnum::GL
    Magnum::Magnum
    Magnum::MeshTools
    Magnum::Primitives
    Magnum::Shaders)

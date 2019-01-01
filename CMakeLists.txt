cmake_minimum_required(VERSION 3.1)

project(MagnumPrimitivesExample)

include(FindPkgConfig)

pkg_search_module(SDL2 REQUIRED sdl2)
find_package(Magnum REQUIRED GL MeshTools Primitives Shaders Sdl2Application)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_executable(magnum-primitives PrimitivesExample.cpp)
target_include_directories(magnum-primitives PRIVATE ${SDL2_INCLUDE_DIRS})
target_link_libraries(magnum-primitives PRIVATE
    Magnum::Application
    Magnum::GL
    Magnum::Magnum
    Magnum::MeshTools
    Magnum::Primitives
    Magnum::Shaders)

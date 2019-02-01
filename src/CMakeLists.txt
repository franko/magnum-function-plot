find_package(Magnum REQUIRED GL MeshTools Shaders Sdl2Application)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

add_executable(function-plot Units.cpp FunctionPlot.cpp)
target_link_libraries(function-plot PRIVATE
    Magnum::Application
    Magnum::GL
    Magnum::Magnum
    Magnum::MeshTools
    Magnum::Shaders)

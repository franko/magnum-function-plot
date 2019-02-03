find_package(Magnum REQUIRED GL MeshTools Shaders Text Sdl2Application)

set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

corrade_add_resource(FunctionPlot_RESOURCES resources.conf)

add_executable(function-plot Units.cpp FunctionPlot.cpp ${FunctionPlot_RESOURCES})
target_link_libraries(function-plot PRIVATE
    Magnum::Application
    Magnum::GL
    Magnum::Magnum
    Magnum::MeshTools
    Magnum::Shaders
    Magnum::Text)

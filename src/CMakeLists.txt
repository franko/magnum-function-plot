find_package(PkgConfig)
pkg_check_modules(FOX REQUIRED fox)

find_package(Magnum REQUIRED GL MeshTools Shaders Text Sdl2Application)

# With the pedantic flags above we have a lot of warnings for FOX header files
# using old style casts and similar.
# set_directory_properties(PROPERTIES CORRADE_USE_PEDANTIC_FLAGS ON)

corrade_add_resource(FunctionPlot_RESOURCES resources.conf)

add_executable(function-plot Units.cpp FunctionPlot.cpp ${FunctionPlot_RESOURCES})
target_link_libraries(function-plot PRIVATE
    Magnum::Application
    Magnum::GL
    Magnum::Magnum
    Magnum::MeshTools
    Magnum::Shaders
    Magnum::Text)

link_directories(${FOX_LIBRARY_DIRS})

# The flags below need to be set manually. The WIN32 flag may be not
# required.
# Normally it would be better to have these flags provided by pkg-config
# or a cmake find package function.
add_definitions(-DWIN32=1 -DHAVE_GLU_H=1 -DHAVE_GL_H=1)

set(CMAKE_EXE_LINKER_FLAGS ${FOX_LDFLAGS} ${FOX_LDFLAGS_OTHER})

add_executable(gltest gltest.cpp)
add_executable(fox-glplot fox-glplot.cpp)

# Does not seems to actually work.
# set_target_properties(gltest PROPERTIES CXX_STANDARD 98)
# target_compile_features(gltest PUBLIC cxx_std_98)

target_link_libraries(gltest ${FOX_STATIC_LIBRARIES})
target_include_directories(gltest PUBLIC ${FOX_INCLUDE_DIRS})

target_link_libraries(fox-glplot ${FOX_STATIC_LIBRARIES})
target_include_directories(fox-glplot PUBLIC ${FOX_INCLUDE_DIRS})

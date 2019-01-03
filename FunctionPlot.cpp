#include <cmath>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Trade/MeshData3D.h>

using namespace Magnum;
using namespace Magnum::Math::Literals;

template <typename F>
Trade::MeshData3D mathFunctionMeshDataTriangles(float x_min, float x_max, float y_min, float y_max, F evalF) {
    const int nx = 8, ny = 8;
    std::vector<Vector3> positions{};
    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i <= nx; i++) {
            const float x = x_min + i * (x_max - x_min) / nx;
            const float y = y_min + j * (y_max - y_min) / ny;
            const float z = evalF(x, y);
            positions.push_back(Vector3{x, y, z});
        }
    }

    auto getPositionAtIndex = [&](int i, int j) { return positions[j * (nx + 1) + i]; };

    std::vector<Vector3> positionsDirect{};
    for (int j = 0; j <= ny - 1; j++) {
        for (int i = 1; i <= nx; i++) {
            positionsDirect.push_back(getPositionAtIndex(i - 1, j));
            positionsDirect.push_back(getPositionAtIndex(i - 1, j + 1));
            positionsDirect.push_back(getPositionAtIndex(i, j));

            positionsDirect.push_back(getPositionAtIndex(i - 1, j + 1));
            positionsDirect.push_back(getPositionAtIndex(i, j + 1));
            positionsDirect.push_back(getPositionAtIndex(i, j));
        }
    }

    return Trade::MeshData3D{MeshPrimitive::Triangles,
        {}, {positionsDirect}, {}, {}, {}, nullptr};
}

class MyApp: public Platform::Application {
    public:
        explicit MyApp(const Arguments& arguments);

    private:
        void drawEvent() override;
        void mousePressEvent(MouseEvent& event) override;
        void mouseReleaseEvent(MouseEvent& event) override;
        void mouseMoveEvent(MouseMoveEvent& event) override;

        GL::Buffer _vertexBuffer;
        GL::Mesh _mesh;
        Shaders::MeshVisualizer _shader;

        Matrix4 _transformation, _projection;
        Vector2i _previousMousePosition;
        Color3 _color;
};

MyApp::MyApp(const Arguments& arguments):
    Platform::Application{arguments, Configuration{}.setTitle("Function Plot App"), GLConfiguration{}.setSampleCount(8)},
    _shader{Shaders::MeshVisualizer::Flag::Wireframe|Shaders::MeshVisualizer::Flag::NoGeometryShader}
{
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    // Multisampling is enabled by default.
    // GL::Renderer::enable(GL::Renderer::Feature::Multisampling);

    const float gauss_c = 5.0f;
    auto f = [=](float x, float y) {
        return std::exp(-gauss_c * (x*x + y*y));
    };
#if 0
    auto df = [=](float x, float y) {
        const float e = std::exp(-gauss_c * (x*x + y*y));
        return Vector2{-2*gauss_c*x*e, -2*gauss_c*y*e};
    };
#endif
    const Trade::MeshData3D functionMeshData = mathFunctionMeshDataTriangles(-1.0, 1.0, -1.0, 1.0, f);

    _vertexBuffer.setData(functionMeshData.positions(0));

    _mesh.setPrimitive(functionMeshData.primitive())
        .setCount(functionMeshData.positions(0).size())
        .addVertexBuffer(_vertexBuffer, 0, Shaders::MeshVisualizer::Position{});

    _transformation = Matrix4::rotationX(-60.0_degf);
    _projection =
        Matrix4::perspectiveProjection(
            35.0_degf, Vector2{windowSize()}.aspectRatio(), 0.01f, 100.0f)*
        Matrix4::translation(Vector3::zAxis(-8.0f));
    _color = Color3::fromHsv(195.0_degf, 1.0f, 1.0f);
}

void MyApp::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Color|GL::FramebufferClear::Depth);

    _shader
        .setColor(0x2f83cc_rgbf)
        .setWireframeColor(0xdcdcdc_rgbf)
        .setTransformationProjectionMatrix(_projection * _transformation);

    _mesh.draw(_shader);

    swapBuffers();
}

void MyApp::mousePressEvent(MouseEvent& event) {
    if(event.button() != MouseEvent::Button::Left) return;

    _previousMousePosition = event.position();
    event.setAccepted();
}

void MyApp::mouseReleaseEvent(MouseEvent& event) {
    event.setAccepted();
}

void MyApp::mouseMoveEvent(MouseMoveEvent& event) {
    if(!(event.buttons() & MouseMoveEvent::Button::Left)) return;

    const Vector2 delta = 3.0f*
        Vector2{event.position() - _previousMousePosition}/
        Vector2{GL::defaultFramebuffer.viewport().size()};

    _transformation =
        Matrix4::rotationX(Rad{delta.y()})*
        _transformation*
        Matrix4::rotationY(Rad{delta.x()});

    _previousMousePosition = event.position();
    event.setAccepted();
    redraw();
}

MAGNUM_APPLICATION_MAIN(MyApp)

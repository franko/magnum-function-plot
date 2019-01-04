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

template <typename F, typename dF>
Trade::MeshData3D mathFunctionMeshData(float x_min, float x_max, float y_min, float y_max, F evalF, dF evaldF) {
    const int nx = 16, ny = 16;
    std::vector<Vector3> positions{};
    std::vector<Vector3> normals{};
    for (int j = 0; j <= ny; j++) {
        for (int i = 0; i <= nx; i++) {
            const float x = x_min + i * (x_max - x_min) / nx;
            const float y = y_min + j * (y_max - y_min) / ny;
            const float z = evalF(x, y);
            const Vector2 df = evaldF(x, y);
            const float nf = std::sqrt(df.x() * df.x() + df.y() * df.y() + 1.0f);
            positions.push_back(Vector3{x, y, z});
            normals.push_back(Vector3{-df.x() / nf, -df.y() / nf, 1 / nf});
        }
    }

    auto index = [=](int i, int j) { return j * (nx + 1) + i; };

    std::vector<UnsignedInt> indices{};
    for (int j = 0; j <= ny - 1; j++) {
        for (int i = 1; i <= nx; i++) {
            indices.push_back(index(i - 1, j));
            indices.push_back(index(i - 1, j + 1));
            indices.push_back(index(i, j));

            indices.push_back(index(i - 1, j + 1));
            indices.push_back(index(i, j + 1));
            indices.push_back(index(i, j));
        }
    }

    return Trade::MeshData3D{MeshPrimitive::Triangles,
        indices, {positions}, {normals}, {}, {}, nullptr};
}

class MyApp: public Platform::Application {
    public:
        explicit MyApp(const Arguments& arguments);

    private:
        void drawEvent() override;
        void mousePressEvent(MouseEvent& event) override;
        void mouseReleaseEvent(MouseEvent& event) override;
        void mouseMoveEvent(MouseMoveEvent& event) override;

        GL::Buffer _indexBuffer, _vertexBuffer;
        GL::Mesh _mesh;
        Shaders::MeshVisualizer _shader;

        Matrix4 _transformation, _projection;
        Vector2i _previousMousePosition;
        Color3 _color;
};

MyApp::MyApp(const Arguments& arguments):
    Platform::Application{arguments, Configuration{}.setTitle("Function Plot App"), GLConfiguration{}.setSampleCount(8)},
    _shader{Shaders::MeshVisualizer::Flag::Wireframe}
{
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    // Multisampling is enabled by default.
    // GL::Renderer::enable(GL::Renderer::Feature::Multisampling);

    const float gauss_c = 5.0f;
    auto f = [=](float x, float y) {
        return std::exp(-gauss_c * (x*x + y*y));
    };
    auto df = [=](float x, float y) {
        const float e = std::exp(-gauss_c * (x*x + y*y));
        return Vector2{-2*gauss_c*x*e, -2*gauss_c*y*e};
    };
    const Trade::MeshData3D functionMeshData = mathFunctionMeshData(-1.0, 1.0, -1.0, 1.0, f, df);

    _vertexBuffer.setData(functionMeshData.positions(0));
    _indexBuffer.setData(functionMeshData.indices());

    _mesh.setPrimitive(functionMeshData.primitive())
        .setCount(functionMeshData.indices().size())
        .addVertexBuffer(_vertexBuffer, 0, Shaders::MeshVisualizer::Position{})
        .setIndexBuffer(_indexBuffer, 0, GL::MeshIndexType::UnsignedInt);

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
        .setViewportSize(Vector2{framebufferSize()})
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

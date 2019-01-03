#include <cmath>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Trade/MeshData3D.h>

using namespace Magnum;
using namespace Magnum::Math::Literals;

template <typename F, typename dF>
Trade::MeshData3D mathFunctionMeshData(float x_min, float x_max, float y_min, float y_max, F evalF, dF evaldF) {
    const int i_min = -10, i_max = 10;
    const int j_min = -10, j_max = 10;
    std::vector<Vector3> positions{};
    std::vector<Vector3> normals{};
    for (int j = j_min; j <= j_max; j++) {
        for (int i = i_min; i <= i_max; i++) {
            const float x = x_min + (i - i_min) * (x_max - x_min) / (i_max - i_min);
            const float y = y_min + (j - j_min) * (y_max - y_min) / (j_max - j_min);
            const float z = evalF(x, y);
            const Vector2 df = evaldF(x, y);
            const float nf = std::sqrt(df.x() * df.x() + df.y() * df.y() + 1.0f);
            positions.push_back(Vector3{x, y, z});
            normals.push_back(Vector3{-df.x() / nf, -df.y() / nf, 1 / nf});
        }
    }

    auto index = [&](int i, int j) { return (j - j_min) * (i_max - i_min + 1) + (i - i_min); };

    std::vector<UnsignedInt> indices{};
    for (int j = j_min; j <= j_max - 1; j++) {
        for (int i = i_min + 1; i <= i_max; i++) {
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
        Shaders::Phong _shader;

        Matrix4 _transformation, _projection;
        Vector2i _previousMousePosition;
        Color3 _color;
};

MyApp::MyApp(const Arguments& arguments):
    Platform::Application{arguments, Configuration{}.setTitle("Function Plot App"), GLConfiguration{}.setSampleCount(16)}
{
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    // Multisampling is enabled by default.
    // GL::Renderer::enable(GL::Renderer::Feature::Multisampling);

    const float gauss_c = 5.0f;
    auto f = [&](float x, float y) {
        return std::exp(-gauss_c * (x*x + y*y));
    };
    auto df = [=](float x, float y) {
        const float e = std::exp(-gauss_c * (x*x + y*y));
        return Vector2{-2*gauss_c*x*e, -2*gauss_c*y*e};
    };
    const Trade::MeshData3D functionMeshData = mathFunctionMeshData(-1.0, 1.0, -1.0, 1.0, f, df);

    _vertexBuffer.setData(MeshTools::interleave(functionMeshData.positions(0), functionMeshData.normals(0)));

    Containers::Array<char> indexData;
    MeshIndexType indexType;
    UnsignedInt indexStart, indexEnd;
    std::tie(indexData, indexType, indexStart, indexEnd) =
        MeshTools::compressIndices(functionMeshData.indices());
    _indexBuffer.setData(indexData);

    _mesh.setPrimitive(functionMeshData.primitive())
        .setCount(functionMeshData.indices().size())
        .addVertexBuffer(_vertexBuffer, 0, Shaders::Phong::Position{},
                                           Shaders::Phong::Normal{})
        .setIndexBuffer(_indexBuffer, 0, indexType, indexStart, indexEnd);

    _transformation = Matrix4::rotationX(-60.0_degf);
    _projection =
        Matrix4::perspectiveProjection(
            35.0_degf, Vector2{windowSize()}.aspectRatio(), 0.01f, 100.0f)*
        Matrix4::translation(Vector3::zAxis(-8.0f));
    _color = Color3::fromHsv(195.0_degf, 1.0f, 1.0f);
}

void MyApp::drawEvent() {
    GL::defaultFramebuffer.clear(
        GL::FramebufferClear::Color|GL::FramebufferClear::Depth);

    _shader.setLightPosition({7.0f, 5.0f, 2.5f})
        .setLightColor(Color3{1.0f})
        .setDiffuseColor(_color)
        .setAmbientColor(Color3::fromHsv(_color.hue(), 1.0f, 0.3f))
        .setTransformationMatrix(_transformation)
        .setNormalMatrix(_transformation.rotationScaling())
        .setProjectionMatrix(_projection);
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

#include <cmath>
#include <utility>
#include <Magnum/GL/Buffer.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/Trade/MeshData3D.h>
#include "CartesianGrid.h"

using namespace Magnum;
using namespace Magnum::Math::Literals;

template <typename F, typename Grid>
std::vector<Vector3> mathFunctionPositions(const Grid& grid, F evalF) {
    std::vector<Vector3> positions{};
    const int nPoints = grid.pointsCount();
    positions.reserve(nPoints);
    auto evalPointF = [&](const Vector2 p) {
        const float z = evalF(p.x(), p.y());
        positions.push_back(Vector3{p.x(), p.y(), z});
    };
    grid.spanGridPoints(evalPointF);
    return positions;
}

template <typename F, typename dF, typename Grid>
std::pair<std::vector<Vector3>, std::vector<Vector3>> mathFunctionPositionsNormals(const Grid& grid, F evalF, dF evaldF) {
    std::vector<Vector3> positions{};
    std::vector<Vector3> normals{};
    const int nPoints = grid.pointsCount();
    positions.reserve(nPoints);
    normals.reserve(nPoints);
    auto evalPointF = [&](const Vector2 p) {
        const float z = evalF(p.x(), p.y());
        positions.push_back(Vector3{p.x(), p.y(), z});
        const Vector2 df = evaldF(p.x(), p.y());
        const float nf = std::sqrt(df.x() * df.x() + df.y() * df.y() + 1.0f);
        normals.push_back(Vector3{-df.x() / nf, -df.y() / nf, 1 / nf});
    };
    grid.spanGridPoints(evalPointF);
    return std::make_pair(positions, normals);
}

template <typename F, typename dF, typename Grid>
Trade::MeshData3D mathFunctionMeshData(const Grid& grid, F evalF, dF evaldF) {
    std::vector<Vector3> positions{}, normals{};
    std::tie(positions, normals) = mathFunctionPositionsNormals(grid, evalF, evaldF);
    std::vector<UnsignedInt> indices{};
    indices.reserve(grid.trianglesCount());
    auto evalTriangleF = [&](int indexA, int indexB, int indexC) {
        indices.push_back(indexA);
        indices.push_back(indexB);
        indices.push_back(indexC);
    };
    grid.spanTrianglesIndices(evalTriangleF);
    return Trade::MeshData3D{MeshPrimitive::Triangles,
        indices, {positions}, {normals}, {}, {}, nullptr};
}

template <typename F, typename Grid>
Trade::MeshData3D mathFunctionLinesData(const Grid& grid, F evalF) {
    std::vector<Vector3> positions = mathFunctionPositions(grid, evalF);
    std::vector<UnsignedInt> indices{};
    indices.reserve(grid.linesCount());
    auto evalLinesF = [&](int indexA, int indexB) {
        indices.push_back(indexA);
        indices.push_back(indexB);
    };
    grid.spanLinesIndices(evalLinesF);
    return Trade::MeshData3D{MeshPrimitive::Lines,
        indices, {positions}, {}, {}, {}, nullptr};
}

class MyApp: public Platform::Application {
    public:
        explicit MyApp(const Arguments& arguments);

    private:
        void drawEvent() override;
        void mousePressEvent(MouseEvent& event) override;
        void mouseReleaseEvent(MouseEvent& event) override;
        void mouseMoveEvent(MouseMoveEvent& event) override;

        GL::Buffer _indexBuffer, _indexLinesBuffer, _vertexBuffer, _vertexLinesBuffer;
        GL::Mesh _mesh, _lines;
        Shaders::Phong _shader;
        Shaders::Flat3D _flatShader;

        Matrix4 _transformation, _projection;
        Vector2i _previousMousePosition;
        Color3 _color;
};

MyApp::MyApp(const Arguments& arguments):
    Platform::Application{arguments, Configuration{}.setTitle("Function Plot App"), GLConfiguration{}.setSampleCount(16)}
{
    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::PolygonOffsetFill);
    GL::Renderer::setPolygonOffset(1.0f, 1.0f);
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
    // To obtain a cartesianGrid use CartesianGrid<NoTransform, 4>
    const CartesianGrid<CircleTransform, 4> grid{Vector2{-1.0f, -1.0f}, Vector2{1.0f, 1.0f}, 40, 40};
    const Trade::MeshData3D functionMeshData = mathFunctionMeshData(grid, f, df);
    const Trade::MeshData3D functionLines = mathFunctionLinesData(grid, f);

    _vertexBuffer.setData(MeshTools::interleave(functionMeshData.positions(0), functionMeshData.normals(0)));

    Containers::Array<char> indexData;
    MeshIndexType indexType;
    UnsignedInt indexStart, indexEnd;
    std::tie(indexData, indexType, indexStart, indexEnd) =
        MeshTools::compressIndices(functionMeshData.indices());
    _indexBuffer.setData(indexData);

    _vertexLinesBuffer.setData(functionLines.positions(0));
    _indexLinesBuffer.setData(functionLines.indices());

    _mesh.setPrimitive(functionMeshData.primitive())
        .setCount(functionMeshData.indices().size())
        .addVertexBuffer(_vertexBuffer, 0, Shaders::Phong::Position{},
                                           Shaders::Phong::Normal{})
        .setIndexBuffer(_indexBuffer, 0, indexType, indexStart, indexEnd);

    _lines.setPrimitive(functionLines.primitive())
        .setCount(functionLines.indices().size())
        .addVertexBuffer(_vertexLinesBuffer, 0, Shaders::Flat3D::Position{})
        .setIndexBuffer(_indexLinesBuffer, 0, MeshIndexType::UnsignedInt, 0, functionLines.indices().size());

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

    _flatShader.setColor(Color3{0.9f})
        .setTransformationProjectionMatrix(_projection*_transformation);
    _lines.draw(_flatShader);

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

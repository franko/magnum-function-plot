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

static Trade::MeshData3D gridData(Vector2 llc, Vector2 urc, float stepX, float stepY, float z) {
    std::vector<UnsignedInt> indices{};
    std::vector<Vector3> positions;
    const float x0 = llc.x(), y0 = llc.y();
    const float x1 = urc.x(), y1 = urc.y();
    positions.push_back(Vector3{x0, y0, z});
    positions.push_back(Vector3{x0, y1, z});
    positions.push_back(Vector3{x1, y1, z});
    positions.push_back(Vector3{x1, y0, z});
    for (UnsignedInt i : {0, 1, 1, 2, 2, 3, 3, 0}) {
        indices.push_back(i);
    }
    int currentIndex = 4;
    const float xEpsilon = stepX / 40, yEpsilon = stepY / 40;
    const int i0 = ceil((x0 + xEpsilon) / stepX), i1 = floor((x1 - xEpsilon) / stepX);
    for (int i = i0; i <= i1; i++) {
        positions.push_back(Vector3{i * stepX, y0, z});
        positions.push_back(Vector3{i * stepX, y1, z});
        indices.push_back(currentIndex++);
        indices.push_back(currentIndex++);
    }
    const int j0 = ceil((y0 + yEpsilon) / stepY), j1 = floor((y1 - yEpsilon) / stepY);
    for (int j = j0; j <= j1; j++) {
        positions.push_back(Vector3{x0, j * stepY, z});
        positions.push_back(Vector3{x1, j * stepY, z});
        indices.push_back(currentIndex++);
        indices.push_back(currentIndex++);
    }
    return Trade::MeshData3D{MeshPrimitive::Lines,
        indices, {positions}, {}, {}, {}, nullptr};
}

class MyApp: public Platform::Application {
    public:
        explicit MyApp(const Arguments& arguments);

        Matrix4 transformation() const {
            return _rotation * _model;
        }

    private:
        void drawEvent() override;
        void mousePressEvent(MouseEvent& event) override;
        void mouseReleaseEvent(MouseEvent& event) override;
        void mouseMoveEvent(MouseMoveEvent& event) override;

        GL::Buffer _indexBuffer, _indexLinesBuffer, _vertexBuffer, _vertexLinesBuffer;
        GL::Buffer _indexGridBuffer, _vertexGridBuffer;
        GL::Mesh _mesh, _lines;
        GL::Mesh _baseGrid;
        Shaders::Phong _shader;
        Shaders::Flat3D _flatShader;

        Matrix4 _rotation, _model, _projection;
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

    const float gauss_c = 0.05f;
    auto f = [=](float x, float y) {
        const float xc = (x - 20.0f), yc = (y - 20.0f);
        return std::exp(-gauss_c * (xc*xc + yc*yc));
    };
    auto df = [=](float x, float y) {
        const float xc = (x - 20.0f), yc = (y - 20.0f);
        const float e = std::exp(-gauss_c * (xc*xc + yc*yc));
        return Vector2{-2*gauss_c*xc*e, -2*gauss_c*yc*e};
    };
    // To obtain a cartesianGrid use CartesianGrid<NoTransform, 4>
    const CartesianGrid<NoTransform, 4> grid{Vector2{10.0f, 10.0f}, Vector2{30.0f, 30.0f}, 40, 40};
    const Trade::MeshData3D functionMeshData = mathFunctionMeshData(grid, f, df);
    const Trade::MeshData3D functionLines = mathFunctionLinesData(grid, f);
    const Trade::MeshData3D bottomGrid = gridData(Vector2{10.0f, 10.0f}, Vector2{30.0f, 30.0f}, 3.0f, 3.0f, -0.5f);

    _vertexBuffer.setData(MeshTools::interleave(functionMeshData.positions(0), functionMeshData.normals(0)));

    Containers::Array<char> indexData;
    MeshIndexType indexType;
    UnsignedInt indexStart, indexEnd;
    std::tie(indexData, indexType, indexStart, indexEnd) =
        MeshTools::compressIndices(functionMeshData.indices());
    _indexBuffer.setData(indexData);

    _vertexLinesBuffer.setData(functionLines.positions(0));
    _indexLinesBuffer.setData(functionLines.indices());

    _vertexGridBuffer.setData(bottomGrid.positions(0));
    _indexGridBuffer.setData(bottomGrid.indices());

    _mesh.setPrimitive(functionMeshData.primitive())
        .setCount(functionMeshData.indices().size())
        .addVertexBuffer(_vertexBuffer, 0, Shaders::Phong::Position{},
                                           Shaders::Phong::Normal{})
        .setIndexBuffer(_indexBuffer, 0, indexType, indexStart, indexEnd);

    _lines.setPrimitive(functionLines.primitive())
        .setCount(functionLines.indices().size())
        .addVertexBuffer(_vertexLinesBuffer, 0, Shaders::Flat3D::Position{})
        .setIndexBuffer(_indexLinesBuffer, 0, MeshIndexType::UnsignedInt, 0, functionLines.indices().size());

    _baseGrid.setPrimitive(bottomGrid.primitive())
        .setCount(bottomGrid.indices().size())
        .addVertexBuffer(_vertexGridBuffer, 0, Shaders::Flat3D::Position{})
        .setIndexBuffer(_indexGridBuffer, 0, MeshIndexType::UnsignedInt, 0, bottomGrid.indices().size());

    _model = Matrix4::scaling({0.1f, 0.1f, 1.0f}) * Matrix4::translation({-20.0f, -20.0f, 0.0f});

    _rotation = Matrix4::rotationX(-60.0_degf);
    _projection =
        Matrix4::perspectiveProjection(
            35.0_degf, Vector2{windowSize()}.aspectRatio(), 0.01f, 100.0f)*
        Matrix4::translation(Vector3::zAxis(-8.0f));
    _color = 0x2f83cc_rgbf;
}

void MyApp::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Depth).clearColor(Color4{1.0f});

    _shader.setLightPosition({7.0f, 5.0f, 2.5f})
        .setLightColor(Color3{1.0f})
        .setDiffuseColor(_color)
        .setAmbientColor(Color3::fromHsv(_color.hue(), 1.0f, 0.6f))
        .setTransformationMatrix(transformation())
        .setNormalMatrix(transformation().rotationScaling().inverted().transposed())
        .setProjectionMatrix(_projection);
    _mesh.draw(_shader);

    _flatShader.setColor(0xdcdcdc_rgbf)
        .setTransformationProjectionMatrix(_projection*transformation());
    _lines.draw(_flatShader);

    _flatShader.setColor(Color4{0.7f});
    _baseGrid.draw(_flatShader);

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

    _rotation =
        Matrix4::rotationX(Rad{delta.y()})*
        _rotation*
        Matrix4::rotationY(Rad{delta.x()});

    _previousMousePosition = event.position();
    event.setAccepted();
    redraw();
}

MAGNUM_APPLICATION_MAIN(MyApp)

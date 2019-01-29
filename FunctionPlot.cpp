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
#include <Magnum/Math/Algorithms/GaussJordan.h>
#include "CartesianGrid.h"
#include "Units.h"

using namespace Magnum;
using namespace Magnum::Math::Literals;

template <typename Grid, typename F>
std::pair<float, float> mathFunctionFindExtrema(const Grid& grid, F f) {
    float zMin = 1, zMax = -1;
    grid.spanGridPoints([&](const Vector2& p) {
        float z = f(p.x(), p.y());
        if (zMin > zMax) {
            zMin = z;
            zMax = z;
        } else {
            if (z < zMin) {
                zMin = z;
            }
            if (z > zMax) {
                zMax = z;
            }
        }
    });
    return {zMin, zMax};
}

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
        // The normal vector does not need to be normalized.
        normals.push_back(Vector3{-df.x(), -df.y(), 1});
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

template <typename Grid>
Trade::MeshData3D unitsGridData(const Grid& grid, const Units& unitsX, const Units& unitsY, const float z) {
    std::vector<UnsignedInt> indices{};
    std::vector<Vector3> positions;
    const float x0 = grid.llc().x(), y0 = grid.llc().y();
    const float x1 = grid.urc().x(), y1 = grid.urc().y();
    positions.push_back(Vector3{x0, y0, z});
    positions.push_back(Vector3{x0, y1, z});
    positions.push_back(Vector3{x1, y1, z});
    positions.push_back(Vector3{x1, y0, z});
    for (UnsignedInt i : {0, 1, 1, 2, 2, 3, 3, 0}) {
        indices.push_back(i);
    }
    int currentIndex = 4;

    UnitsIterator xIter{unitsX}, yIter{unitsY};
    double xval;
    const char *xlabel;
    while (xIter.next(xval, xlabel)) {
        const float xvalf = float(xval);
        if (xvalf < x0 || xvalf > x1) continue;
        positions.push_back(Vector3{xvalf, y0, z});
        positions.push_back(Vector3{xvalf, y1, z});
        indices.push_back(currentIndex++);
        indices.push_back(currentIndex++);
    }

    double yval;
    const char *ylabel;
    while (yIter.next(yval, ylabel)) {
        const float yvalf = float(yval);
        if (yvalf < y0 || yvalf > y1) continue;
        positions.push_back(Vector3{x0, yvalf, z});
        positions.push_back(Vector3{x1, yvalf, z});
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

        Matrix3 normalsTransformation() const {
            // For the normals we need the inverse of the transpose of the matrix, excluding the
            // translation term.
            // As the rotation is an orthogonal matrix its inverse of transpose corresponds
            // to the matrix itself.
            Matrix3 modelInvT = Math::Algorithms::gaussJordanInverted(_model.rotationScaling()).transposed();
            return _rotation.rotationScaling() * modelInvT;
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

    const float plotX1 = 10.0f, plotX2 = 30.0f;
    const float plotY1 = -3.0f, plotY2 = 37.0f;

    // To obtain a cartesianGrid use CartesianGrid<NoTransform, 4>
    const CartesianGrid<NoTransform, 4> grid{Vector2{plotX1, plotY1}, Vector2{plotX2, plotY2}, 40, 40};

    const auto limits = mathFunctionFindExtrema(grid, f);
    fprintf(stderr, "extrema: %f,%f\n", limits.first, limits.second);

    Units xUnits{plotX1, plotX2}, yUnits{plotY1, plotY2};
    Units zUnits = Units{limits.first, limits.second};
    UnitsIterator uIt{zUnits};

    const char *uLabel;
    double uVal;
    while (uIt.next(uVal, uLabel)) {
        fprintf(stderr, "label: %s\n", uLabel);
    }

    float zMin = float(zUnits.mark_value(zUnits.begin()));
    float zMax = float(zUnits.mark_value(zUnits.end()));

    fprintf(stderr, "limits: %f,%f\n", zMin, zMax);

    const Trade::MeshData3D functionMeshData = mathFunctionMeshData(grid, f, df);
    const Trade::MeshData3D functionLines = mathFunctionLinesData(grid, f);
    const Trade::MeshData3D bottomGrid = unitsGridData(grid, xUnits, yUnits, zMin);

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

    _model = Matrix4::scaling({2.0f / grid.xSpan(), 2.0f / grid.ySpan(), 1.0f / (zMax - zMin)}) * Matrix4::translation({-grid.xCenter(), -grid.yCenter(), -zMin});

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
        .setNormalMatrix(normalsTransformation())
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

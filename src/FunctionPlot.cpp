#include <cmath>
#include <utility>
#include <Corrade/Utility/Format.h>
#include <Corrade/PluginManager/Manager.h>
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
#include <Magnum/Shaders/DistanceFieldVector.h>
#include <Magnum/Text/AbstractFont.h>
#include <Magnum/Text/DistanceFieldGlyphCache.h>
#include <Magnum/Text/Renderer.h>
#include "CartesianGrid.h"
#include "PlotConfig.h"
#include "Units.h"

using namespace Magnum;
using namespace Magnum::Math::Literals;

enum {
    AXIS_TICKS_NONE  = 0,
    AXIS_TICKS_X_INF = 1 << 0,
    AXIS_TICKS_X_SUP = 1 << 1,
    AXIS_TICKS_Y_INF = 1 << 2,
    AXIS_TICKS_Y_SUP = 1 << 3,
};

struct TextLabel {
    Containers::Pointer<Text::Renderer2D> textRenderer;
    Vector3 position;
};

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

static void addPlaneGridData(UnsignedInt& currentIndex, std::vector<UnsignedInt>& indices,
        std::vector<Vector3>& positions, const float x0, const float y0, const float x1, const float y1,
        const Units& unitsX, const Units& unitsY,
        const Vector3& xVector, const Vector3& yVector, const float z, UnsignedInt drawTicks, const float tickRelativeSize) {
    const Vector3 zVector = Math::cross(xVector, yVector);
    positions.push_back(x0 * xVector + y0 * yVector + z * zVector);
    positions.push_back(x0 * xVector + y1 * yVector + z * zVector);
    positions.push_back(x1 * xVector + y1 * yVector + z * zVector);
    positions.push_back(x1 * xVector + y0 * yVector + z * zVector);
    for (UnsignedInt i : {0, 1, 1, 2, 2, 3, 3, 0}) {
        indices.push_back(currentIndex + i);
    }
    currentIndex += 4;

    const float y0p = y0 - (drawTicks & AXIS_TICKS_X_INF ? tickRelativeSize * (y1 - y0) : 0.0f);
    const float y1p = y1 + (drawTicks & AXIS_TICKS_X_SUP ? tickRelativeSize * (y1 - y0) : 0.0f);

    UnitsIterator xIter{unitsX}, yIter{unitsY};
    double xval;
    const char *xlabel;
    while (xIter.next(xval, xlabel)) {
        const float xvalf = float(xval);
        if (xvalf < x0 || xvalf > x1) continue;
        positions.push_back(xvalf * xVector + y0p * yVector + z * zVector);
        positions.push_back(xvalf * xVector + y1p * yVector + z * zVector);
        indices.push_back(currentIndex++);
        indices.push_back(currentIndex++);
    }

    const float x0p = x0 - (drawTicks & AXIS_TICKS_Y_INF ? tickRelativeSize * (x1 - x0) : 0.0f);
    const float x1p = x1 + (drawTicks & AXIS_TICKS_Y_SUP ? tickRelativeSize * (x1 - x0) : 0.0f);

    double yval;
    const char *ylabel;
    while (yIter.next(yval, ylabel)) {
        const float yvalf = float(yval);
        if (yvalf < y0 || yvalf > y1) continue;
        positions.push_back(x0p * xVector + yvalf * yVector + z * zVector);
        positions.push_back(x1p * xVector + yvalf * yVector + z * zVector);
        indices.push_back(currentIndex++);
        indices.push_back(currentIndex++);
    }
}

template <typename F>
std::vector<TextLabel> axisLabels(const Units& units, const Vector3& axisVector, const Vector3& origin, const float x0, const float x1, F newTextRenderer, Text::Alignment alignement) {
    UnitsIterator iter{units};
    double xval;
    const char *xlabel;
    std::vector<TextLabel> labels;
    while (iter.next(xval, xlabel)) {
        const float xvalf = float(xval);
        if (xvalf < x0 || xvalf > x1) continue;
        TextLabel label{Containers::Pointer<Text::Renderer2D>{newTextRenderer(alignement)}, origin + xvalf * axisVector};
        label.textRenderer->reserve(12, GL::BufferUsage::DynamicDraw, GL::BufferUsage::StaticDraw);
        label.textRenderer->render(xlabel);
        labels.push_back(std::move(label));
    }
    return labels;
}

class MyApp: public Platform::Application {
    public:
        explicit MyApp(const Arguments& arguments);

    private:
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

        void renderAxisLabels(std::vector<TextLabel>& axisLabels);

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

        PluginManager::Manager<Text::AbstractFont> _manager;
        Containers::Pointer<Text::AbstractFont> _font;
        Text::DistanceFieldGlyphCache _cache;

        Shaders::DistanceFieldVector2D _textShader;
        Matrix3 _textViewportScaling;
        std::vector<TextLabel> _xAxisLabels, _yAxisLabels, _zAxisLabels;

        Matrix4 _rotation, _model, _projection;
        Vector2i _previousMousePosition;
        PlotConfig _plotConfig;
};

MyApp::MyApp(const Arguments& arguments):
    Platform::Application{arguments, Configuration{}.setTitle("Function Plot App"), GLConfiguration{}.setSampleCount(16)}, _cache(Vector2i(2048), Vector2i(512), 12)
{
    /* Load FreeTypeFont plugin */
    _font = _manager.loadAndInstantiate("FreeTypeFont");
    if(!_font) std::exit(1);

    /* Open the font and fill glyph cache */
    Utility::Resource rs("fonts");
    if(!_font->openSingleData(rs.getRaw("DejaVuSans.ttf"), 110.0f)) {
        Error() << "Cannot open font file";
        std::exit(1);
    }

    _font->fillGlyphCache(_cache, "xyzXYZ0123456789.eE+- ");

    GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
    GL::Renderer::enable(GL::Renderer::Feature::PolygonOffsetFill);
    GL::Renderer::setPolygonOffset(1.0f, 1.0f);
    // Multisampling is enabled by default.
    // GL::Renderer::enable(GL::Renderer::Feature::Multisampling);
    GL::Renderer::enable(GL::Renderer::Feature::Blending);
    GL::Renderer::setBlendFunction(GL::Renderer::BlendFunction::SourceAlpha, GL::Renderer::BlendFunction::OneMinusSourceAlpha);
    GL::Renderer::setBlendEquation(GL::Renderer::BlendEquation::Add, GL::Renderer::BlendEquation::Add);

    auto f = [](float x, float y) {
        return std::sin(x + y*y);
    };
    auto df = [](float x, float y) {
        const float c = std::cos(x + y*y);
        return Vector2{c, 2*y*c};
    };

    const float plotX1 = -3.0f, plotX2 = 3.0f;
    const float plotY1 = -2.0f, plotY2 = 2.0f;

    // To obtain a cartesianGrid use CartesianGrid<NoTransform, 4>
    const CartesianGrid<NoTransform, 8> grid{Vector2{plotX1, plotY1}, Vector2{plotX2, plotY2}, 80, 80};

    const auto limits = mathFunctionFindExtrema(grid, f);

    Units xUnits{plotX1, plotX2}, yUnits{plotY1, plotY2};
    Units zUnits = Units{limits.first, limits.second};

    const float zMin = float(zUnits.mark_value(zUnits.begin()));
    const float zMax = float(zUnits.mark_value(zUnits.end()));

    const float tickSize = _plotConfig.axisTickSize;
    auto newTextRenderer = [&](Text::Alignment alignment) { return new Text::Renderer2D(*_font, _cache, 0.02f, alignment); };
    _xAxisLabels = axisLabels(xUnits, Vector3::xAxis(), Vector3{0.0f, plotY1, zMin} - tickSize * (plotY2 - plotY1) * Vector3::yAxis(), plotX1, plotX2, newTextRenderer, Text::Alignment::TopCenter);
    _yAxisLabels = axisLabels(yUnits, Vector3::yAxis(), Vector3{plotX2, 0.0f, zMin} + tickSize * (plotX2 - plotX1) * Vector3::xAxis(), plotY1, plotY2, newTextRenderer, Text::Alignment::MiddleLeft);
    _zAxisLabels = axisLabels(zUnits, Vector3::zAxis(), Vector3{plotX1, plotY1, 0.0f} - tickSize * (plotY2 - plotY1) * Vector3::yAxis(), zMin, zMax, newTextRenderer, Text::Alignment::MiddleRight);

    const Trade::MeshData3D functionMeshData = mathFunctionMeshData(grid, f, df);
    const Trade::MeshData3D functionLines = mathFunctionLinesData(grid, f);

    std::vector<UnsignedInt> gridIndices;
    std::vector<Vector3> gridPositions;
    UnsignedInt gridPositionsIndex = 0;
    addPlaneGridData(gridPositionsIndex, gridIndices, gridPositions, plotX1, plotY1, plotX2, plotY2, xUnits, yUnits, Vector3::xAxis(), Vector3::yAxis(), zMin, AXIS_TICKS_X_INF|AXIS_TICKS_Y_SUP, tickSize);
    addPlaneGridData(gridPositionsIndex, gridIndices, gridPositions, zMin, plotX1, zMax, plotX2, zUnits, xUnits, Vector3::zAxis(), Vector3::xAxis(), plotY2, AXIS_TICKS_NONE, tickSize);
    addPlaneGridData(gridPositionsIndex, gridIndices, gridPositions, plotY1, zMin, plotY2, zMax, yUnits, zUnits, Vector3::yAxis(), Vector3::zAxis(), plotX1, AXIS_TICKS_Y_INF, tickSize);
    const auto gridData = Trade::MeshData3D{MeshPrimitive::Lines, gridIndices, {gridPositions}, {}, {}, {}, nullptr};

    _vertexBuffer.setData(MeshTools::interleave(functionMeshData.positions(0), functionMeshData.normals(0)));

    Containers::Array<char> indexData;
    MeshIndexType indexType;
    UnsignedInt indexStart, indexEnd;
    std::tie(indexData, indexType, indexStart, indexEnd) =
        MeshTools::compressIndices(functionMeshData.indices());
    _indexBuffer.setData(indexData);

    _vertexLinesBuffer.setData(functionLines.positions(0));
    _indexLinesBuffer.setData(functionLines.indices());

    _vertexGridBuffer.setData(gridData.positions(0));
    _indexGridBuffer.setData(gridData.indices());

    _mesh.setPrimitive(functionMeshData.primitive())
        .setCount(functionMeshData.indices().size())
        .addVertexBuffer(_vertexBuffer, 0, Shaders::Phong::Position{},
                                           Shaders::Phong::Normal{})
        .setIndexBuffer(_indexBuffer, 0, indexType, indexStart, indexEnd);

    _lines.setPrimitive(functionLines.primitive())
        .setCount(functionLines.indices().size())
        .addVertexBuffer(_vertexLinesBuffer, 0, Shaders::Flat3D::Position{})
        .setIndexBuffer(_indexLinesBuffer, 0, MeshIndexType::UnsignedInt, 0, functionLines.indices().size());

    _baseGrid.setPrimitive(gridData.primitive())
        .setCount(gridData.indices().size())
        .addVertexBuffer(_vertexGridBuffer, 0, Shaders::Flat3D::Position{})
        .setIndexBuffer(_indexGridBuffer, 0, MeshIndexType::UnsignedInt, 0, gridData.indices().size());

    _model = Matrix4::scaling({2.0f / grid.xSpan(), 2.0f / grid.ySpan(), 0.5f / (zMax - zMin)}) * Matrix4::translation({-grid.xCenter(), -grid.yCenter(), -zMin});

    _rotation = Matrix4::rotationX(-60.0_degf);
    _projection =
        Matrix4::perspectiveProjection(
            35.0_degf, Vector2{windowSize()}.aspectRatio(), 0.01f, 100.0f)*
        Matrix4::translation(Vector3::zAxis(- _plotConfig.cameraDistance));

    _textViewportScaling = Matrix3::scaling(Vector2::yScale(Vector2(GL::defaultFramebuffer.viewport().size()).aspectRatio()));

    _plotConfig.overwriteGrids = true;
}

void MyApp::renderAxisLabels(std::vector<TextLabel>& axisLabels) {
    for (auto& label : axisLabels) {
        Vector3 posNorm = _model.transformPoint(label.position);
        Vector3 textCoord = (_projection * _rotation).transformPoint(posNorm);
        Vector2 projTextCoord{textCoord.x(), textCoord.y()};
        _textShader.setTransformationProjectionMatrix(Matrix3::translation(projTextCoord) * _textViewportScaling);
        label.textRenderer->mesh().draw(_textShader);
    }
}

void MyApp::drawEvent() {
    GL::defaultFramebuffer.clear(GL::FramebufferClear::Depth).clearColor(Color4{1.0f});

    _textShader.bindVectorTexture(_cache.texture());
    _textShader.setColor(Color3{0.3f}).setSmoothness(0.2f);

    renderAxisLabels(_xAxisLabels);
    renderAxisLabels(_yAxisLabels);
    renderAxisLabels(_zAxisLabels);

    _flatShader.setTransformationProjectionMatrix(_projection*transformation());

    if (_plotConfig.drawGrids) {
        _flatShader.setColor(Color4{0.7f});
        _baseGrid.draw(_flatShader);
    }

    if (_plotConfig.overwriteGrids) {
        GL::defaultFramebuffer.clear(GL::FramebufferClear::Depth);
    }

    if (_plotConfig.drawSurface) {
        Color3 color = _plotConfig.surfaceColor;
        _shader.setLightPosition({7.0f, 5.0f, 2.5f})
            .setLightColor(Color3{1.0f})
            .setDiffuseColor(color)
            .setAmbientColor(Color3::fromHsv(color.hue(), 1.0f, 0.4f))
            .setTransformationMatrix(transformation())
            .setNormalMatrix(normalsTransformation())
            .setProjectionMatrix(_projection);
        _mesh.draw(_shader);
    }

    if (_plotConfig.drawSurfaceLines) {
        _flatShader.setColor(0x555555_rgbf);
        _lines.draw(_flatShader);
    }

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

#pragma once

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

struct TextLabel {
    Corrade::Containers::Pointer<Text::Renderer2D> textRenderer;
    Magnum::Vector3 position;
};

class FunctionPlot {
    public:
        explicit FunctionPlot(Magnum::Vector2i windowSize);

        void setPosition(const Magnum::Vector2i& p);
        void setMovePosition(const Magnum::Vector2i& p);
        void draw();

    private:
        Magnum::Matrix4 transformation() const {
            return _rotation * _model;
        }

        Magnum::Matrix3 normalsTransformation() const {
            // For the normals we need the inverse of the transpose of the matrix, excluding the
            // translation term.
            // As the rotation is an orthogonal matrix its inverse of transpose corresponds
            // to the matrix itself.
            Magnum::Matrix3 modelInvT = Magnum::Math::Algorithms::gaussJordanInverted(_model.rotationScaling()).transposed();
            return _rotation.rotationScaling() * modelInvT;
        }

        void renderAxisLabels(std::vector<TextLabel>& axisLabels);

#if 0
        void drawEvent() override;
        void mousePressEvent(MouseEvent& event) override;
        void mouseReleaseEvent(MouseEvent& event) override;
        void mouseMoveEvent(MouseMoveEvent& event) override;
#endif

        Magnum::GL::Buffer _indexBuffer, _indexLinesBuffer, _vertexBuffer, _vertexLinesBuffer;
        Magnum::GL::Buffer _indexGridBuffer, _vertexGridBuffer;
        Magnum::GL::Mesh _mesh, _lines;
        Magnum::GL::Mesh _baseGrid;
        Magnum::Shaders::Phong _shader;
        Magnum::Shaders::Flat3D _flatShader;

        Magnum::PluginManager::Manager<Text::AbstractFont> _manager;
        Magnum::Containers::Pointer<Text::AbstractFont> _font;
        Magnum::Text::DistanceFieldGlyphCache _cache;

        Magnum::Shaders::DistanceFieldVector2D _textShader;
        Magnum::Matrix3 _textViewportScaling;
        std::vector<TextLabel> _xAxisLabels, _yAxisLabels, _zAxisLabels;

        Magnum::Matrix4 _rotation, _model, _projection;
        Magnum::Vector2i _previousMousePosition;
        PlotConfig _plotConfig;
};

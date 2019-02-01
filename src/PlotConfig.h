#pragma once

#include <Magnum/Math/Color.h>

using namespace Magnum;
using namespace Magnum::Math::Literals;

struct PlotConfig {
    PlotConfig(): overwriteGrids(false), surfaceColor(0xff9900_rgbf), drawGrids(true), drawSurface(true), drawSurfaceLines(true) { }

    bool overwriteGrids;
    Color3 surfaceColor;
    bool drawGrids;
    bool drawSurface;
    bool drawSurfaceLines;
};

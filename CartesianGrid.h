#pragma once
#include <cmath>
#include <limits>
#include <Magnum/Math/Vector2.h>

struct NoTransform {
    static inline void transform(float& x, float& y) { }
};

struct CircleTransform {
    static inline void transform(float& x, float& y) {
        const float absX = fabs(x), absY = fabs(y);
        const float norm = sqrt(absX*absX + absY*absY);
        if (norm > std::numeric_limits<float>::epsilon()) {
            const float ratio = (absX > absY ? absX / norm : absY / norm);
            x *= ratio;
            y *= ratio;
        }
    }
};

template <typename Transformer = NoTransform, unsigned int LinesRatio = 1>
class CartesianGrid {
public:
    CartesianGrid(Magnum::Vector2 llc, Magnum::Vector2 urc, int nx, int ny):
        _llc(llc), _urc(urc), _nx(nx), _ny(ny) {
    }

    int pointsCount() const {
        return (_nx + 1) * (_ny + 1);
    }

    template <typename F>
    void spanGridPoints(F eval) const {
        for (int j = 0; j <= _ny; j++) {
            for (int i = 0; i <= _nx; i++) {
                float x = _llc.x() + i * (_urc.x() - _llc.x()) / _nx;
                float y = _llc.y() + j * (_urc.y() - _llc.y()) / _ny;;
                Transformer::transform(x, y);
                eval(x, y);
            }
        }
    }

    int trianglesCount() const {
        return 2 * _nx * _ny;
    }

    template <typename F>
    void spanTrianglesIndices(F eval) const {
        for (int j = 0; j <= _ny - 1; j++) {
            for (int i = 1; i <= _nx; i++) {
                eval(this->pointIndex(i - 1, j), this->pointIndex(i - 1, j + 1), this->pointIndex(i, j));
                eval(this->pointIndex(i - 1, j + 1), this->pointIndex(i, j + 1), this->pointIndex(i, j));
            }
        }
    }

    int linesCount() const {
        return (2 * _nx * (_ny + 1) + 2 * (_nx + 1) * _ny) / LinesRatio;
    }

    template <typename F>
    void spanLinesIndices(F eval) const {
        for (int j = 0; j <= _ny; j += LinesRatio) {
            for (int i = 1; i <= _nx; i++) {
                eval(this->pointIndex(i - 1, j), this->pointIndex(i, j));
            }
        }
        for (int j = 1; j <= _ny; j++) {
            for (int i = 0; i <= _nx; i += LinesRatio) {
                eval(this->pointIndex(i, j - 1), this->pointIndex(i, j));
            }
        }
    }

private:
    int pointIndex(int i, int j) const {
        return j * (_nx + 1) + i;
    }

    const Magnum::Vector2 _llc, _urc;
    const int _nx, _ny;
};

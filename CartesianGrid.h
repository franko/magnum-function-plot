#pragma once

#include <Magnum/Math/Vector2.h>

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
                const float x = _llc.x() + i * (_urc.x() - _llc.x()) / _nx;
                const float y = _llc.y() + j * (_urc.y() - _llc.y()) / _ny;;
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
        return 2 * _nx * (_ny + 1) + 2 * (_nx + 1) * _ny;
    }

    template <typename F>
    void spanLinesIndices(F eval) const {
        for (int j = 0; j <= _ny; j++) {
            for (int i = 1; i <= _nx; i++) {
                eval(this->pointIndex(i - 1, j), this->pointIndex(i, j));
            }
        }
        for (int j = 1; j <= _ny; j++) {
            for (int i = 0; i <= _nx; i++) {
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

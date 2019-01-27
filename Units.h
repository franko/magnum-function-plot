
/* Units.h
 *
 * Copyright (C) 2009, 2010 Francesco Abbate
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#pragma once

struct LabelIterator {
    virtual bool next(double& val, const char*& text) = 0;
    virtual ~LabelIterator() {}
};

class Units {
public:
    enum { label_format_max_size = 16 };
    enum format_e { format_int, format_float, format_invalid };

    Units() { clear(); }
    Units(double min, double max, double spacefact = 4.0) {
        set(min, max, spacefact);
    }

    void clear() {
        _major = 1;
        _order = 0;
        _dmajor = 1;
        _inf = 0;
        _sup = 1;
        _decimals = 0;
    }

    void set(double min, double max, double spacefact = 4.0);

    int begin() const {
        return _inf;
    };
    int end()   const {
        return _sup;
    };

    void limits(int &start, int &fin, double &step) const
    {
        start = _inf;
        fin = _sup;
        step = _dmajor;
    };

    void mark_label (char *label, unsigned size, int mark) const;

    double mark_value (int mark) const {
        return _dmajor * mark;
    };

    double mark_scale(double x);

private:
    int _major;
    int _order;
    double _dmajor; // equal to (_major * 10^order)
    int _inf, _sup; // expressed in the base of (_major * 10^order)
    int _decimals; // Number of decimal units.
};

class UnitsIterator : public LabelIterator {
public:
    UnitsIterator(const Units& u): _units(u) {
        _index = u.begin();
    }

    virtual bool next(double& val, const char*& text)
    {
        if (_index > _units.end())
            return false;

        _units.mark_label(_buffer, 32, _index);

        val = _units.mark_value(_index);
        text = _buffer;
        _index ++;
        return true;
    }

private:
    char _buffer[32];
    int _index;
    const Units& _units;
};

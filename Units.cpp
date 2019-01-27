
/* Units.cpp
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

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Units.h"
#include "printf_check.h"

void Units::set(double yinf, double ysup, double spacefact)
{
    double del;

    if (ysup == yinf)
        ysup = yinf + 1.0;

    del = (ysup - yinf) / spacefact;

    _order = int(floor(log10(del)));

    double expf = pow(10, _order);
    double delr = del / expf;

    if (5 <= delr)
        _major = 5;
    else if (2 <= delr)
        _major = 2;
    else
        _major = 1;

    _inf = int(floor(yinf / (_major * expf) + 1e-5));
    _sup = int(ceil (ysup / (_major * expf) - 1e-5));

    _decimals = (_order < 0 ? -_order : 0);

    _dmajor = _major * expf;
}

void Units::mark_label (char *lab, unsigned size, int mark) const
{
    bool minus = (_inf < 0);
    int asup = (minus ? -_inf : _sup);
    char fmt[16];

    if (size < 16)
        return;

    if (_decimals == 0)
    {
        snprintf (lab, size, "%.0f", mark * _dmajor);
        lab[size-1] = '\0';
    }
    else
    {
        int dec = (_decimals < 10 ? _decimals : 9);
        int base = int(floor(asup * _dmajor));
        int space = dec + (base > 0 ? int(log10(base)) : 0) + 1 + (minus ? 1 : 0) + 1;
        snprintf (fmt, 16, "%%%i.%if", space, dec);
        fmt[15] = '\0';
        snprintf (lab, size, fmt, mark * _dmajor);
        lab[size-1] = '\0';
    }
}

double Units::mark_scale (double x)
{
    double xinf = _inf * _dmajor, xsup = _sup * _dmajor;
    return (x - xinf) / (xsup - xinf);
}

void Units::fmt_label(char* label, unsigned size, format_e tag, const char* fmt, int mark) const
{
    double val = mark_value(mark);
    switch (tag)
    {
    case format_int:
    {
        unsigned nchars = snprintf(label, size, fmt, int(val));
        if (nchars >= size)
            label[size-1] = 0;
        break;
    }
    case format_float:
    {
        unsigned nchars = snprintf(label, size, fmt, val);
        if (nchars >= size)
            label[size-1] = 0;
        break;
    }
    default:
        memcpy(label, "*", 2);
    }
}

Units::format_e Units::parse_label_format(const char* fmt)
{
    if (strlen(fmt) >= label_format_max_size)
        return format_invalid;

    arg_type_e arg_type;
    const char* tail;
    int n = check_printf_argument(fmt, tail, arg_type);
    if (n != 1)
        return format_invalid;
    else
    {
        const char* tt;
        arg_type_e aa;
        if (check_printf_argument(tail, tt, aa) != 0)
            return format_invalid;
    }

    return (arg_type == argument_int ? format_int : format_float);
}

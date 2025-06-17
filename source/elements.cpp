/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * elements.cpp
 * Copyright (C) 2019 Kazimierz Skrobas <kskrobas@unipress.waw.pl>
 *
 * npcl is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * npcl is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "elements.h"


namespace  Elements {

    const std::map<std::string, std::string> mass{
        {"Al", "26.982"},
        {"H",  "1.00794"},
        {"C",  "12.0107"},
        {"Ca", "40.078"},
        {"Cd", "112.411"},
        {"Ce", "140.116"},
        {"Co", "58.933195"},
        {"Cr", "51.996"},
        {"Cu", "63.546"},
        {"Fe", "55.845"},
        {"F",  "18.9984"},
        {"Ga", "69.723"},
        {"N",  "14.007"},
        {"Na", "22.9897"},
        {"Ni", "58.6934"},
        {"Mg", "24.305"},
        {"Mn", "54.94"},
        {"Mo", "95.94"},
        {"O",  "15.9994"},
        {"P",  "30.973762"},
        {"S",  "32.065"},
        {"Se", "78.96"},
        {"Si", "28.0855"},
        {"U",  "238.028913"},
        {"Te", "127.6"},
        {"Ti", "47.867"},
        {"W",  "183.84"},
        {"Xe", "131.29"},
        {"Zn", "65.32"},
        {"Zr", "91.22"},
        {"H2O","18.015"}

    };

/////https://slideplayer.com/slide/6143092/
/// ///https://www.chemie-biologie.uni-siegen.de/ac/be/lehre/ws1213/x-ray_powder_diffraction_(xrpd).pdf
    const std::map<std::string, std::string> xrad{

        {"Ag", "0.561"},
        {"Co", "1.7929"},
        {"Cr", "2.293663"},
        {"Cu", "1.54059"},
        {"Fe", "1.936"},
        {"Mo", "0.7136"},
        {"W",  "0.2106"},

    };
}

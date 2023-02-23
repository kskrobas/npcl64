#include "colormsg.h"

/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * colormsg.cpp
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

void errMsg(const string msg)
{
    cerr<<cred<<"ERROR: "<<cdef<<msg<<endl;
}


void warnMsg(const std::string msg)
{
    cerr<<cblue<<"WARNING: "<<cdef<<msg<<endl;
}


void infoMsg(const string msg)
{
    cerr<<cyellow<<"INFO: "<<cdef<<msg<<endl;
}

void logMsg(const string &msg)
{
    cerr<<cgreen<<"     "<<msg<<cdef<<endl;
}


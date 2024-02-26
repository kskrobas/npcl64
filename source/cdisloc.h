/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cdisloc.h
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
#ifndef CDISLOC_H
#define CDISLOC_H

#include <vector>
#include "nanograin.h"


using namespace std;

typedef const string cstring;

class Cdisloc{
private:


    void clearData();

    void insertLoop();
    void rotateLoop();
    void rpyLoop();
    void validateDisloc(NanoGrain::vatoms &dislocAtoms);

public:
    enum Edisstatus{OK, ERR_NOFPARAMS};

vector<NanoGrain::StAtomType> atomTypes;
string axis,axispos,rangeR,rangeA;
string angle,projh;
string mindist,scatter,mode;
string saveopt;
vector<string> vrpy;

string fileNameIn,fileNameHeader;
vector<string> fileNameOut;

NanoGrain::StNanoGrain *grain;
stdumap *ptr_uvar;




    Cdisloc();

    int findAtomName(cstring &aname__);
    void addAtomName(cstring &aname__);

    bool parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__);
    void calc();






};




#endif // CDISLOC_H

/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cavepdh.h
 * Copyright (C) 2019 Kazimierz Skrobas <Kazimierz.Skrobas@unipress.waw.pl
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

#ifndef CAVEPDH_H
#define CAVEPDH_H



#include <iostream>
#include <string>


#include "scriptanalyser.h"
#include "cpdh.h"



class Cavepdh
{
private:
    stdumap *ptr_uvar;
    std::string openFiles,fileNameOut;
    string weight,numberOfFiles2Ave;
    size_t numberOfAtomTypes;
    vector<string> atomTypes;
    double binWidth,minBin,maxBin;


    size_t numOfAtoms;
    vector<dataPdh> dataYnn;
    //vector<string> randFileNames;

    bool printprm;

    void saveFilePdhs();
    void saveFileDat();
    void saveFileLhs();

    void printPrm();

    void exportToCpdh();
public:
    enum EAvePdhStatus{OK,RUN,ERR_OPEN,ERR_LIST} status;

    Cavepdh();

    Cpdh * pdh;

     bool parseCommands(vcmdlist &cmd,size_t &index, stdumap *uvar__);

     void calc();
     //void calc0();

     void clearData();

};

#endif // CAVEPDH_H

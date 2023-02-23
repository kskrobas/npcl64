/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cdiff.h
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

#ifndef CGR_H
#define CGR_H

#include "cdiff.h"


#ifdef __linux__

//#elif _WIN32 or _WIN64
#else
    # define M_PI		3.14159265358979323846	/* pi */
    # define M_1_PI	0.318309886183790671537767526745028724L /* 1/pi */
#endif



typedef StDataGeneric<position> dataGr;
//=============================================================================
class CWindowFunction;
//=============================================================================

class Cgr
{
private:

     stdumap *ptr_uvar;
     std::time_t startTime,stopTime;


    void saveResults();
    void saveDatFile(const std::string &fileName);
    void saveGrFile (const std::string &fileName);
    void saveBinFile(const std::string &fileName);

    enum EGr{OK,ERR_FILEOPEN,ERR_EMPTY_DIFF};

    void dispPrm();
    void normGr();

    friend ostream & operator<< (ostream &, Cgr &cgr) ;

    void runCalc();
    void allocMem(const size_t &size);
    void clearMem();

#ifdef __linux__
    void plotFigure();
    void saveFigure(const std::string &fileName);

    void plotData(FILE *plotHandle);
#endif

public:
    Cdiff *diff;
    Cgr();

    std::string range,wf,comment;
    std::string threads;
    std::string correctionFactor;    
    vector<std::string> fileNameOut;
    bool diffTime,printPrm,norm;

    dataGr dataX,dataY;


#ifdef __linux__
    std::string plotPrm;
    std::string figFileName;
    std::string figWidth,figHeight;
#endif


    bool parseCommands(vcmdlist &cmd,size_t &index, stdumap *uvar__);
    void calc();
    void presetData()
    {
        range="0 0.125";
        wf="box";
        comment.clear();
        threads="1";
        correctionFactor="1";


        fileNameOut.clear();
        diffTime=false;
        printPrm=false;
        norm=false;

#ifdef __linux__
        plotPrm.clear();
        figFileName.clear();
#endif
    }

};

//=============================================================================
class CWindowFunction{
private:

public:
        virtual double operator()(const double &x)=0;
        virtual void setPrm(const double &p)=0;

        virtual ~CWindowFunction(){ }
};
//'''''''''''''''''''''''''''''''''''''''''''''''''''''''
class CWindowLorch: public CWindowFunction{

private:
        double k0;
public:
        CWindowLorch(){ }
        void setPrm(const double &p){k0=M_PI/p;}
        double operator()(const double &x){ const double arg= k0*x; return sin(arg)/arg;}
};
//'''''''''''''''''''''''''''''''''''''''''''''''''''''''
class CWindowBoxCar: public CWindowFunction{
public:
        CWindowBoxCar(){ }
        void setPrm(const double &){ }
        double operator()(const double &x){ (void) x;  return 1.0;}

};

#endif // CGR_H

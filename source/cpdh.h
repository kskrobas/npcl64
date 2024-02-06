/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cpdh.h
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


#ifndef CPDH_H
#define CPDH_H


#include "nanograin.h"
#include "stdatagen.h"

#include <ctime>


typedef StDataGeneric<position> dataPdh;
typedef const size_t csize;

class Cpdh {

private:
    bool printPrm,diffTime;
    position maxDist;
    stdumap *ptr_uvar;


    std::string fileName,fileNameIn;
    std::string comment;
    std::string title;
    std::string range;


    position getBinWidth();


    void (Cpdh::*ptrCalc)();
    struct StMultiPdhPrm{
        size_t ifrom,ito,jfrom,jto,size;
        position wbin,ibin;
    };

    void monoLatticePdh();    
    //void biLatticePdh();
    void multiLatticePdh();
	
	void pdhByCppThreads();

    bool ijLatticePdh(StMultiPdhPrm *prm,vector<dataPdh>::iterator & buffPdh,const std::string &);


    void saveResults();
    void openFile();


    void saveDatFile();
    void savePdhlFile();
    void savePdhsFile();
    void saveLhsFile();
    void saveExtFile();
    void saveBinFile();
    void saveBinFilePLS();


    void openPdhDat();
    void openPdhls();
    //void openPdhs();
    void openPdhLhs();
    void openBinFilePLS();



    std::time_t startTime,stopTime;

    void dispParams();

    friend ostream & operator<< (ostream &, Cpdh &cpdh) ;

#ifdef __linux__
    void plotFigure();
    void saveFigure(const std::string &fileName);

    void plotData(FILE *plotHandle);
#endif


public:
    Cpdh();

    enum Epdhstatus{OK,RUN,ERR_EMPTYGRAIN,ERR_FILEOPEN} status;

    NanoGrain::StNanoGrain *grain;
    dataPdh dataX;
    dataPdh dataYii;
    vector<dataPdh> dataYnn;

#ifdef __linux__
    std::string plotPrm;
    std::string figFileName;
    std::string figWidth,figHeight;
#endif

    std::string threads;
    std::string bin;

    enum EBinMode{SINGLE,DOUBLE};
    EBinMode binMode=SINGLE;
	
	enum EMthMode{OPENMP,STDTHREAD};
	EMthMode mthmode=OPENMP;


     bool parseCommands(vcmdlist &cmd,size_t &index, stdumap *uvar__);
     void calc();
     void import(dataPdh &dataX,dataPdh &dataY);


     size_t dataSize(){return dataX.size; }


     void clearData(){
         bin="0.03125";
         binMode=SINGLE;
         threads="1";
         fileName="";
         fileNameIn.clear();
         range.clear();
         comment="";
         diffTime=false;
         printPrm=false;

         dataX.freeMem();
         dataYii.freeMem();

         startTime=0;
         stopTime=0;                  

         title="results of PDH calculations";

#ifdef __linux__
        plotPrm.clear();
        figFileName.clear();
#endif

     }


};


//struct Cpdh::StMultiPdhPrm{
//size_t ifrom,ito,jfrom,jto,size;
//position wbin,ibin;
//};


#endif // CPDH_H

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



#ifndef CDIFF_H
#define CDIFF_H


#include "cpdh.h"
#include <bitset>
#include <complex>




typedef StDataGeneric<position> dataDiff;
typedef std::complex<position> cxposition ;
typedef std::vector< std::vector< std::vector <cxposition> >> cx3dspace;

class Cdiff
{
private:

    bool scfact;
    position maxDist;
    stdumap *ptr_uvar;

    enum Ediffstatus{OK,RUN,ERR_EMPTYPDH,ERR_SCATCOEFF,
                    ERR_RANGE,ERR_FILEOPEN,ERR_ATYPE_NUM,ERR_FFORMAT,ERR_FOPEN} diffstatus;

    std::time_t startTime,stopTime;

    // debyea diffraction (powder)
    void debDiff_mono();
    void debDiff_multi();

    // laue diffraction
    void laueDiff();

    //vector< vector<position> > pdhNZBdist;

    struct StNonZeroBins{
    vector<position> dataX;
    vector<position> dataY;

            void reserveMemory(const size_t &size)
            {
                    dataX.reserve(size);
                    dataY.reserve(size);
            }


            void push_values(const position *x__,const position *y__)
            {
                    dataX.push_back(*x__);
                    dataY.push_back(*y__);
            }

            void shrink()
            {
                    dataX.shrink_to_fit();
                    dataY.shrink_to_fit();
            }

    };

    vector<StNonZeroBins> pdhNZB;
    void buildNonZeroBinsPdh();

    struct StSaveOpt{
    enum EAXIS{X,Y,Z} axis;
    enum EFORMAT{SHORT,LONG} format;
    enum EBITSAVE{BT,BQ,BI,BS};
    std::bitset<4> dataToSave;
    int layer;
    void set(){dataToSave.set(); axis=Z; format=SHORT; layer=-1;}

    } saveopt;

    void saveResults();
    void saveDatFile(const string &fileName);
    void saveDiffFile(const string &fileName);
    void saveBinFile(const string &fileName);
    void saveLaueFile(const string &fileName);
    void saveLaueFileBin(const string &fileName);

    void openFile();
    void openDiffFile();

    void allocMem(const size_t &size);
    void clearMem();

    bool loadScattFactors(const string &aname,const int index);
    union Uscattfactors{
        struct{float a1,b1,a2,b2,a3,b3,a4,b4,c,ncs,d1,e1,d2,e2,d3,e3,d4,e4,d5,e5;};
        float sf_array[20];
    } *scfactors;



    double (Cdiff::*ptr_scattfactor)(const double &q, const size_t index);
    double noScatFact(const double &q, const size_t index){ (void) q; (void) index; return 1.0;}
    double xrayScatEqu(const double &q, const size_t index);
    double electronScatEqu(const double &q, const size_t index);
    double neutronCrossSec(const double &q, const size_t index);



    void polarizationI();
    void debyeaWaller();
    void noiseSQ();
    void normIS();


    void dispParams();
    void plot();

    char radiationToChar();

    friend ostream & operator<< (ostream &, Cdiff &cdiff) ;

public:
    Cdiff();
    ~Cdiff(){ if(scfactors) delete [] scfactors;}

    Cpdh *pdh;

    bool printPrm,diffTime,extrapolate,theta;
    bool fastsinc,norm,polarization,u2dendl;
    std::string range,lambda,dw,ki;

    std::string threads;
    std::string noise;
    std::string fileNameIn,fileSft;
    vector<std::string> fileNameOut;
    std::string comment;
    std::string plotPrm;
    enum EMODE {debyea,laue} mode;
    enum ERADIATION {OFF,XRAY,NEUTR,ELECT} radiation;

    dataDiff t,q;
    dataDiff I,S;

    vector<position> kpoints;
    cx3dspace laueDiffData;


    bool parseCommands(vcmdlist &cmd,size_t &index, stdumap *uvar__);
    void calc();
    void presetData()
    {

        threads="1";

        if(scfactors)
        delete [] scfactors;

        extrapolate=true;
        printPrm=false;
        diffTime=false;
        theta=true;
        fastsinc=false;
        norm=false;
        polarization=false;
        u2dendl=true;


        fileNameIn.clear();
        fileNameOut.clear();
        comment.clear();
        fileSft.clear();
        pdhNZB.clear();
        noise.clear();
        dw.clear();
        ki=std::string("0 0 1");
        plotPrm.clear();


        startTime=0;
        stopTime=0;

        saveopt.set(); //save all: 2θ,Q,I,S

        radiation=XRAY;
        mode=debyea;

    }


};

#endif // CDIFF_H

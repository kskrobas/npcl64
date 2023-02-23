/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * Cmergepdh.h
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


#ifndef CMERGEPDH_H
#define CMERGEPDH_H

#include "cpdh.h"

//-----------------------------------------------------------------------------
struct StNameRadii{
    string name;double radii;
    static double maxRadii;
    StNameRadii(const string &__name):name(__name) { }

    void num2radii(const size_t N)
    {
        radii=std::pow(N,1.0/3.0)/0.901182;///wartosc pobrana z fitu N(r)
        if (radii>maxRadii)maxRadii=radii;
    }

};

//-----------------------------------------------------------------------------

class Cmergepdh
{
private:
    stdumap *ptr_uvar;
    std::string openFiles,fileNameOut;
    std::string weight,flimit;
    std::string from,step,to;
    dataPdh sumBinij;
    size_t totalNumberOfAtoms;
    size_t numberOfAtomTypes;
    double totalMinBin,totalMaxBin,binWidth;
    vector<StNameRadii> pdhfiles;


    bool printprm;

    bool getPdhsStatistic(vector < vector<StNameRadii> > &vcells, vector<StNameRadii> &pdhfiles);
    bool getPlsStatistic(vector < vector<StNameRadii> > &vcells, vector<StNameRadii> &pdhfiles);
    void mergePdhsBlockPdhs(const vector<StNameRadii> &pdhsBlock,const vector<int> & id );
    void mergePdhsBlockPls(const vector<StNameRadii> &pdhsBlock,const vector<int> & id );

    void saveFilePdhs();
    void exportToCpdh();
public:
    enum EmergePdhStatus{OK,RUN,ERR_OPEN,ERR_WEIGHT,ERR_RANGE} status;

    Cmergepdh();
    ~Cmergepdh();

    Cpdh * pdh;

     bool parseCommands(vcmdlist &cmd,size_t &index, stdumap *uvar__);

     void calc();
     void clearData();
};

//-----------------------------------------------------------------------------


#endif // CMERGEPDH_H

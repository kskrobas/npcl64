/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cavepdh.cpp
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

#include "cavepdh.h"
#include "cprogress.h"
#include "crandom.h"
#include "createdir.h"

#include <iomanip>
#include <ctime>
#include <stdexcept>
#include <chrono>
#include <thread>
#include <future>

using namespace std;
typedef  const double cdouble ;


#ifdef DEBUG
#define DB true
#else
#define DB false
#endif

#ifdef __linux__

#else
#define M_SQRT2  1.414213562373095	
#endif

//===================================================================================

Cavepdh::Cavepdh()
{


}

//===================================================================================

void Cavepdh::saveFilePdhs()
{
fstream fout(fileNameOut,ios::out);

                if(!fout){
                    cerr<<"ERROR: couldn't save data to file: "<<fileNameOut<<endl;
                return;
                }

const size_t dataSize=dataYnn.front().size;
csize numOfAtomTypes= numberOfAtomTypes;

csize mpdhSize=(numOfAtomTypes+1)*numOfAtomTypes/2;
std::time_t datetime = std::time(nullptr);
size_t i,j;
csize nprec=11;
csize colwh=12;

const int numOfDigits=static_cast<int> (std::log10(dataSize)+2);
std::streampos strpos;
size_t numOfNonZeroBins=0;
const double iBinWidth=1.0/binWidth;

                //cout<<"numOfDigits: "<<numOfDigits<<endl;

                fout<<"#ver: 1"<<endl;
                fout<<"#title: avePdh"<<endl;
                fout<<"#sizeRC: "; strpos=fout.tellp();fout<<setw(numOfDigits)<<" "<<(1+mpdhSize)<<endl;
                fout<<"#date: "<<std::asctime(std::localtime(&datetime));
                fout<<"#atomTypes: "<<numOfAtomTypes; for(auto &aname : atomTypes) fout<<" "<<aname;

                fout<<endl;

                fout<<"#atomsNumber: "<<numOfAtoms<<endl;
                fout<<"#atomsXYZ: ignore"<<endl;
                fout<<"#binWidth: "<<binWidth<<endl;
                fout<<"#comment: WARNING: number of atoms is approximated; (minBin,maxBin)=( ";
                    fout<<setprecision(nprec)<<setw(colwh)<<minBin<<", ";
                    fout<<setprecision(nprec)<<setw(colwh)<<maxBin<<")"<<endl;
                fout<<"#";


                /// wypisuje nazwy column
                fout<<setw(colwh)<<'X';

                for (i=0;i<numOfAtomTypes;i++){
                    for(j=i;j<numOfAtomTypes;j++)
                            fout<<" "<<setw(colwh)<<std::string(atomTypes[i]+"-"+atomTypes[j]);

                }

                /// koniec naglowka zaznaczony #
                fout<<endl;
                fout.fill('#');
                fout<<setw((colwh+1)*(mpdhSize+1))<<'#'<<endl;///13=set(12)+space
                fout.fill(' ');


                if(numberOfAtomTypes==1){
                double x;                
                auto  &sumBinij=dataYnn[0];

                    for(i=0;i<dataSize;i++){                        
                        if(sumBinij[i]>0){
                            x=minBin+i*iBinWidth;
                            fout<<" "<<setprecision(nprec)<<setw(colwh)<<x<<" "<<setprecision(nprec)<<setw(colwh)<<sumBinij[i]<<endl;
                            numOfNonZeroBins++;
                        }
                    }
                }
                else{
                double x;
                const size_t numOfBinTypes=dataYnn.size();
                bool nonZeroBins=false;

                        for(i=0;i<dataSize;i++){

                            for(j=0;j<numOfBinTypes;j++)
                                if(dataYnn[j][i]>0){
                                    nonZeroBins=true;
                                    break;
                                }


                            if(nonZeroBins){
                                nonZeroBins=false;
                                 x=minBin+i*iBinWidth;
                                 fout<<" "<<setprecision(nprec)<<setw(colwh)<<x;

                                 for(j=0;j<numOfBinTypes;j++)
                                     fout<<" "<<setprecision(nprec)<<setw(colwh)<<dataYnn[j][i];

                                 fout<<endl;
                                 numOfNonZeroBins++;
                            }
                        }

                }//else's end

                fout.seekp(strpos);
                fout<<numOfNonZeroBins<<"\t";

                fout.close();
}

//===================================================================================
void Cavepdh::saveFileDat()
{
fstream fout(fileNameOut,ios::out);

                if(!fout){
                    cerr<<"ERROR: couldn't save data to file: "<<fileNameOut<<endl;
                return;
                }


double x;
size_t i,j;
const size_t dataSize=dataYnn.front().size;
csize numOfAtomTypes= numberOfAtomTypes;
const double iBinWidth=1.0/binWidth;


                for(i=0;i<dataSize;i++){
                    x=minBin+i*iBinWidth;
                    fout<<x<<"   ";

                    for(j=0;j<numOfAtomTypes;j++)
                        fout<<dataYnn[j][i]<<endl;
                }


                fout.close();
return;
}

//===================================================================================
void Cavepdh::saveFileLhs()
{
fstream fout(fileNameOut,ios::out);

                if(!fout){
                    cerr<<"ERROR: couldn't save data to file: "<<fileNameOut<<endl;
                return;
                }

const size_t dataSize=dataYnn.front().size;
csize numOfAtomTypes= numberOfAtomTypes;

csize mpdhSize=(numOfAtomTypes+1)*numOfAtomTypes/2;
std::time_t datetime = std::time(nullptr);
size_t i,j;
csize nprec=11;
csize colwh=12;

const int numOfDigits=static_cast<int> (std::log10(dataSize)+2);
std::streampos strpos;
size_t numOfNonZeroBins=0;
const double iBinWidth=1.0/binWidth;

                    fout<<"#date: "<<std::asctime(std::localtime(&datetime));
                    fout<<"#sizeRC: "; strpos=fout.tellp(); fout<<setw(numOfDigits)<<" "<<(1+mpdhSize)<<endl;
                    fout<<"#latticeParam: 1"<<endl;
                    fout<<"#grainRadius: "<<0.5*maxBin<<endl;
                    fout<<"#structure: hcp"<<endl;
                    fout<<"#sphere: 1"<<endl;
                    fout<<"#modified: 1"<<endl;
                    fout<<"#atoms: "; for(auto &aname : atomTypes) fout<<" "<<aname; fout<<endl;
                    fout<<"#numOfatoms: "<<numOfAtoms<<endl;
                    fout<<"#binWidth: "<<setw(colwh)<<setprecision(nprec)<<1.0/binWidth<<endl;
                    fout<<"#comment: WARNING: number of atoms and radius is approximated"<<endl;

                    /// wypisuje nazwy column
                    fout<<"#"<<setw(colwh)<<'X';

                    for (i=0;i<numOfAtomTypes;i++){
                        for(j=i;j<numOfAtomTypes;j++)
                                fout<<" "<<setw(colwh)<<std::string(atomTypes[i]+"-"+atomTypes[j]);
                    }

                    /// koniec naglowka zaznaczony #
                    fout<<endl;


                    if(numberOfAtomTypes==1){
                    double x;
                    auto  &sumBinij=dataYnn[0];

                        for(i=0;i<dataSize;i++){
                            if(sumBinij[i]>0){
                                x=minBin+i*iBinWidth;
                                fout<<" "<<setprecision(nprec)<<setw(colwh)<<x<<" "<<setprecision(nprec)<<setw(colwh)<<sumBinij[i]<<endl;
                                numOfNonZeroBins++;
                            }
                        }
                    }
                    else{
                    double x;
                    const size_t numOfBinTypes=dataYnn.size();

                            for(i=0;i<dataSize;i++){
                                x=minBin+i*iBinWidth;
                                fout<<" "<<setprecision(nprec+2)<<setw(colwh+2)<<x;

                                for(j=0;j<numOfBinTypes;j++)
                                    fout<<" "<<setprecision(nprec)<<setw(colwh)<<dataYnn[j][i];

                                numOfNonZeroBins++;

                                fout<<endl;
                            }

                    }//else's end

                    fout.seekp(strpos);
                    fout<<numOfNonZeroBins<<"\t";

                    fout.close();


}
//===================================================================================
void Cavepdh::printPrm()
{
                cout<<"#openFiles: "<<openFiles<<endl;
                cout<<"#weight: "<<weight<<endl;
                cout<<"#numberOfAtomTypes: "<<numberOfAtomTypes<<endl;
                cout<<"#binWidth: "<< binWidth<<endl;
                cout<<"#minmaxBin: "<<minBin<<" "<<maxBin<<endl;


               // cout<<"#filesNumber: "<<fileNames.size()<<endl;
                cout<<"#filesListing:"<<endl;

             //   for(auto &file: randFileNames)
               //     cout<<"\t"<<file<<endl;

}
//===================================================================================
void Cavepdh::exportToCpdh()
{

                if(numberOfAtomTypes==1){
                const double iBinWidth=1.0/binWidth;
                auto  &sumBinij=dataYnn[0];

                        pdh->grain->resetPrms();
                        pdh->clearData();

                        pdh->dataX.allocMem(sumBinij.size);

                        for(unsigned i=0;i<sumBinij.size;i++)
                            pdh->dataX[i]=minBin+i*iBinWidth;

                        pdh->dataYii=std::move(sumBinij);

                        pdh->grain->atomNamesNumber.resize(1);
                        pdh->grain->atomNamesNumber[0]=numOfAtoms;

                        pdh->grain->atomTypes.resize(1);
                        pdh->grain->atomTypes[0].name=atomTypes[0];
                }

}

//===================================================================================
bool Cavepdh::parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__)
{
 const std::string send("end");


            clearData();
            ptr_uvar=uvar__;

            while(cmd[index]!=send){


                if(cmd[index]=="files2ave"){
                    numberOfFiles2Ave=cmd[index++][1];
                continue;
                }


                if(cmd[index]=="open"){
                    openFiles=cmd[index++][1];
                    Script::replaceVars(ptr_uvar,openFiles);
                continue;
                }


                if(cmd[index]=="save"){
                    fileNameOut=cmd[index++][1];
                    Script::replaceVars(ptr_uvar,fileNameOut);
                continue;
                }

                if(cmd[index]=="weight"){
                    weight=cmd[index][1];
                    Script::replaceVars(ptr_uvar,weight);
                    index++;
                continue;
                }

                if(cmd[index]=="printprm"){
                    printprm=true;
                    index++;
                continue;
                }

                cerr<<"Error: cavepdh.cpp  unknown command "<< cmd[index]<<" line "<<__LINE__<<endl;
            return false;
            }

            try {

                if(openFiles.empty())  throw Cavepdh::ERR_OPEN;

            return true;

            }
            catch(EAvePdhStatus status){

                switch (status) {
                case ERR_OPEN : cerr<<"ERROR: empty set, did you forgot set 'open' command"; break;
                default: cerr<<"ERROR: Unknown "<<__FILE__<<":"<<__LINE__<<endl;

                }

                cout<<endl;


            }
            catch(...){
                cerr<<"ERROR: unknown type "<<__FILE__<<":"<<__LINE__<<endl;
            }



return false;
}
//===================================================================================


/*
 *
 *
#!/bin/bash

flist=$(ls pdh_1??.pdhs)
echo $flist | wc -w >pdhfile.tmp
echo $flist | xargs head -n 12 | grep -v "#" | egrep "^[[:space:]]*[0-9]" | sort -n -k1 | head -n 1 >>pdhfile.tmp
echo $flist | xargs tail -n 1 | egrep "^[[:space:]]*[0-9]+" | sort -n -k1 | tail -n 1 >>pdhfile.tmp
echo $flist >>pdhfile.tmp



 *
 *
 *
 */


//===================================================================================
void Cavepdh::calc()
{
                    status=Cavepdh::RUN;

const std::string tmpPdhFiles("pdhfiles.tmp");
const string lsFiles("rm -f "+tmpPdhFiles+"; "+
                     "flist=$("+openFiles+");"+
                     "echo $flist|wc -w>"+tmpPdhFiles+";"+
                     "echo $flist|xargs head -n 12 |grep -v \"#\"|egrep \"^[[:space:]]*[0-9]\"|sed 's/^\\s*//'|sort -n -k1 |head -n 1 >>"+tmpPdhFiles+";"+
                     "echo $flist|xargs tail -n 1 |egrep \"[[:space:]]*[0-9]+\"|sed 's/^\\s*//'|sort -n -k1 |tail -n 1 >>"+tmpPdhFiles+";"+
                     "echo $flist|xargs -n 1 >>"+tmpPdhFiles
                     );

            //if(DB) { cout<<"lsFiles "<<lsFiles<<":"<<lsFiles.length()<<endl; }
            
            //if(DB) {    cout<<std::system(lsFiles.c_str())<<endl;}
            //else std::system(lsFiles.c_str())<<endl;

            std::system(lsFiles.c_str());

fstream fileList(tmpPdhFiles,ios::in);

            try{
                    //................................................................

                    if(!fileList){
                        cerr<<"ERROR: files listing failure"<<endl;
                    throw 0;
                    }

                    //..................... READ HEADER ...................................
            size_t numOfFiles;

                     fileList>>numOfFiles;

                     if(!fileList){
                         cerr<<"ERROR: number of files is equal  0"<<endl;
                     throw 1;
                     }

                     if(!numberOfFiles2Ave.empty()){
                     csize numOfFiles2Ave=std::stoi(numberOfFiles2Ave);

                            if(numOfFiles<numOfFiles2Ave){
                                cerr<<"ERROR: number of listed files less than 'files2ave'"<<endl;
                            throw 1;
                            }

                            numOfFiles=numOfFiles2Ave;
                     }


                    fileList>>maxBin>>minBin;

                    if(fileNameOut.find(".lhs")!=string::npos){
                        minBin=0;
                    }


                    if(DB) cout<<" min, max: "<<minBin<<" , "<<maxBin<<endl;

                    while(fileList.get()!='\n' && !fileList.eof() ) ;

                    if(!fileList.good()){
                        status=Cavepdh::ERR_LIST;
                        cerr<<"ERROR: fileList corrupted"<<endl;
                    return;
                    }

                     //......................... PREPARE STATISTIC .................................


            string fileName;
            vector<string> fileNames;

                     //cout<<"prob "<<prob<<", (%)="<<100.0*prob/numOfFiles<<endl;

                    fileNames.reserve(numOfFiles);

                    //........................ READ FILENAMES OF RANDOM SELECTED FILES ..........................

                    for(size_t i=0;i<numOfFiles;i++){
                        fileList>>fileName;
                        if(DB)cout<<fileName<<endl;
                        fileNames.push_back(fileName);
                    }

                    fileList.close();

                    fileNames.shrink_to_fit();

                    //........................... READ HEADER OF A FIRST FILE ........................


                    fstream file0(fileNames.front(),ios::in);

                            if(!file0){ 
                                cerr<<"ERROR: couldn't open a file "<<fileNames.front()<<endl;
                                cerr<<"       see: "<<__FILE__<<":"<<__LINE__<<endl;
                                throw 2;
                            }

                    string cmd;
                    size_t rows=0,cols=0;


                            binWidth=0;
                            numOfAtoms=0;

                            while(file0.peek()=='#' && !file0.eof() ){
                                file0>>cmd;

                                if(cmd.find("sizeRC")!=string::npos){
                                    file0>>rows>>cols;

                                    if(! rows*cols) {cerr<<"ERROR: rows*cols==0"<<endl; throw 0;}
                                }

                                if(cmd.find("atomTypes")!=string::npos){
                                    file0>>numberOfAtomTypes;

                                    atomTypes.resize(numberOfAtomTypes);
                                    for(string &aname: atomTypes)
                                        file0>>aname;

                                }


                                if(cmd.find("binWidth")!=string::npos){
                                    file0>>binWidth;

                                    if(binWidth<1) binWidth=1.0/binWidth;

                                    break;
                                }

                                /// go to the EOL
                                while(file0.get()!='\n' && !file0.eof()) ;

                            }


                            file0.close();
                            if(binWidth==0){
                                cerr<<"ERROR: wrong format (binwidth is equal 0 or not set)"<<endl;
                            return;
                            }
                            if(numberOfAtomTypes==0){
                                cerr<<"ERROR: wrong format (binwidth is equal 0 or not set)"<<endl;
                            return;

                            }

                   //................... PREPARE BUFFER FOR DATA .......................
                    const size_t binSize=static_cast<size_t>(std::ceil(maxBin-minBin+1)*binWidth);
                    double x,y,totalSumBins;
                    size_t numOfRows;
                    unsigned binIndex;
                    cdouble inumOfRandFiles=1.0/fileNames.size();
                    const size_t numOfBinTypes=cols-1;

                            dataYnn.resize(numOfBinTypes);


                    if(numberOfAtomTypes==1){
                    CProgress progress;
                    auto  &sumBinij=dataYnn[0];
                            sumBinij.allocMem(binSize,0);
                            totalSumBins=0;

                            progress.title=std::string(" ave pdh ");
                            progress.start(fileNames.size());

                            //................. READ DATA FROM FILES AND FILL UP THE BUFFER ...............

                            for(string &fname : fileNames){
                            fstream fin(fname,ios::in);

                                numOfRows=0;

                                //........... READ HEADER .................
                                if(!fin){
                                    cerr<<"ERROR: couldn't open a file "<<fname<<endl;
                                    cerr<<"       see: "<<__FILE__<<":"<<__LINE__<<endl;
                                return;
                                }

                                fin.exceptions(ios::failbit | ios::badbit | ios::eofbit);

                                while(fin.peek()=='#' && !fin.eof()){

                                    fin>>cmd;


                                    if(cmd.find("sizeRC")!=string::npos){
                                        fin>>numOfRows;
                                    }

                                    while(fin.get()!='\n' && !fin.eof()) ;

                                }

                                if(!fin.good() || fin.eof()){
                                    cerr<<"ERROR: failure during reading the file: "<<fname<<endl;
                                    cerr<<"     good(), eof(): "<<fin.good()<<", "<<fin.eof()<<endl;
                                return;
                                }

                                if(!numOfRows){
                                    cerr<<"ERROR: unknown number of rows, file: "<<fname<<endl;
                                return;
                                }

                                //.......... READ DATA ........................................

                                if(DB){cout<<__FILE__<<":"<<__LINE__<<fname<<endl;}

                                for(size_t j=0;j<numOfRows;j++){
                                    fin>>x>>y;

                                    binIndex=static_cast<unsigned>(binWidth*(x-minBin));

                                    if(binIndex<binSize){
                                        sumBinij[binIndex]+=y*inumOfRandFiles;
                                        totalSumBins+=y;
                                    }
                                    else {
                                        cerr<<"ERROR: bin index out of range "<<binIndex<<", maxBin: "<<binSize<<endl;
                                        cerr<<"   x "<<x<<endl;
                                    }
                                }

                                fin.close();
                                progress++;
                            }
                    }/// end : numberOfatoms==1
                    else {
                    CProgress progress;

                        progress.title=std::string(" ave pdh ");
                        progress.start(fileNames.size());


                        for(auto &bins: dataYnn){
                            bins.allocMem(binSize,0);
                        }

                        totalSumBins=0;

                        for(string &fname : fileNames){
                        fstream fin(fname,ios::in);

                                numOfRows=0;

                                //........... READ HEADER .................
                                if(!fin){
                                    cerr<<"ERROR: couldn't open a file "<<fname<<endl;
                                    cerr<<"       see: "<<__FILE__<<":"<<__LINE__<<endl;
                                return;
                                }

                                fin.exceptions(ios::failbit | ios::badbit | ios::eofbit);

                                while(fin.peek()=='#' && !fin.eof()){

                                    fin>>cmd;


                                    if(cmd.find("sizeRC")!=string::npos){
                                        fin>>numOfRows;
                                    }

                                    while(fin.get()!='\n' && !fin.eof()) ;

                                }

                                if(!fin.good() || fin.eof()){
                                    cerr<<"ERROR: failure during reading the file: "<<fname<<endl;
                                    cerr<<"     good(), eof(): "<<fin.good()<<", "<<fin.eof()<<endl;
                                return;
                                }

                                if(!numOfRows){
                                    cerr<<"ERROR: unknown number of rows, file: "<<fname<<endl;
                                return;
                                }

                                //.......... READ DATA ........................................

                                if(DB){cout<<__FILE__<<":"<<__LINE__<<"   "<<fname<<endl;}

                                for(size_t j=0;j<numOfRows;j++){
                                    fin>>x;

                                    binIndex=static_cast<unsigned>(binWidth*(x-minBin));

                                    if(binIndex<binSize){

                                        for(size_t k=0;k<numOfBinTypes;k++){
                                            fin>>y;
                                            dataYnn[k][binIndex]+=y*inumOfRandFiles;
                                            totalSumBins+=y;
                                        }
                                    }
                                    else {
                                        cerr<<"ERROR: bin index out of range "<<binIndex<<", maxBin: "<<binSize<<endl;
                                        cerr<<"   x "<<x<<endl;

                                        while(fin.get()!='\n' && !fin.eof())
                                            ;

                                        if(fin.eof()) { cerr<<"ERROR: unexptected end of file"<<endl; throw 0;}
                                    }

                                }

                                fin.close();
                                progress++;
                            }

                    }

                    totalSumBins*=inumOfRandFiles;

                    //  WARNING:  the numOfAtoms is an approximation
                    numOfAtoms= static_cast<size_t>(  std::sqrt(totalSumBins)*M_SQRT2);


                    if(!fileNameOut.empty()){

                        if(fileNameOut.find(".pdhs")!=string::npos) saveFilePdhs();
                        else
                            if(fileNameOut.find(".lhs")!=string::npos)
                                saveFileLhs();
                            else
                                saveFileDat();

                    }

                    exportToCpdh();

                    if(printprm) printPrm();

                    status=Cavepdh::OK;
            }
            catch(std::ifstream::failure &e ){
                cerr<<" ERROR "<<e.what()<<endl;
                status=Cavepdh::ERR_OPEN;

            }
            catch(...){
                status=Cavepdh::ERR_OPEN;
                fileList.close();
            }


}
//===================================================================================
//===================================================================================
void Cavepdh::clearData()
{
        openFiles.clear();
        fileNameOut.clear();
        //weight.clear();
        dataYnn.clear();
        atomTypes.clear();
        //sumBinij.freeMem();
        //randFileNames.clear();
        numberOfFiles2Ave.clear();
        printprm=false;
}

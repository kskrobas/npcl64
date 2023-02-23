/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cmergepdh.cpp
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


#include "cmergepdh.h"
#include "scriptanalyser.h"
#include "crandom.h"

#include<string>
#include<fstream>
//#include <filesystem>
#include <bits/stdc++.h>



//namespace fs = std::filesystem;


#if DEBUG
#define DB 1
#else
#define	DB 0
#endif



double StNameRadii::maxRadii=0;


Cmergepdh::Cmergepdh()
{

}

Cmergepdh::~Cmergepdh()
{

}
//--------------------------//--------------------------//--------------------------
bool Cmergepdh::parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__)
{

const std::string send("end");


           clearData();
           ptr_uvar=uvar__;

            while(cmd[index]!=send){

                if(cmd[index]=="flimit"){
                    flimit=cmd[index++][1];
                continue;
                }


                if(cmd[index]=="open"){
                   openFiles=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,openFiles);
                continue;
                }


                if(cmd[index]=="range"){
                   from=cmd[index][1];
                   step=cmd[index][2];
                     to=cmd[index][3];

                    index++;
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


               cerr<<"Error: cmergepdh.cpp  unknown command "<< cmd[index]<<" line "<<__LINE__<<endl;
           return false;
           }



           try {

               if(openFiles.empty())  throw Cmergepdh::ERR_OPEN;
               if(weight.empty())   throw Cmergepdh::ERR_WEIGHT;
               if(from.empty())     throw Cmergepdh::ERR_RANGE;

           return true;

           }
           catch(EmergePdhStatus status){

                switch (status) {
                case ERR_OPEN : cerr<<"ERROR: empty set, did you forgot set 'open' command"; break;
                case ERR_WEIGHT : cerr<<"ERROR: weight distribution not given"; break;
                case ERR_RANGE : cerr<<"ERROR: range not given"; break;
                case OK:
                case RUN:
                        break;

                }

               cout<<endl;

           }
           catch(...){
               cerr<<"ERROR: unknown type "<<__FILE__<<":"<<__LINE__<<endl;
           }



return true;
}
//--------------------------//--------------------------//--------------------------

bool getPdhFiles(const string & openFiles,vector<StNameRadii> &pdhfiles)
{
const std::string sysCmdRmPipe("rm -f npclfifo*");
const std::string sysCmd(sysCmdRmPipe+"; mkfifo npclfifoA; "+openFiles+
                             " | wc -l >npclfifoA & mkfifo npclfifoB; "+openFiles+">npclfifoB &") ;
                    system(sysCmd.c_str());
string buff;
fstream npclFifo;
int numOfPdhFiles=0;

           //----------------------------//----------------------------
           try {

                  npclFifo.open("npclfifoA",ios::in);

                  if(npclFifo)
                      std::getline(npclFifo,buff,'\n');
                  else {
                      cerr<<"ERROR: "<<__FILE__<<":"<<__LINE__<<endl;
                      return false;
                  }

                  npclFifo.close();

                   numOfPdhFiles=std::stoi(buff);

                    cout<<"numOfpdhfiles:  "<<numOfPdhFiles<<endl;

                   if(!numOfPdhFiles){
                       cerr<<"ERROR: file listing is empty"<<endl;
                   return false;
                   }


                   pdhfiles.reserve(numOfPdhFiles);


                   npclFifo.open("npclfifoB",ios::in);

                   //cout<<__LINE__<<endl;

                   for(int i=0;i<numOfPdhFiles;i++){
                      std::getline(npclFifo,buff,'\n');
                      pdhfiles.push_back(StNameRadii(buff));
                   }

                   pdhfiles.shrink_to_fit();
                   npclFifo.close();

                   system(sysCmdRmPipe.c_str());

           }
           catch (int e) {
                   switch(e){
                   case 0: cerr<<"ERROR: number of listed files is equal 0";break;

                   }

                   cerr<<endl;

                   npclFifo.close();
                   system(sysCmdRmPipe.c_str());

           return false;
           }
           catch (...) {
                   cerr<<" ERROR: Unknown type "<<__FILE__<<":"<<__LINE__<<endl;

                   npclFifo.close();
                   system(sysCmdRmPipe.c_str());
           return false;
           }

return true;
}


//--------------------------//--------------------------//--------------------------
void Cmergepdh::calc()
{

vector < vector<StNameRadii> > vcells;
bool modePdhs=true;
                    status=Cmergepdh::RUN;


                    //--------------------------------------------------------

                    if(pdhfiles.empty()){
                        if(!getPdhFiles(openFiles,pdhfiles))
                            return;
                    }
                    else {
                        cout<<"WARNING: pdh files taken from the previous run"<<endl;
                    }

                    //--------------------------------------------------------

                    cout<<" data loading"<<endl;

                    if(openFiles.find(".pdhs")!=string::npos){
                        if(!getPdhsStatistic(vcells,pdhfiles) ) return ;
                    }
                    else {
                        if(!getPlsStatistic(vcells,pdhfiles) )  return ;

                        modePdhs=false;
                    }

                    cout<<"\r data loaded"<<endl;



                    if(DB) cout<<" total : minBin, maxBin "<<totalMinBin<<", "<<totalMaxBin<<endl;

                    //--------------------------------------------------------



csize numberOfDataBins=static_cast<size_t>(1+std::ceil((totalMaxBin-totalMinBin)*( (binWidth>1) ? binWidth : 1.0/binWidth)));
vector<string> weightPrm(split<string>(weight," "));
cdouble a=std::stod(weightPrm[1]);  ///  μ or m or from
cdouble b=std::stod(weightPrm[2]);  ///  σ or s or to
const double amp=std::stod(flimit);
double (*fprob)(cdouble &, cdouble &, cdouble &);
const double  minR=std::stod(from);
const double stepR=std::stod(step);
double r=minR+stepR*0.5;
double k0;
size_t numOfRandFiles,iter=1;
vector<int> randIdFile;
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

typedef void (Cmergepdh::*ptrMergeFun)(const vector<StNameRadii> &,const vector<int> & );
ptrMergeFun mergeFun=(modePdhs) ? &Cmergepdh::mergePdhsBlockPdhs : &Cmergepdh::mergePdhsBlockPls;


                    sumBinij.allocMem(numberOfDataBins);

                    #pragma omp parallel for
                    for(size_t i=0;i<sumBinij.size;i++)
                        sumBinij[i]=0;


                    if(weightPrm[0]=="uniform"){
                        k0=amp;
                        fprob=&uniform;
                    }
                    else {
                        if(weightPrm[0]=="normal"){
                        cdouble &rpeak=a;
                            fprob=&normal;
                            k0=1.0*amp/fprob(rpeak,a,b);
                        }
                        else{
                        cdouble rpeak=std::exp(a-b*b);
                            fprob=&lognormal;
                            k0=1.0*amp/fprob(rpeak,a,b);
                        }
                    }

                    for(auto &cell: vcells){

                        cout<<iter++<<"/"<<vcells.size()<<"  "<<setw(6)<<setprecision(7)<<r;
                        cout.flush();

                        if(cell.empty()) {cout<<" ignored "<<endl; r+=stepR; continue;}

                    std::uniform_int_distribution<int> intDistr(0,cell.size()-1);

                        randIdFile.clear();
                        numOfRandFiles=static_cast<size_t> (k0*fprob(r,a,b));

                        if(DB) cout<<"bin (r,height,prob) : "<<setw(6)<<setprecision(7)<<r<<", "<<vcells.size()<<", "<<numOfRandFiles<<endl;

                        if(numOfRandFiles<1){
                            cerr<<"WARNING: parameters of random distribution for r="<<r<<" give number of random files less than 0"<<endl;
                        continue;
                        }

                        randIdFile.resize(numOfRandFiles);

                        for(size_t i=0;i<numOfRandFiles;i++)
                            randIdFile[i]=intDistr(generator);

                        if(DB) cout<<"next block "<<cell[0].radii<<endl;

                        (this->*mergeFun)(cell,randIdFile);

                        cout<<"   numOfRandFiles totalNumberOfAtoms: "<<numOfRandFiles<<" "<<totalNumberOfAtoms<<endl;

                        r+=stepR;
                    }


            
                    status=Cmergepdh::OK;

                    if(!fileNameOut.empty()) saveFilePdhs();

                    exportToCpdh();
}



//--------------------------//--------------------------//--------------------------
bool Cmergepdh::getPdhsStatistic(vector<vector<StNameRadii>> &vcells, vector<StNameRadii> &pdhfiles)
{
fstream fin;
const size_t numOfPdhFiles=pdhfiles.size();
const size_t aveSizeOfBin=pdhfiles.size();
size_t iter=0,bin,cignore;
int rows=0;
double val;
string buff;
const double  minR=std::stod(from);
const double stepR=std::stod(step);
const double  maxR=std::stod(to);
const size_t numberOfCells=static_cast<size_t>(std::ceil(maxR-minR)/stepR); /// WARNING: possible round off errors
bool binMinMaxTag;


                vcells.resize(numberOfCells);

                for(auto &vcell: vcells)
                    vcell.reserve(aveSizeOfBin);

                totalMinBin=std::numeric_limits<double>::infinity();
                totalMaxBin=0;
                binWidth=0;
                numberOfAtomTypes=1; //more types to be implemented

                for(auto &pdhfile: pdhfiles){
                fstream fin(pdhfile.name,ios::in);

                        cout<<"\r"<<iter++<<"/"<<numOfPdhFiles;
                        cout.flush();

                        if(!fin){
                            cerr<<"error: "<<pdhfile.name<<" :"<<__FILE__<<":"<<__LINE__<<endl;
                            cerr.flush();

                        break;
                        }

                        cignore=0;
                        binMinMaxTag=false;




                        while(fin.peek()=='#' && !fin.eof()){
                            std::getline(fin,buff,'\n');


                            if(cignore==4) continue;



                            if(buff.find("#atomsNumber")!=string::npos){
                            const string val(buff.substr(13));
                            const size_t N=std::stoi(val);

                                    pdhfile.num2radii(N);

                                    if(pdhfile.radii>=minR && pdhfile.radii<maxR){
                                        bin=(pdhfile.radii-minR)/stepR;
                                        vcells[bin].emplace_back(pdhfile);
                                    }

                                    cignore++;
                            continue;
                            }

                            if(buff.find("#sizeRC")!=string::npos){
                            vector<string> tokens(split<string>(buff," "));
                                    rows=std::stoi(tokens[1]);

                                    cignore++;
                            continue;
                            }


                            if(buff.find("#binWidth")!=string::npos){
                            vector<string> tok(split<string>(buff," "));
                                    cignore++;

                                if(binWidth==0)
                                    binWidth=std::stod(tok[1]);
                                else{
                                    if(binWidth!=std::stod(tok[1])){
                                        cerr<<"ERROR: inconsistent binWidth, file: "<<pdhfile.name<<endl;
                                    return false;
                                    }

                                }
                            continue;
                            }

                            /// pdhs  version: 1
                            if(buff.find("#binMinMax")!=string::npos){
                            vector<string> tok(split<string>(buff," \t"));
                                    cignore++;
                                    binMinMaxTag=true;

                            const double binMin=std::stod(tok[1]);
                            const double binMax=std::stod(tok[2]);

                                    if(binMin<totalMinBin) totalMinBin=binMin;
                                    if(binMax>totalMaxBin) totalMaxBin=binMax;

                            }

                        }

                        if(fin.eof()){
                            cerr<<"ERROR: wrong format of file: "<<pdhfile.name<<endl;
                            fin.close();
                        return false;
                        }

                        /// pdhs version 0
                        if(!binMinMaxTag){

                                //-------------- minBin -----------------------
                                fin>>val;

                                if(val<totalMinBin)
                                    totalMinBin=val;

                                //-------------- ignore rows --------------------

                                for(int i=rows;i!=1;i--)
                                    std::getline(fin,buff,'\n');

                                //-------------- maxBin -----------------------
                                fin>>val;

                                if(val>totalMaxBin)
                                    totalMaxBin=val;

                        }

                        //---------------------------------------------


                        fin.close();

                        if(!fin.good()){
                            cerr<<"ERROR: wrong format of file: "<<pdhfile.name<<endl;
                            fin.close();
                        return false;
                        }

                }

                if(binWidth<1) binWidth=1.0/binWidth;

                return true;
}
//--------------------------//--------------------------//--------------------------
bool Cmergepdh::getPlsStatistic(vector<vector<StNameRadii> > &vcells, vector<StNameRadii> &pdhfiles)
{
fstream fin;
const size_t numOfPdhFiles=pdhfiles.size();
const size_t aveSizeOfBin=pdhfiles.size();
size_t iter=0,bin;
string buff;
const double  minR=std::stod(from);
const double stepR=std::stod(step);
const double  maxR=std::stod(to);
const size_t numberOfCells=static_cast<size_t>(std::ceil(maxR-minR)/stepR); /// WARNING: possible round off errors

size_t numOfAtoms,numOfBins;
double vald;

                vcells.resize(numberOfCells);

                for(auto &cell: vcells)
                    cell.reserve(aveSizeOfBin);

                totalMinBin=std::numeric_limits<double>::infinity();
                totalMaxBin=0;
                binWidth=0;
                numberOfAtomTypes=1; //more types to be implemented

                for(auto &pdhfile: pdhfiles){
                fstream fin(pdhfile.name,ios::in);

                        cout<<"\r"<<iter++<<"/"<<numOfPdhFiles;
                        cout.flush();

                        if(!fin){
                            cerr<<"error: "<<pdhfile.name<<" :"<<__FILE__<<":"<<__LINE__<<endl;
                            cerr.flush();

                        break;
                        }

                        //**********************************
                        fin.seekg(20,ios_base::beg);
                        fin.read( reinterpret_cast<char *>(&numOfAtoms),sizeof(size_t));

                        pdhfile.num2radii(numOfAtoms);

                        if(pdhfile.radii>=minR && pdhfile.radii<maxR){
                            bin=(pdhfile.radii-minR)/stepR;
                            vcells[bin].emplace_back(pdhfile);
                        }
                        //**********************************


                        fin.read( reinterpret_cast<char *>(&binWidth) ,sizeof(double));
                        fin.read( reinterpret_cast<char *>(&numOfBins),sizeof(size_t));


                        fin.read( reinterpret_cast<char *>(&vald),sizeof(double));
                        if(vald<totalMinBin) totalMinBin=vald;

                        fin.read( reinterpret_cast<char *>(&vald),sizeof(double));
                        if(vald>totalMaxBin) totalMaxBin=vald;

                        fin.close();

                }


return true;
}
//--------------------------//--------------------------//--------------------------
void Cmergepdh::mergePdhsBlockPdhs(const vector<StNameRadii> &pdhsBlock, const vector<int> &id)
{
string cmd,buff;
double x,y;
size_t numOfAtoms,binIndex;
size_t numOfRows=0;

                try{

                for(size_t i: id){
                const string fileName(pdhsBlock[i].name);
                fstream fin(fileName,ios::in);


                        ///------------- read header ---------------
                        while (fin.peek()=='#' && !fin.eof() ){
                            fin>>cmd;

                            if(cmd.find("sizeRC")!=string::npos){
                                fin>>numOfRows;
                                while(fin.get()!='\n' && !fin.eof()) ;

                            }

                            if(cmd.find("atomsNumber")!=string::npos){
                                fin>>numOfAtoms;
                                totalNumberOfAtoms+=numOfAtoms;

                                while(fin.get()!='\n' && !fin.eof()) ;

                            break;
                            }


                            while(fin.get()!='\n' && !fin.eof()) ;
                        }

                        while (fin.peek()=='#' && !fin.eof())
                            std::getline(fin,buff,'\n');


                        ///-------------- read data ----------------
                        for(size_t i=0;i<numOfRows;i++){

                                fin>>x>>y;

                                binIndex=static_cast<size_t>(binWidth*(x-totalMinBin));

                                 sumBinij[binIndex]+=y;

                               /* if(binIndex<sumBinij.size){
                                    sumBinij[binIndex]+=y;


                                }
                                else {
                                    cerr<<"ERROR: bin index out of range "<<binIndex<<", maxBin: "<<sumBinij.size<<endl;
                                    cerr<<"   x "<<x<<endl;
                                }*/
                        }
                        ///-------------------------------------------
                        ///


                        fin.close();

                        if(!fin.good()) {
                            cout<<"error: file : "<<fileName<<endl;
                        throw 1;
                        }
                    }

                }

                catch (const std::out_of_range& e) {
                   std::cerr << "ERROR: Out of Range error: " << e.what() << endl;
                 }
                catch(...){

                    cerr<<"ERROR "<<__FILE__<<":"<<__LINE__<<endl;

                }



}
//--------------------------//--------------------------//--------------------------
void Cmergepdh::mergePdhsBlockPls(const vector<StNameRadii> &pdhsBlock, const vector<int> &id)
{
string cmd,buff;
double x,y,binWidth;
size_t numOfAtoms,binIndex;
size_t numOfBins=0;



                try{

                        for(size_t i: id){
                        const string fileName(pdhsBlock[i].name);
                        fstream fin(fileName,ios::in | ios::binary);


                                fin.seekg(20,ios_base::beg);
                                fin.read( reinterpret_cast<char *>(&numOfAtoms),sizeof(size_t));
                                fin.read( reinterpret_cast<char *>(&binWidth)  ,sizeof(double));
                                fin.read( reinterpret_cast<char *>(&numOfBins) ,sizeof(size_t));
                                fin.seekg(20,ios_base::cur);

                                totalNumberOfAtoms+=numOfAtoms;

                                for(size_t i=0;i<numOfBins;i++){

                                    fin.read( reinterpret_cast<char *>(&x),sizeof(double));
                                    fin.read( reinterpret_cast<char *>(&y),sizeof(double));

                                    binIndex=static_cast<size_t>(binWidth*(x-totalMinBin));

                                    #pragma message (" if program works flawlessly  switch off binIndex<sumbinij ")
                                    if(binIndex<sumBinij.size)
                                        sumBinij[binIndex]+=y;
                                    else {
                                        cout<<" it shouldn't happened: index out of range "<<binIndex<<">="<<sumBinij.size<<endl;
                                    }
                                }


                                if(!fin.good()) {
                                    cout<<"error: file : "<<fileName<<endl;
                                throw 1;
                                }
                        }

                }
                catch (const std::out_of_range& e) {
                   std::cerr << "ERROR: Out of Range error: " << e.what() << endl;
                 }
                catch(...){

                    cerr<<"ERROR "<<__FILE__<<":"<<__LINE__<<endl;

                }

}
//--------------------------//--------------------------//--------------------------
void Cmergepdh::saveFilePdhs()
{
fstream fout(fileNameOut,ios::out);

                if(!fout){
                    cerr<<"ERROR: couldn't save data to file: "<<fileNameOut<<endl;
                return;
                }


const size_t dataSize=sumBinij.size;
csize numOfAtomTypes= numberOfAtomTypes;

csize mpdhSize=(numOfAtomTypes+1)*numOfAtomTypes/2;
std::time_t datetime = std::time(nullptr);
size_t i,j;
csize nprec=13;
csize colwh=14;

const int numOfDigits=static_cast<int> (std::log10(dataSize)+2);
std::streampos strpos;
size_t numOfNonZeroBins=0;


                fout<<"#ver: 0"<<endl;
                fout<<"#title: avePdh"<<endl;
                fout<<"#sizeRC: "; strpos=fout.tellp();
                fout<<setw(numOfDigits)<<" "<<(1+mpdhSize)<<endl;
                fout<<"#date: "<<std::asctime(std::localtime(&datetime));
                fout<<"#atomTypes: "<<numOfAtomTypes<<"\tC";

                //for(i=0;i<numOfAtomTypes;i++)
                //    fout<<"\t"<<grain->atomTypes[i].name;

                fout<<endl;

                fout<<"#atomsNumber: "<<totalNumberOfAtoms<<endl;
                fout<<"#atomsXYZ: ignore"<<endl;
                fout<<"#binWidth: "<<binWidth<<endl;
                fout<<"#comment: WARNING: number of atoms is approximated; (minBin,maxBin)=( ";
                    fout<<setprecision(nprec)<<setw(colwh)<<totalMinBin<<",";
                    fout<<setprecision(nprec)<<setw(colwh)<<totalMaxBin<<")"<<endl;
                fout<<"#";

                /// wypisuje nazwy column
                fout<<setw(colwh)<<'X';

                for (i=0;i<numOfAtomTypes;i++){
                    for(j=i;j<numOfAtomTypes;j++)
                            fout<<" "<<setw(colwh)<<std::string("C-C");

                   //fout<<" "<<setw(colwh)<<std::string( (*ptrAtomNames)[i]+"-"+(*ptrAtomNames)[j]);
                }

                /// koniec naglowka zaznaczony #
                fout<<endl;
                fout.fill('#');
                fout<<setw((colwh+1)*(mpdhSize+1))<<'#'<<endl;///13=set(12)+space
                fout.fill(' ');


                if(numberOfAtomTypes==1){
                double x;
                const double iBinWidth=1.0/binWidth;

                    for(i=0;i<dataSize;i++){
                        x=totalMinBin+i*iBinWidth;
                        if(sumBinij[i]>0){
                            fout<<" "<<setprecision(nprec)<<setw(colwh)<<x<<" "<<setprecision(nprec)<<setw(colwh)<<sumBinij[i]<<endl;
                            numOfNonZeroBins++;
                        }
                    }
                }
                else{

                    #pragma message( "WARNING: cavepdh.cpp: multitype systems saving is not implemented" )

                }//else's end

                fout.seekp(strpos);
                fout<<numOfNonZeroBins<<"\t";

                fout.close();

}
//--------------------------//--------------------------//--------------------------
void Cmergepdh::exportToCpdh()
{

            if(numberOfAtomTypes==1){
            const double iBinWidth=1.0/binWidth;

                    pdh->grain->resetPrms();
                    pdh->clearData();

                    pdh->dataX.allocMem(sumBinij.size);

                    for(unsigned i=0;i<sumBinij.size;i++)
                        pdh->dataX[i]=totalMinBin+i*iBinWidth;

                    pdh->dataYii=std::move(sumBinij);

                    pdh->grain->atomNamesNumber.resize(1);
                    pdh->grain->atomNamesNumber[0]=totalNumberOfAtoms;

                    pdh->grain->atomTypes.resize(1);
                    pdh->grain->atomTypes[0].name="C";
            }

}


//--------------------------//--------------------------//--------------------------
void Cmergepdh::clearData()
{
            printprm=false;
            totalNumberOfAtoms=0;
            numberOfAtomTypes=0;
            sumBinij.freeMem();
            from.clear();
            step.clear();
            to.clear();
            weight.clear();
}
//--------------------------//--------------------------//--------------------------

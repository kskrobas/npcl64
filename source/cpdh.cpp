/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cpdh.cpp
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


#include "cpdh.h"
#include "cprogress.h"
#include "colormsg.h"
#include "createdir.h"

#include <omp.h>
#include <iomanip>
#include <chrono>
#include <thread>
#include <future>
#include <ctime>
#include <assert.h>



inline position sqr(const position &x){return x*x;}
inline double sqrd (const position &x)  {return x*x;}



typedef std::string str;

#ifdef __linux__

//#elif _WIN32 or _WIN64
#else
    # define M_PI		3.14159265358979323846	/* pi */
    # define M_1_PIl	0.318309886183790671537767526745028724L /* 1/pi */
#endif

#if DEBUG
#define DB 1
#else
#define	DB 0
#endif	
//=============================================================================




Cpdh::Cpdh()
{

}

//===================================================================================
position Cpdh::getBinWidth()
{
const position wbin=std::stod(bin);

return (wbin>1)? 1.0/wbin : wbin;
}
//===================================================================================

class CBinSelector{
private:

    //-----------------------------------------
    void binFast(NanoGrain::StAtom *patom1)
    {
    const size_t bin=static_cast<size_t> (sqrt(  sqrd(patom1->x-patom0->x)+ sqrd(patom1->y-patom0->y)+ sqrd(patom1->z-patom0->z)));

            bins[bin]++;
    }
    //-----------------------------------------
    void binFromTo(NanoGrain::StAtom *patom1)
    {
    const size_t bin=static_cast<size_t> (sqrt(  sqrd(patom1->x-patom0->x)+ sqrd(patom1->y-patom0->y)+ sqrd(patom1->z-patom0->z)));

            if(bin>=from && bin<=to)
                bins[bin]++;

    }
    //-----------------------------------------

    void (CBinSelector::*fbin)(NanoGrain::StAtom *);

public:
    NanoGrain::StAtom *patom0;
    position from,to;
    size_t  *bins;

    CBinSelector(){
        fbin=&CBinSelector::binFast;
    }

    enum EMODE {FAST,FROMTO};

    void setMODE(EMODE mode){ fbin=(mode==FAST) ? &CBinSelector::binFast  : &CBinSelector::binFromTo; }

    void operator()(NanoGrain::StAtom *patom){ return  (this->*fbin)(patom);}

};

//===================================================================================
void Cpdh::monoLatticePdh()
{
            if(DB) cout<<"monolattice PDH computations"<<endl;


position wbin,ibin;
size_t size;

            if(range.empty()){
                wbin=getBinWidth();
                ibin=1.0/wbin;
                size=static_cast<size_t>(std::ceil(grain->maxR*2.5*ibin));
            }
            else {
            vector<string> stNst(split<string>(range," "));//start N stop
                size=static_cast<size_t>(std::stoi(stNst[1]));
                wbin=(std::stod(stNst[2])-std::stod(stNst[0]))/(size-1);
                ibin=1.0/wbin;
            }



const size_t numOfatoms=grain->atoms.size();
const size_t numOfatomsM1=numOfatoms-1;
const size_t numOfatomsHalf=numOfatoms/2;

const int numOfthreads=stoi(threads);


            if(!size || !numOfatoms){
                cerr<<" Error: PDH not possible, reason: grain not initialized"<<endl;
                status=Cpdh::ERR_EMPTYGRAIN;
            return ;
            }

            dataX.allocMem(size);
            dataYii.allocMem(size);

            for(size_t i=0;i<size;i++){
                dataX[i]=i*wbin;
                dataYii[i]=0;
            }
			
			/// optimalization on
			for(size_t i=0;i<numOfatoms;i++){
                grain->atoms[i]*=ibin;
            }


CProgress progress;

			progress.title=(string(" PDH "));
			progress.start(numOfatomsM1);

            startTime = std::time(nullptr);
            status=Cpdh::RUN;                        

            omp_set_num_threads(numOfthreads);

            if(binMode==Cpdh::SINGLE){

                    #pragma omp parallel
                    {
                    size_t pi,pj,bin;
                    NanoGrain::StAtom *patom0,*patom1;
                    NanoGrain::StAtom *patoms=grain->atoms.data();
                    size_t  *p_dataYii=new size_t[size];
                    const int thNum{omp_get_thread_num()};

                            for(pi=0;pi<size;pi++)
                                p_dataYii[pi]=0;

                            #pragma omp for nowait
                            for(pi=0;pi<numOfatomsHalf;pi++){
                                patom0=patoms+pi;

                                for(pj=pi+1,patom1=patom0+1;pj<numOfatoms;pj++,patom1++){
                                    bin=static_cast<size_t> (sqrt(  sqrd(patom1->x-patom0->x)+
                                                                    sqrd(patom1->y-patom0->y)+
                                                                    sqrd(patom1->z-patom0->z)));

                                    if(DB){
                                        if(!bin) cerr<<"ERROR: bin=0"<<endl;
                                    }

                                    p_dataYii[bin]++;
                                }
                                #pragma omp critical
                                progress++;
                            }

                            #pragma omp for nowait
                            for(pi=numOfatomsM1-1;pi>=numOfatomsHalf;pi--){
                                patom0=patoms+pi;

                                for(pj=pi+1,patom1=patom0+1;pj<numOfatoms;pj++,patom1++){
                                    bin=static_cast<size_t> (sqrt(  sqrd(patom1->x-patom0->x)+
                                                                    sqrd(patom1->y-patom0->y)+
                                                                    sqrd(patom1->z-patom0->z)));

                                    if(DB){
                                        if(!bin) cerr<<"ERROR: bin=0"<<endl;
                                    }


                                    p_dataYii[bin]++;
                                }
                                #pragma omp critical
                                progress++;
                            }

                            #pragma omp critical
                            {
                                for(pi=0;pi<size;pi++)
                                    dataYii[pi]+=p_dataYii[pi];
                            }

                            delete [] p_dataYii;
                    }
            }
            else{
                    #pragma omp parallel
                    {
                    size_t pi,pj,bin;
                    NanoGrain::StAtom *patom0,*patom1;
                    NanoGrain::StAtom *patoms=grain->atoms.data();
                    size_t  *p_dataYii=new size_t[size];
                    double r,beta;

                            for(pi=0;pi<size;pi++)
                                p_dataYii[pi]=0;

                            #pragma omp for nowait
                            for(pi=0;pi<numOfatomsM1;pi++){
                                patom0=patoms+pi;
                                for(pj=pi+1,patom1=patom0+1;pj<numOfatoms;pj++,patom1++){

                                    r=(size_t)sqrt( sqr(patom1->x-patom0->x)+
                                                         sqr(patom1->y-patom0->y)+
                                                         sqr(patom1->z-patom0->z));
                                    bin=(size_t) r;
                                    beta=r-bin;
                                    p_dataYii[bin+1]+=beta;
                                    p_dataYii[bin]+=1-beta;

                                }
                                #pragma omp critical
                                progress++;
                            }

                            #pragma omp critical
                            {
                                for(pi=0;pi<size;pi++)
                                    dataYii[pi]+=p_dataYii[pi];
                            }

                            delete [] p_dataYii;
                   }//pragma end
              }// else end

            stopTime = std::time(nullptr);
            status=Cpdh::OK;
			
			/// optimalization off
			for(size_t i=0;i<numOfatoms;i++){
                grain->atoms[i]*=wbin;
            }
}
//===================================================================================

bool Cpdh::ijLatticePdh(StMultiPdhPrm *prm, vector<dataPdh>::iterator &buffPdh,const string &progressTitle)
{

//cpos ibin=prm->ibin;
bool err=false;
CProgress progress;

			progress.title=std::move(progressTitle);
			

            if(binMode==Cpdh::SINGLE){

                    #pragma omp parallel
                    {
                    size_t pi,pj,bin;
                    size_t numOfatomsI,numOfatomsJ;
                    NanoGrain::StAtom *patom0,*patom1;
                    NanoGrain::StAtom *patomsI;//=grain->atoms.data()+prm->ifrom;
                    NanoGrain::StAtom *patomsJ;//=grain->atoms.data()+prm->jfrom;
                    csize size=prm->size;
                    size_t  *p_dataYii=new size_t[size];

                            for(pi=0;pi<size;pi++)
                                p_dataYii[pi]=0;


                            if(prm->ifrom == prm->jfrom){ // monoatomic lattice : type I == type J

                                patomsI=grain->atoms.data()+prm->ifrom;
                                numOfatomsI=prm->ito-prm->ifrom-1;
                                numOfatomsJ=prm->ito-prm->ifrom;

                            const size_t numOfatomsHalf=numOfatomsI/2;

                                #pragma omp master
                                {
                                    progress.start(numOfatomsI);
                                }


                                //// first half
                                 #pragma omp for nowait
                                 for(pi=0;pi<numOfatomsHalf;pi++){
                                     patom0=patomsI+pi;
                                     for(pj=pi+1;pj<numOfatomsJ;pj++){
                                         patom1=patomsI+pj;

                                         bin=static_cast<size_t>(sqrt(   sqrd(patom1->x-patom0->x)+
                                                                         sqrd(patom1->y-patom0->y)+
                                                                         sqrd(patom1->z-patom0->z)) );

                                         p_dataYii[bin]++;
                                     }

                                     #pragma omp critical
                                     progress++;
                                 }//for outer


                                 //// second half
                                 #pragma omp for nowait
                                 for(pi=numOfatomsI-2;pi>=numOfatomsHalf;pi--){
                                     patom0=patomsI+pi;

                                     for(pj=pi+1;pj<numOfatomsJ;pj++){
                                         patom1=patomsI+pj;

                                         bin=static_cast<size_t>(sqrt(  sqrd(patom1->x-patom0->x)+
                                                                        sqrd(patom1->y-patom0->y)+
                                                                        sqrd(patom1->z-patom0->z)) );
                                         p_dataYii[bin]++;
                                     }
                                     #pragma omp critical
                                     progress++;
                                 }


                            }// 312 if
                            else{ // biatomic lattice
                                patomsI=grain->atoms.data()+prm->ifrom;
                                patomsJ=grain->atoms.data()+prm->jfrom;

                                numOfatomsI=prm->ito-prm->ifrom;
                                numOfatomsJ=prm->jto-prm->jfrom;

                            const size_t numOfatomsHalf=numOfatomsI/2;

                                #pragma omp master
                                {
                                    progress.start(numOfatomsI);
                                }

                                /// first half
                               #pragma omp for nowait
                               for(pi=0;pi<numOfatomsHalf;pi++){
                                   patom0=patomsI+pi;

                                   for(pj=0;pj<numOfatomsJ;pj++){
                                       patom1=patomsJ+pj;

                                       bin=(size_t)sqrt(  sqrd(patom1->x-patom0->x)+
                                                               sqrd(patom1->y-patom0->y)+
                                                               sqrd(patom1->z-patom0->z));

                                       p_dataYii[bin]++;
                                   }

                                   #pragma omp critical
                                   progress++;
                               }// for outer


                               #pragma omp for nowait
                               for(pi=numOfatomsI-1;pi>=numOfatomsHalf;pi--){
                                   patom0=patomsI+pi;

                                   for(pj=0;pj<numOfatomsJ;pj++){
                                       patom1=patomsJ+pj;

                                       bin=(size_t)sqrt(  sqrd(patom1->x-patom0->x)+
                                                               sqrd(patom1->y-patom0->y)+
                                                               sqrd(patom1->z-patom0->z));

                                       p_dataYii[bin]++;
                                   }

                                   #pragma omp critical
                                   progress++;
                               }// for outer
                            }//else

                            #pragma omp critical
                            {
                                for(pi=0;pi<size;pi++)
                                    (*buffPdh)[pi]+=p_dataYii[pi];
                            }

                            delete [] p_dataYii;
                    }
            }
           if((*buffPdh)[0]){
                cerr<<" ERROR:  detected atoms at the same positions "<<endl;
                err=true;
            }
			
			progress.stop();
			cout<<endl;
return err;
}
//===================================================================================
void Cpdh::multiLatticePdh()
{
using namespace NanoGrain;
csize numOfAtomTypes=grain->atomTypes.size();
csize numOfAtoms=grain->atoms.size();
csize mpdhSize=(numOfAtomTypes+1)*numOfAtomTypes/2;
const position wbin=getBinWidth();
const position ibin=1.0/wbin;
csize size=std::ceil(grain->maxR*2.25*ibin);
size_t iAtomType,jAtomType;
vector<size_t> atypesPos(numOfAtomTypes+1);
vector<size_t>::iterator iter_atp;
StMultiPdhPrm multiPdhPrm;
vector<dataPdh>::iterator iter_dataynn;
const int numOfthreads=stoi(threads);

            if( !size ){
                cerr<<"ERROR: pdh computation is not possible, grain->maxR*2.125*ibin=0"<<endl;
            return;
            }


            omp_set_num_threads(numOfthreads);

            //// multiPdh results, prepare all output data buffers

            dataX.allocMem(size);
            dataYnn.resize(mpdhSize);

            for(size_t i=0;i<mpdhSize;i++)
                dataYnn[i].allocMem(size);


            for(size_t i=0;i<size;i++){
                dataX[i]=i*wbin;
                for(size_t j=0;j<mpdhSize;j++)
                    dataYnn[j][i]=0;
            }

            iter_dataynn=dataYnn.begin();

          //  (*iter_dataynn)[0];


            //// prepare pdh common parameters

            multiPdhPrm.wbin=wbin;
            multiPdhPrm.ibin=ibin;
            multiPdhPrm.size=size;


            //// prepare positions of first atom for each type

            iAtomType=0;
            atypesPos[0]=0;

            iter_atp=atypesPos.begin()+1;

             // **first atom position in the array for each type
            for(size_t i=0;i<numOfAtoms;i++){
                if(grain->atoms[i].atype!=iAtomType){
                    iAtomType=grain->atoms[i].atype;
                    (*iter_atp++)=i;
                }

                /// for optimalization purposes , multiply original positions by ibin.
                /// Recreation of positions at the end of current procedure
                grain->atoms[i]*=ibin;
            }

                // ** the last atom position is the numberofAtoms
                // ** it makes easy the 'for ... for ...' iterations
            *iter_atp=numOfAtoms;

            ////  main PDH task
            ///
int progressMultiPDHValue=1;

            startTime = std::time(nullptr);
            status=Cpdh::RUN;

            //runPDH=true;//?

            for(iAtomType=0;iAtomType<numOfAtomTypes;iAtomType++){

                multiPdhPrm.ifrom=atypesPos[iAtomType];
                multiPdhPrm.ito=atypesPos[iAtomType+1];
				
                for(jAtomType=iAtomType;jAtomType<numOfAtomTypes;jAtomType++,progressMultiPDHValue++){

                    multiPdhPrm.jfrom=atypesPos[jAtomType];
                    multiPdhPrm.jto=atypesPos[jAtomType+1]; 
                                       
					std::string mpdhProgressTitle(" PDH "+std::to_string(progressMultiPDHValue)+"/"+std::to_string(mpdhSize)+" ");

                    if(ijLatticePdh(&multiPdhPrm,iter_dataynn,mpdhProgressTitle)){
                    StAtom *atomI=grain->atoms.data()+multiPdhPrm.ifrom;
                    StAtom *atomJ=grain->atoms.data()+multiPdhPrm.jfrom;

                        cerr<<"\t\t types: "<<grain->atomTypes[atomI->atype].name<<
                                    ", "    <<grain->atomTypes[atomJ->atype].name<<endl;

                    }

                    iter_dataynn++;
                }
            }


            // ** recreating original atom positions
            for(size_t i=0;i<numOfAtoms;i++){
                grain->atoms[i]*=wbin;
            }

        stopTime = std::time(nullptr);
        status=Cpdh::OK;
}
//===================================================================================
void Cpdh::saveResults()
{

        if(!createDirsIfDontExist(fileName)){
            errMsg("couldn't create nested directories for "+fileName);
        throw Cpdh::ERR_FILEOPEN;
        }

        if(fileName.rfind(".dat")!=string::npos){
            saveDatFile();
        return;
        }

        if(fileName.rfind(".lhs")!=string::npos){
            saveLhsFile();
        return;
        }

        if(fileName.rfind(".pdhl")!=string::npos){
            savePdhlFile();
        return;
        }

        if(fileName.rfind(".pdhs")!=string::npos){
            savePdhsFile();
        return;
        }

        if(fileName.rfind(".plb")!=string::npos){
            saveBinFile();
        return;
        }

        if(fileName.rfind(".pls")!=string::npos){
            saveBinFilePLS();
        return;
        }

        cerr<<"Warning:  PDH results not saved"<<endl;

}
//===================================================================================
void Cpdh::openFile()
{

        if(fileNameIn.find(".dat")!=string::npos){
            openPdhDat();
        return;
        }


        //// both types: pdhs,pdhl
        if(fileNameIn.find(".pdh")!=string::npos){
            openPdhls();
        return;
        }


        if(fileNameIn.find(".lhs")!=string::npos){
            openPdhLhs();
        return;
        }

        ///
        if(fileNameIn.find(".pls")!=string::npos){
            openBinFilePLS();
        return;
        }

        cerr<<"ERROR: unrecognized PDH file format"<<endl;

        status=Cpdh::ERR_FILEOPEN;


}
//===================================================================================
void Cpdh::saveDatFile()
{
fstream fout(fileName,ios::out);

        if(!fout){
            cerr<<"ERROR: couldn't open file for saving, PDH results will be lost"<<endl;
            status=Cpdh::ERR_FILEOPEN;
        return;
        }


const size_t pdhsize=dataSize();
size_t i;


            if(grain->atomTypes.size()==1){
                for(i=0;i<pdhsize;i++)
                    fout<<setprecision(9)<<setw(9)<<dataX[i]<<"\t"<<setprecision(12)<<setw(11)<<dataYii[i]<<endl;
            }
            else{

                    for(size_t i=0;i<pdhsize;i++){
                        fout<<setprecision(9)<<setw(9)<<dataX[i];

                        for(dataPdh & pdh : dataYnn)
                            fout<<"\t"<<setprecision(12)<<setw(11)<<pdh[i];

                        fout<<endl;
                    }
            }


        fout.close();
}
//===================================================================================
void Cpdh::savePdhlFile()
{
fstream fout(fileName,ios::out);

        if(!fout){
            cerr<<"ERROR: couldn't open file for saving, PDH results will be lost"<<endl;
            status=Cpdh::ERR_FILEOPEN;
        return;
        }


const size_t dataSize=dataX.size;
csize numOfAtomTypes= grain->atomTypes.size();
csize numOfAtoms=   grain->atoms.size();
csize mpdhSize=(numOfAtomTypes+1)*numOfAtomTypes/2;
std::time_t datetime = std::time(nullptr);
size_t i,j;

            fout<<"#ver: 0"<<endl;
            fout<<"#title: "<<title<<endl;
            fout<<"#sizeRC: "<<dataSize<<"\t"<<(1+mpdhSize)<<endl;
            fout<<"#date: "<<std::asctime(std::localtime(&datetime));
            fout<<"#atomTypes: "<<numOfAtomTypes<<"\t";

            for(i=0;i<numOfAtomTypes;i++)
                fout<<"\t"<<grain->atomTypes[i].name;

            fout<<endl;

            fout<<"#atomsNumber: "<<numOfAtoms<<endl;
            fout<<"#atomsXYZ: "<<((grain->fileNameIn.empty()) ? "npcl" : grain->fileNameIn)<<endl;
            fout<<"#binWidth: "<<bin<<endl;
            fout<<"#comment: "<<comment<<endl;
            fout<<"#";

csize nprec=11;
csize colwh=12;

            /// wypisuje nazwy column
            fout<<setw(colwh)<<'X';

            for (i=0;i<numOfAtomTypes;i++){
                for(j=i;j<numOfAtomTypes;j++)                    
                        fout<<" "<<setw(colwh)<<std::string(grain->atomTypes[i].name+"-"+grain->atomTypes[j].name);

               //fout<<" "<<setw(colwh)<<std::string( (*ptrAtomNames)[i]+"-"+(*ptrAtomNames)[j]);
            }

            /// koniec naglowka zaznaczony #
            fout<<endl;
            fout.fill('#');
            fout<<setw((colwh+1)*(mpdhSize+1))<<'#'<<endl;///13=set(12)+space
            fout.fill(' ');


            if(grain->atomTypes.size()==1){
                for(i=0;i<dataSize;i++)
                    fout<<" "<<setprecision(nprec)<<setw(colwh)<<dataX[i]<<" "<<setprecision(nprec)<<setw(colwh)<<dataYii[i]<<endl;
            }
            else{

                for(i=0;i<dataSize;i++){

                        fout<<" "<<setprecision(nprec)<<setw(colwh)<<dataX[i];

                        for(dataPdh & pdh : dataYnn)
                          fout<<" "<<setprecision(nprec)<<setw(colwh)<<pdh[i];

                        fout<<endl;
                }//for's end
            }//else's end

            fout.close();
}
//===================================================================================
void Cpdh::savePdhsFile()
{
fstream fout(fileName,ios::out);

            if(!fout){
                cerr<<"ERROR: couldn't open file for saving, PDH results will be lost"<<endl;
                status=Cpdh::ERR_FILEOPEN;
            return;
            }


const size_t dataSize=dataX.size;
csize numOfAtomTypes=grain->atomTypes.size();
csize numOfAtoms=grain->atoms.size();
csize mpdhSize=(numOfAtomTypes+1)*numOfAtomTypes/2;
std::time_t datetime = std::time(nullptr);
size_t i,j,nonZeroBins=0;
std::streampos sizeRCpos,binMinMaxpos;
csize nprec=18;
csize colwh=nprec+3;
double binMin,binMax;



            binMin=std::numeric_limits<double>::infinity();
            binMax=0;

            fout<<"#ver: 1"<<endl;
            fout<<"#title: results of PDH calculations"<<endl;
            fout<<"#sizeRC: ";
            sizeRCpos=fout.tellp();
            fout<<setw(20)<<' '<<endl;

            fout<<"#date: "<<std::asctime(std::localtime(&datetime));
            fout<<"#atomTypes: "<<numOfAtomTypes<<"\t";

            for(i=0;i<numOfAtomTypes;i++)
                fout<<"\t"<<grain->atomTypes[i].name;

            fout<<endl;

            fout<<"#atomsNumber: "<<numOfAtoms<<endl;
            fout<<"#atomsXYZ: "<<((grain->fileNameIn.empty()) ? "npcl" : grain->fileNameIn)<<endl;
            fout<<"#binWidth: "<<bin<<endl;
            fout<<"#binMinMax: ";
            binMinMaxpos=fout.tellp();
            fout<<setw(2*colwh+1)<<' '<<endl;
            fout<<"#comment: "<<comment<<endl;
            fout<<"#";



            /// wypisuje nazwy column
            fout<<setw(colwh)<<'X';

            for (i=0;i<numOfAtomTypes;i++){
                for(j=i;j<numOfAtomTypes;j++)
                        fout<<" "<<setw(colwh)<<std::string(grain->atomTypes[i].name+"-"+grain->atomTypes[j].name);
            }

            /// koniec naglowka zaznaczony #
            fout<<endl;
            fout.fill('#');
            fout<<setw((colwh+1)*(mpdhSize+1))<<'#'<<endl;///13=set(12)+space
            fout.fill(' ');


            if(grain->atomTypes.size()==1){
                for( i=0;i<dataSize;i++){
                    if(dataYii[i]>0){
                        fout<<" "<<setprecision(nprec)<<setw(colwh)<<dataX[i]<<" "<<setprecision(nprec)<<setw(colwh)<<dataYii[i]<<endl;
                        nonZeroBins++;

                        if(dataX[i]<binMin) binMin=dataX[i];

                        binMax=dataX[i];
                    }
                }
            }
            else{
            bool nzb;

                for(i=0;i<dataSize;i++){

                     nzb=false;
                     for(dataPdh & pdh : dataYnn)
                        if(pdh[i]>0){
                            nzb=true;
                            nonZeroBins++;
                            break;
                        }


                      if(nzb){
                        fout<<" "<<setprecision(nprec)<<setw(colwh)<<dataX[i];

                        for(dataPdh & pdh : dataYnn)
                          fout<<" "<<setprecision(nprec)<<setw(colwh)<<pdh[i];

                        fout<<endl;
                      }
                }
            }


            fout.seekp(sizeRCpos,ios_base::beg);
            fout<<nonZeroBins<<"\t"<<(1+mpdhSize);

            fout.seekp(binMinMaxpos,ios_base::beg);
            fout<<binMin<<"\t"<<binMax;

            fout.close();
}
//===================================================================================
void Cpdh::saveLhsFile()
{
fstream fout(fileName,ios::out);

        if(!fout){
            cerr<<"ERROR: couldn't open file for saving, PDH results will be lost"<<endl;
            status=Cpdh::ERR_FILEOPEN;
        return;
        }


const size_t dataSize=dataX.size;
csize numOfAtomTypes=grain->atomTypes.size();
csize numOfAtoms=grain->atoms.size();
csize mpdhSize=(numOfAtomTypes+1)*numOfAtomTypes/2;
std::time_t datetime = std::time(nullptr);
size_t i,j;


            fout<<"#: AVE pdh"<<endl;
            fout<<"#date: "<<std::asctime(std::localtime(&datetime));
            fout<<"#sizeRC: "<<dataSize<<"\t"<<(1+mpdhSize)<<endl;
            fout<<"#latticeParam: 1"<<endl;
            fout<<"#grainRadius: "<<dataX[dataSize-1]*0.5<<endl;
            fout<<"#structure: hcp"<<endl;
            fout<<"#sphere: 0"<<endl;
            fout<<"#modified: 1"<<endl;
            fout<<"#atoms:";   for(i=0;i<numOfAtomTypes;i++) fout<<"\t"<<grain->atomTypes[i].name; fout<<endl;
            fout<<"#numOfatoms: "<<numOfAtoms<<endl;
            fout<<"#binWidth: "<<bin<<endl;
            fout<<"#comment: "<<endl;
            fout<<"#X";
            for (i=0;i<numOfAtomTypes;i++){
                for(j=i;j<numOfAtomTypes;j++)
                        fout<<"\t"<<std::string(grain->atomTypes[i].name+"-"+grain->atomTypes[j].name);
            }
            fout<<endl;

            if(grain->atomTypes.size()==1){
                for(i=0;i<dataSize;i++)
                    fout<<" "<<dataX[i]<<" "<<dataYii[i]<<endl;
            }
            else{

                for(i=0;i<dataSize;i++){

                        fout<<" "<<dataX[i];

                        for(dataPdh & pdh : dataYnn)
                          fout<<" "<<pdh[i];

                        fout<<endl;
                }//for's end
            }//else's end

            fout.close();
}
//===================================================================================
void Cpdh::saveExtFile()
{

}
//===================================================================================
void Cpdh::saveBinFile()
{
fstream fout(fileName,ios::out | ios::binary);


         if(!fout){
            cerr<<"ERROR: couldn't open file for saving, PDH results will be lost"<<endl;
            status=Cpdh::ERR_FILEOPEN;
        return;
        }


const string version("v00");
const string data("DATA");
const size_t pdhsize=dataSize();
const size_t atypSize=grain->atomTypes.size();
size_t i;

        position *datx=dataX.values;
        position *daty=dataYii.values;
        //for(i=0;i<numOfAtomTypes;i++) fout<<"\t"<<grain->atomTypes[i].name;

            fout.write( (char *) version.c_str(), version.length()*sizeof(char));
            fout.write( (char *) &atypSize,sizeof(size_t));

            if(atypSize==1){
            const char *atype=grain->atomTypes[0].name.c_str();
            const size_t atypeLength=grain->atomTypes[0].name.length();

                fout.write( (char*) &atypeLength,sizeof(size_t));
                fout.write( atype, sizeof(char)*atypeLength);


                fout.write((char *) &pdhsize,sizeof(size_t));
                fout.write( (char *) data.c_str(),    data.length()*sizeof(char) );


                for(i=0;i<pdhsize;i++,datx++,daty++){
                    fout.write((char *) datx , sizeof(position));
                    fout.write((char *) daty , sizeof(position));
                }

            }


        fout.close();

}
//===================================================================================
void Cpdh::saveBinFilePLS()
{
fstream fout(fileName,ios::out | ios::binary);

             if(DB)cerr<<"\ndata writing "<<fileName<<endl;


             if(!fout){
                cerr<<"ERROR: couldn't open file for saving, PDH results will be lost"<<endl;
                status=Cpdh::ERR_FILEOPEN;
            return;
            }


const string version("v01");
const string data("DATA");
const size_t pdhsize=dataSize();
const size_t atypSize=grain->atomTypes.size();
const size_t numOfBinTypes=dataYnn.size();
const double dbin=std::stod(bin);
csize numOfAtoms=grain->atoms.size();
size_t numOfBins;
double binMin,binMax;
std::streampos sizeRCpos;


            numOfBins=0;
            binMin=INFINITY;
            binMax=0;



            fout.write(reinterpret_cast<const char*>(version.c_str()), version.length()*sizeof(char));
            fout.write(reinterpret_cast<const char*> (&atypSize),sizeof(size_t));  /// liczba typów atomów


            if(atypSize==1){
            const char *atype=grain->atomTypes[0].name.c_str();
            const size_t atypeLength=grain->atomTypes[0].name.length();
            position *datx=dataX.values;
            position *daty=dataYii.values;


                fout.write( atype, sizeof(char)*atypeLength);
                fout.write(reinterpret_cast<const char*> (&atypeLength),sizeof(size_t));


                fout.write(reinterpret_cast<const char*> (&numOfAtoms), sizeof(size_t));
                fout.write(reinterpret_cast<const char*> (&dbin)      , sizeof(double));

                sizeRCpos=fout.tellp();
                fout.write(reinterpret_cast<const char*> (&numOfBins) , sizeof(size_t));
                fout.write(reinterpret_cast<const char*> (&binMin)    , sizeof(double));
                fout.write(reinterpret_cast<const char*> (&binMax)    , sizeof(double));

                /// insert 'DATA' marker
                fout.write(reinterpret_cast<const char*> (data.c_str()), data.length()*sizeof(char) );


                for(size_t i=0;i<pdhsize;i++,datx++,daty++){

                    if(*daty>0){
                        fout.write(reinterpret_cast<const char*> (datx) , sizeof(position));
                        fout.write(reinterpret_cast<const char*> (daty) , sizeof(position));

                        numOfBins++;
                        if(*datx<binMin) binMin=*datx;
                        binMax=*datx;
                    }
                }

                fout.seekp(sizeRCpos,ios_base::beg);
                fout.write(reinterpret_cast<const char*> (&numOfBins),sizeof(size_t));
                fout.write(reinterpret_cast<const char*> (&binMin)  , sizeof(double));
                fout.write(reinterpret_cast<const char*> (&binMax)  , sizeof(double));
            }
            else {

                    for (auto &atype : grain->atomTypes){
                        size_t atypeLength=atype.name.length();
                        fout.write(reinterpret_cast<const char *>(&atypeLength),sizeof(size_t));
                        fout.write(reinterpret_cast<const char *>(atype.name.c_str()),sizeof(char)*atype.name.length());
                    }

                    fout.write(reinterpret_cast<const char*> (&numOfAtoms), sizeof(size_t));
                    fout.write(reinterpret_cast<const char*> (&dbin)      , sizeof(double));

                    sizeRCpos=fout.tellp();
                    fout.write(reinterpret_cast<const char*> (&numOfBins) , sizeof(size_t));
                    fout.write(reinterpret_cast<const char*> (&binMin)    , sizeof(double));
                    fout.write(reinterpret_cast<const char*> (&binMax)    , sizeof(double));

                    /// insert 'DATA' marker
                    fout.write(reinterpret_cast<const char*> (data.c_str()), data.length()*sizeof(char) );

                bool nzb;
                position *datax=dataX.values;


                    for(size_t i=0;i<pdhsize;i++,datax++){

                        nzb=false;
                        for(dataPdh & pdh : dataYnn)
                           if(pdh[i]>0){
                               nzb=true;
                               numOfBins++;
                            break;
                           }


                        if(nzb) {
                            fout.write(reinterpret_cast<const char*> (datax) , sizeof(position));

                            if(*datax<binMin) binMin=*datax;
                            binMax=*datax;

                            for(size_t j=0;j<numOfBinTypes;j++){
                                fout.write(reinterpret_cast<char *>(&dataYnn[j][i]),sizeof(position));
                            }
                        }
                    }

                    fout.seekp(sizeRCpos,ios_base::beg);
                    fout.write(reinterpret_cast<const char*> (&numOfBins),sizeof(size_t));
                    fout.write(reinterpret_cast<const char*> (&binMin)  , sizeof(double));
                    fout.write(reinterpret_cast<const char*> (&binMax)  , sizeof(double));

                    if(DB)cout<<"OK"<<endl;

            }


        fout.close();

}
//===================================================================================
void Cpdh::openPdhDat()
{

}
//===================================================================================
void Cpdh::openPdhls()
{
ifstream fin(fileNameIn);


            if(!fin){
                cerr<<"ERROR: file couldn't be opened"<<endl;
            return;
            }

const size_t lineSize=63;
char line[lineSize];
string cmd,fver;
size_t rows,cols,atomTypes,atomsNumber,headerLines=0;
vector<string> atomNames;
vector<size_t> &atomNamesNumber=grain->atomNamesNumber;

            atomNamesNumber.clear();

            while(fin.peek()=='#' && !fin.eof() ){

                fin>>cmd;
                if(DB) {cout<<cmd<<endl;}


                if(cmd.find("#ver:")!=string::npos) { fin>>fver;}
                if(cmd.find("#sizeRC:")!=string::npos){fin>>rows>>cols;}

                if(cmd.find("#atomTypes:")!=string::npos){
                string aname;

                    fin>>atomTypes;
                    atomNames.reserve(atomTypes);

                    for(size_t i=0;i<atomTypes;i++){
                        fin>>aname;

                        if(!fin.good()) break;

                        atomNames.emplace_back(aname);
                    }

                    if(atomNames.size()!=atomTypes){
                        cerr<<"ERROR: number of atom names is not consistent with number of loaded values"<<endl;
                    break;
                    }


                }


                if(cmd.find("#atomsNumber:")!=string::npos){
                    fin>>atomsNumber;
                    //cout<<"#atomsNumber: "<<atomsNumber<<endl;
                }



                if(cmd.find("#binWidth:")!=string::npos){fin>>bin;}

                //pobiera all do konca linii
                while(fin.get()!='\n');

                 headerLines++;
             }

                if(atomTypes==1) {
                    atomNamesNumber.push_back(atomsNumber);
                }
                else{
                    //  ..... to be done soon as possible ........

                    #pragma message( "WARNING: cpdh.cpp: reading multitypes systems is not implemented" )
                }


const size_t pdhcols=(atomTypes*(atomTypes+1)/2);


             if(!fin.good () ||
                !( fver=="0" && headerLines==11) ||
                atomNames.size()!=atomTypes||  pdhcols!=(cols-1)){
                cerr<<"ERROR: the header of input file is corrupted: "<< fileNameIn<<endl;
                fin.close();
             return;
             }


size_t i,j;

            try{

                dataX.allocMem(rows);
                //fin.exceptions(ios::badbit | ios::failbit );

                    if(cols==2){

                            dataYii.allocMem(rows);

                            for(i=0;i<rows;i++){
                                fin>>dataX[i]>>dataYii[i];

                                if(!fin.good () )
                                throw 1;

                                //consistency data checking; the end of line should be empty
                                fin.getline(line,lineSize);
                                if(fin.gcount()>1){
                                string str(line);
                                    if(!ltrim(str).empty())
                                        throw 2;
                                }


                            }
                    }
                    else{
                            dataYnn.resize(pdhcols);

                            for(i=0;i<pdhcols;i++)
                                dataYnn[i].allocMem(rows);

                            for(i=0;i<rows;i++){
                                fin>>dataX[i];

                                for(j=0;j<pdhcols;j++)
                                    fin>>dataYnn[j][i];


                                if(!fin.good () ){
                                    if( fin.eof() && i==rows-1)
                                        fin.clear();
                                    else
                                    throw 1;
                                 }

                                //consistency data checking; the end of line should be empty
                                fin.getline(line,lineSize);
                                if(fin.gcount()>1){
                                string str(line);
                                    if(!ltrim(str).empty())
                                        throw 2;
                                }
                            }
                    }

                    title="data taken from file: "+fileNameIn;
                    status=Cpdh::OK;

                    grain->atomTypes.clear();
                    grain->atomTypes.assign(atomNames.begin(),atomNames.end());


                    //try to guess the number of atoms for each type
                    if(atomNamesNumber.empty()){
                    size_t nn,N,sumTot;
                    double sum;

                            grain->atomNamesNumber.clear();
                            grain->atomNamesNumber.resize(cols);


                            for(i=0,nn=0,sumTot=0;i<atomTypes;i++){

                                sum=0;
                                for(j=0;j<rows;j++)
                                    sum+=dataYnn[nn][j];

                                N=0.5+sqrt(0.25+2*sum);
                                grain->atomNamesNumber[i]=N;
                                sumTot+=N;

                                nn+=atomTypes-i;
                            }

                            cout<<"cpdh.cpp: "<<__LINE__<<" N="<<N<<endl;


                            if(sumTot!=atomsNumber){
                                cerr<<"WARNING: discrepancy between calculated number of atoms for each type and number of atoms given in header ";
                                cerr<<(float)(sumTot)/(float)(atomsNumber)<<endl;

                            }
                    }



            }
            catch (std::exception e){

                    if(fin.eof()){
                        if(i<rows-1){
                            cerr<<"ERROR: wrong number of rows "<<endl;
                            dataX.freeMem();
                            status=Cpdh::ERR_FILEOPEN;
                         }
                         else //it's ok
                            title="data taken from file: "+fileNameIn;
                    }
                    else{
                        title="data are corrupted "+fileNameIn;
                        dataX.freeMem();
                        status=Cpdh::ERR_FILEOPEN;
                    }


            }
            catch (int e){

                    if(e==1 && i<rows-1){
                        cerr<<"ERROR: wrong number of rows "<<endl;
                            dataX.freeMem();
                            status=Cpdh::ERR_FILEOPEN;
                    }
                    else{
                        if(e==2)
                            cerr<<"ERROR: number of columns is inconsistent with header declaration"<<endl;
                    }

                title="data are corrupted "+fileNameIn;
                dataX.freeMem();
                status=Cpdh::ERR_FILEOPEN;
            }




        fin.close();

}
//===================================================================================
void Cpdh::openPdhLhs()
{
ifstream fin(fileNameIn);


            if(!fin){
                cerr<<"ERROR: file couldn't be opened"<<endl;
            return;
            }

const size_t lineSize=64;
char line[lineSize];
string cmd,fver;
size_t rows=0,cols=0,atomTypes=0,atomsNumber;
vector<string> atomNames;
vector<size_t> &atomNamesNumber=grain->atomNamesNumber;
size_t testHeader=0;

            atomNamesNumber.clear();

            while(fin.peek()=='#' && !fin.eof() ){

                fin>>cmd;
                if(DB) {cout<<cmd<<endl;}


                //if(cmd.find("#ver:")!=string::npos) { fin>>fver;}
                if(cmd.find("#sizeRC:")!=string::npos)    {testHeader+=1;  fin>>rows>>cols;}
                if(cmd.find("#numOfatoms:")!=string::npos){testHeader+=10; fin>>atomsNumber;}
                if(cmd.find("#binWidth:")!=string::npos)  {testHeader+=100 ;fin>>bin;}


                if(cmd.find("#atoms:")!=string::npos){
                    testHeader+=1000;


                    fin.getline(line,lineSize,'\n');

                vector<string> tok(split<string>(std::string(line)," \t"));


                    atomTypes=tok.size();
                    atomNames.reserve(atomTypes);

                    for(size_t i=0;i<atomTypes;i++){
                        atomNames.push_back(tok[i]);
                    }


                continue;
                }

                //pobiera all do konca linii
                while(fin.get()!='\n');


             }

                if(atomTypes==1) {
                    atomNamesNumber.push_back(atomsNumber);
                }
                else{
                    //  ..... to be done soon as possible ........

                    #pragma message( "WARNING: cpdh.cpp: reading multitype systems is not implemented" )
                }


const size_t pdhcols=(atomTypes*(atomTypes+1)/2);


             if(! fin.good () ||
                !(testHeader==1111) ||
                atomNames.size()!=atomTypes ||  pdhcols!=(cols-1) || !rows){
                cerr<<"ERROR: the header of input file is corrupted: "<< fileNameIn<<endl;
                cerr<<"       testHeader="<<testHeader<<endl;
                fin.close();
             return;
             }


size_t i,j;

            try{

                dataX.allocMem(rows);
                //fin.exceptions(ios::badbit | ios::failbit );

                    if(cols==2){

                            dataYii.allocMem(rows);

                            for(i=0;i<rows;i++){
                                fin>>dataX[i]>>dataYii[i];

                                if(!fin.good () )
                                throw 1;

                                //consistency data checking; the end of line should be empty
                                fin.getline(line,lineSize);
                                if(fin.gcount()>1){
                                string str(line);
                                    if(!ltrim(str).empty())
                                        throw 2;
                                }


                            }
                    }
                    else{
                            dataYnn.resize(pdhcols);

                            for(i=0;i<pdhcols;i++)
                                dataYnn[i].allocMem(rows);

                            for(i=0;i<rows;i++){
                                fin>>dataX[i];

                                for(j=0;j<pdhcols;j++)
                                    fin>>dataYnn[j][i];


                                if(!fin.good () ){
                                    if( fin.eof() && i==rows-1)
                                        fin.clear();
                                    else
                                    throw 1;
                                 }

                                //consistency data checking; the end of line should be empty
                                fin.getline(line,lineSize);
                                if(fin.gcount()>1){
                                string str(line);
                                    if(!ltrim(str).empty())
                                        throw 2;
                                }
                            }
                    }

                    title="data taken from file: "+fileNameIn;
                    status=Cpdh::OK;

                    grain->atomTypes.clear();
                    grain->atomTypes.assign(atomNames.begin(),atomNames.end());


                    //try to guess the number of atoms for each type
                    if(atomNamesNumber.empty()){
                    size_t nn,N,sumTot;
                    double sum;

                            grain->atomNamesNumber.clear();
                            grain->atomNamesNumber.resize(cols);


                            for(i=0,nn=0,sumTot=0;i<atomTypes;i++){

                                sum=0;
                                for(j=0;j<rows;j++)
                                    sum+=dataYnn[nn][j];

                                N=0.5+sqrt(0.25+2*sum);
                                grain->atomNamesNumber[i]=N;
                                sumTot+=N;

                                nn+=atomTypes-i;
                            }

                            cout<<"cpdh.cpp: "<<__LINE__<<" N="<<N<<endl;


                            if(sumTot!=atomsNumber){
                                cerr<<"WARNING: discrepancy between calculated number of atoms for each type and number of atoms given in header ";
                                cerr<<(float)(sumTot)/(float)(atomsNumber)<<endl;

                            }
                    }



            }
            catch (std::exception e){

                    if(fin.eof()){
                        if(i<rows-1){
                            cerr<<"ERROR: wrong number of rows "<<endl;
                            dataX.freeMem();
                            status=Cpdh::ERR_FILEOPEN;
                         }
                         else //it's ok
                            title="data taken from file: "+fileNameIn;
                    }
                    else{
                        title="data are corrupted "+fileNameIn;
                        dataX.freeMem();
                        status=Cpdh::ERR_FILEOPEN;
                    }


            }
            catch (int e){

                    if(e==1 && i<rows-1){
                        cerr<<"ERROR: wrong number of rows "<<endl;
                            dataX.freeMem();
                            status=Cpdh::ERR_FILEOPEN;
                    }
                    else{
                        if(e==2)
                            cerr<<"ERROR: number of columns is inconsistent with header declaration"<<endl;
                    }

                title="data are corrupted "+fileNameIn;
                dataX.freeMem();
                status=Cpdh::ERR_FILEOPEN;
            }




        fin.close();


}
//===================================================================================
void Cpdh::openBinFilePLS()
{
fstream fin(fileNameIn,ios::in | ios::binary);


         if(!fin){
            cerr<<"ERROR: couldn't open file for saving, PDH results will be lost"<<endl;
            status=Cpdh::ERR_FILEOPEN;
        return;
        }

//fout.write( (char *) version.c_str(), version.length()*sizeof(char));

char ver[]={0,0,0,0};

            fin.read( (char *) ver,3);

            cout<<ver<<endl;

size_t atypSize;
            fin.read( (char *) &atypSize, sizeof(size_t));
            cout<<"sizeof(size_t)="<<sizeof(size_t)<<", "<<atypSize<<endl;

size_t aa;
            fin.read( (char *) &aa, sizeof(size_t));
            cout<<aa<<endl;

char atype[]={0,0,0,0};



            fin.read( atype, 1);

            cout<<atype<<":"<<fin.tellg()<<endl;

size_t numOfBins;

            fin.read( (char *)&numOfBins,sizeof(size_t));

            cout<<numOfBins<<endl;

double minBin,maxBin;

             fin.read( (char *) &minBin, sizeof(double));
             fin.read( (char *) &maxBin,sizeof(double));

             cout<<"sizeof(double)="<<sizeof(double)<<"    "<<minBin<<",   "<<maxBin<<endl;


}
//===================================================================================
void Cpdh::dispParams()
{
        cout<<"*** pdh parameters ***"<<endl;

        if(fileNameIn.empty())
            cout<<(*this);
        else{
            cout<<"input file : "<<fileNameIn<<endl;
            cout<<"X size "<<dataX.size<<endl;
            cout<<"Y size "<< ((dataYnn.empty()) ? 1 : dataYnn.size())<<endl;
            cout<<"bin width "<<bin<<endl;
        }


}

//===================================================================================

#ifdef __linux__
void Cpdh::plotFigure()
{
FILE *plotHandle = NULL;

            if(plotHandle == NULL){
              plotHandle = popen("gnuplot -persist 2>/dev/null", "w");
              assert(plotHandle != NULL);
            }

            fprintf(plotHandle,"set xlabel \"r\"\n");
            fprintf(plotHandle,"set ylabel \"frequency\"\n");
            fprintf(plotHandle,"set xticks font \",16\"\n");
            fprintf(plotHandle,"set yticks font \",16\"\n");
            fprintf(plotHandle,"set grid\n");

            plotData(plotHandle);
            pclose(plotHandle);
}
//===================================================================================
void Cpdh::saveFigure(const string &fileName)
{
FILE *plotHandle = NULL;

            if(plotHandle == NULL){
              plotHandle = popen("gnuplot -persist 2>/dev/null", "w");
              assert(plotHandle != NULL);
            }

const std::string cmdSetTerm("set terminal png size "+figWidth+", "+figHeight+"\n");
const std::string cmdSetOutput("set output '"+fileName+"'\n");


            fprintf(plotHandle,cmdSetTerm.c_str());
            fprintf(plotHandle,cmdSetOutput.c_str());
            fprintf(plotHandle,"set xlabel \"r\"\n");
            fprintf(plotHandle,"set ylabel \"frequency\"\n");
            fprintf(plotHandle,"set xticks font \",16\"\n");
            fprintf(plotHandle,"set yticks font \",16\"\n");
            fprintf(plotHandle,"set grid\n");



            plotData(plotHandle);
            pclose(plotHandle);
}
//===================================================================================
void Cpdh::plotData(FILE *plotHandle)
{
vector<string> tokens(split<string> (plotPrm,","));
const std::string sizeType{(sizeof(position)==sizeof(double))? "double" : "float"};

const double *x,*y;
const size_t dataSize=dataX.size;
position *plotData=new position[2*dataSize];
size_t i,j;
std::string plotcmd("plot ");



             if(grain->atomTypes.size()==1){
             vector<string> subtokens(split<string>(tokens[0],":"));
             vector<string> sub0(split<string>(subtokens[0]," \t"));
             vector<string> sub1(split<string>(subtokens[1]," \t"));

                 plotcmd+="'-' binary format=\"%%2"+sizeType+"\" record="+std::to_string(dataSize)+" endian=little u 1:2 w imp title 'PDH for "+grain->atomTypes[0].name +"-"+grain->atomTypes[0].name+"'\n";
                 fprintf(plotHandle,plotcmd.c_str());

                 try{

                 const int aX(std::stoi(sub0.back()));
                 const int aY(std::stoi(sub1[0]));

                     if(aX!=1){
                         errMsg(" wrong number of first column: "+sub0.back()+". The first column must be 1");
                     throw 1;
                     }
                     if(aY!=2 ){
                         errMsg(" wrong number of second column: "+sub1[0]+". The second column must be 2 ");
                     throw 2;
                     }

                     x=dataX.values;
                     y=dataYii.values;

                     for(i=0,j=0;i<dataSize;i++,x++,y++){
                         plotData[j++]=(*x);
                         plotData[j++]=(*y);
                     }

                     fwrite(plotData,sizeof(position),2*dataSize,plotHandle);
                 }
                 catch(...){

                 }


             }
             else{
             csize numOfAtomTypes=grain->atomTypes.size();
             csize mpdhSize=(numOfAtomTypes+1)*numOfAtomTypes/2;
             vector <string> atypes;
                            atypes.reserve(mpdhSize);

                for (i=0;i<numOfAtomTypes;i++){
                     for(j=i;j<numOfAtomTypes;j++){
                     string atomAB("PDH for "+(grain->atomTypes[i].name)+"-"+grain->atomTypes[j].name);
                            atypes.emplace_back(atomAB);
                     }
                }



                    for(auto &token : tokens){
                        vector<string> subtokens(split<string>(token,":"));
                        //vector<string> sub0(split<string>(subtokens[0]," \t"));
                        vector<string> sub1(split<string>(subtokens[1]," \t"));

                        plotcmd+="'-' binary format=\"%%2"+sizeType+"\" record="+std::to_string(dataSize)+" endian=little u 1:2 w "+sub1[2]+" title '"+atypes[std::stoi(sub1[0])-2]+"',";
                    }

                    plotcmd.back()='\n';

                    fprintf(plotHandle,plotcmd.c_str());

                    for(auto &token : tokens){
                    vector<string> subtokens(split<string>(token,":"));
                    vector<string> sub0(split<string>(subtokens[0]," \t"));
                    vector<string> sub1(split<string>(subtokens[1]," \t"));

                    const int aX(std::stoi(sub0.back()));
                    const int aY(std::stoi(sub1[0]));

                        if(aX!=1){
                            errMsg(" wrong number of first column: "+sub0.back()+". The first column must be 1");
                        continue;
                        }

                        if(aY<2 || aY>=(dataYnn.size()+2) ){
                            errMsg(" wrong number of second column: "+sub1[0]+". The second column must be greater than 1 and less than "+std::to_string(dataYnn.size()-2));
                        continue;
                        }

                        x=dataX.values;
                        y=dataYnn[aY-2].values;

                        for(i=0,j=0;i<dataSize;i++,x++,y++){
                            plotData[j++]=(*x);
                            plotData[j++]=(*y);
                        }

                        fwrite(plotData,sizeof(position),2*dataSize,plotHandle);
                    }

             }

            delete [] plotData;

}


#endif
//===================================================================================
ostream & operator<<(ostream &o, Cpdh &cpdh)
{
        o<<"#bin: "<<cpdh.bin<<endl;
        o<<"#mode: "<<( (cpdh.binMode)?"double":"single")<<endl;
        o<<"#threads: "<<cpdh.threads<<endl;

        if(!cpdh.comment.empty())
            o<<"#comment: "<<cpdh.comment<<endl;

return o;
}
//===================================================================================
bool Cpdh::parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__)
{
 const str send("end");

            clearData();
            ptr_uvar=uvar__;

            while(cmd[index]!=send){


                if(cmd[index]=="bin"){
                    bin=cmd[index++][1];
                continue;
                }

                if(cmd[index]=="mode"){
                    binMode=(cmd[index++][1]=="single")?SINGLE:DOUBLE;
                continue;
                }
				
				if(cmd[index]=="mthmode"){
					mthmode=(cmd[index++][1]=="openmp")?OPENMP:STDTHREAD;
				continue;
				}

                if(cmd[index]=="open"){
                    fileNameIn=cmd[index++][1];
                    Script::replaceVars(ptr_uvar,fileNameIn);
                continue;
                }

                if(cmd[index]=="rangeN"){
                const double from=std::stod(cmd[index][1]);
                const int    N=   std::stoi(cmd[index][2]);
                const double to=  std::stod(cmd[index][3]);
                const double step=(to-from)/(N-1);

                        range=cmd[index][1]+" "+std::to_string(step)+" "+cmd[index][3];

                    index++;

                continue;
                }


                if(cmd[index]=="threads"){
                    threads=cmd[index++][1];
                continue;
                }

                if(cmd[index]=="save"){
                    fileName=cmd[index++][1];
                    Script::replaceVars(ptr_uvar,fileName);
                continue;
                }

                if(cmd[index]=="saveFigure"){
                string fileName(cmd[index][1]);
                        Script::replaceVars(ptr_uvar,fileName);
                        figFileName=fileName;

                        if(cmd[index].numOfKeyValues()==4){
                            figWidth=cmd[index][2];
                            figHeight=cmd[index][3];
                        }
                        else{
                            figWidth="800";
                            figHeight="600";
                        }

                        index++;
                continue;
                }

                if( cmd[index].isKey("difftime")){
                    diffTime=true;
                    index++;
                continue;
                }

                if ( cmd[index].isKey("plot")){
                    plotPrm=cmd[index][1];

                    // if gnuplot statement is empty add a space
                    if(plotPrm.empty()) plotPrm+=" ";

                    index++;
                continue;
                }

                if( cmd[index].isKey("printprm") || cmd[index].isKey("printPrm")){
                    printPrm=true;
                    index++;
                continue;
                }

                if(cmd[index].isKey("comment")){
                    comment=cmd[index++][1];
                continue;
                }


                cerr<<"Error: cpdh.cpp  unknown command "<< cmd[index]<<" line "<<__LINE__<<endl;
            return false;
            }

			//cout<<"mthmode : "<<mthmode<<endl;
            //ptrCalc=(grain->atomTypes.size()==1) ? &Cpdh::monoLatticePdh : &Cpdh::multiLatticePdh;
			if(mthmode==Cpdh::OPENMP)
				ptrCalc=(grain->atomTypes.size()==1) ? &Cpdh::monoLatticePdh : &Cpdh::multiLatticePdh;
			else
				ptrCalc=&Cpdh::pdhByCppThreads;
return true;
}
//===================================================================================
void Cpdh::calc()
{

        if(fileNameIn.empty()){

                cout<<" start PDH computation"<<endl;
                if(printPrm)
                    dispParams();


                (this->*ptrCalc)();


                if(diffTime)
                    cout<<" computation time (seconds): "<<std::difftime(stopTime,startTime)<<endl;
        }
        else{
            status=Cpdh::RUN;
            openFile();

            if(printPrm)
                    dispParams();
        }

        if(!fileName.empty() && status==Cpdh::OK)
            saveResults();


#ifdef __linux__
        if(!plotPrm.empty())plotFigure();
        if(!figFileName.empty()) saveFigure(figFileName);
#endif

}
//===================================================================================
void Cpdh::import(dataPdh &inputX, dataPdh &inputY)
{
        dataX=std::move(inputX);
        dataYii=std::move(inputY);

}
//===================================================================================
struct StRegPdh
{
static NanoGrain::StNanoGrain *grain;
static position ibin;
static CProgress *progress;
static size_t pdhsize;
size_t iStart,iStop,ith;

vector<position> pdh;


	void operator ()()
	{
		
		//cout<<ith<<"  "<<iStart<<" : "<<iStop<<endl;
		
	NanoGrain::StAtom *patom0,*patom1;
	NanoGrain::StAtom *patoms=grain->atoms.data();
	size_t bin,pi,pj;	
	const size_t numOfatoms=grain->atoms.size();
		
		pdh.resize(pdhsize);
		for(auto &v : pdh)
			v=0;
		
		
		for (pi=iStart;pi<iStop;pi++){
			patom0=patoms+pi;
			
			/// to do : mutex
			//progress->next();
			(*progress)++;
			
			for(pj=pi+1,patom1=patom0+1;pj<numOfatoms;pj++,patom1++){

				bin=(size_t)     sqrt(  sqrd(patom1->x-patom0->x)+
										sqrd(patom1->y-patom0->y)+
										sqrd(patom1->z-patom0->z));
				

				(*(pdh.data()+bin))++;
				
			}	
		}

	}
		
};

NanoGrain::StNanoGrain *StRegPdh::grain=nullptr;
position StRegPdh::ibin=1;
CProgress *StRegPdh::progress=nullptr;
size_t StRegPdh::pdhsize=0;



void Cpdh::pdhByCppThreads()
{
const position wbin=getBinWidth();
const position ibin=1.0/wbin;
const size_t numOfatoms=grain->atoms.size();
const size_t numOfthreads=stoi(threads);
const size_t pdhsize=std::ceil(grain->maxR*2.25*ibin);
CProgress progress;


			StRegPdh::ibin=ibin;
			StRegPdh::grain=grain;
			StRegPdh::pdhsize=pdhsize;
			StRegPdh::progress=&progress;
			
			
StRegPdh  vregpdh[numOfthreads];			
const size_t frac=0.75*numOfatoms/numOfthreads;	
size_t i;
			vregpdh[0].iStart=0;
			vregpdh[0].ith=0;

			for(i=1;i<numOfthreads;i++){
				
				vregpdh[i-1].iStop=i*frac;
								
				vregpdh[ i ].iStart = vregpdh[i-1].iStop ;
				vregpdh[ i ].ith=i;
				
			}
			
			vregpdh[ i-1 ].iStop=numOfatoms-1;
			
			/// optimalization on
			for(size_t i=0;i<numOfatoms;i++){
                grain->atoms[i]*=ibin;
            }
			
			
std::thread vthread[numOfthreads];

			progress.start(numOfatoms-1);
			
			startTime = std::time(NULL);
            status=Cpdh::RUN;

			for(i=0;i<numOfthreads;i++){
				vthread[i]=std::thread(std::ref(vregpdh[i]));
				//vthread[i].join();
				
			}
					
			
			dataX.allocMem(pdhsize);
            dataYii.allocMem(pdhsize);


            for(size_t i=0;i<pdhsize;i++){
                dataX[i]=i*wbin;
                dataYii[i]=0;
            }			
			
			
			for(i=0;i<numOfthreads;i++){
				
				vthread[i].join();
								
				for(size_t j=0;j<pdhsize;j++){
					dataYii[j]+=vregpdh[i].pdh[j];
				}
				
			}
			
			
			stopTime = std::time(NULL);
			progress.stop();
		
		
		
			/// optimalization off
			for(size_t i=0;i<numOfatoms;i++){
                grain->atoms[i]*=wbin;
            }
			
			status=Cpdh::OK;
}






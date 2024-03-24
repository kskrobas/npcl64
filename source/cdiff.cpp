/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cdiff.cpp
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



#include "cdiff.h"
#include "cprogress.h"
#include "crandom.h"
#include "colormsg.h"
#include "createdir.h"

#include <omp.h>
#include <iomanip>
//#include <chrono>
#include <thread>
#include <future>
#include <cstdlib>
#include <regex>
#include <assert.h>

//#include <random>

inline position sqr(const position &x){return x*x;}
inline position cub(const position &x){return x*x*x;}
inline double sqrd(const double &x){return x*x;}



#ifdef __linux__

#else
    # define M_PI		3.14159265358979323846	/* pi */
    # define M_1_PI	0.318309886183790671537767526745028724L /* 1/pi */
#endif
//#elif _WIN32 or _WIN64

#if DEBUG
#define DB 1
#else
#define	DB 0
#endif



position sincfast(const position &x)
{return (x>318*M_PI) ? 0.0 : std::sin(x)/x;}

position sinc(const position &x)
{return sin(x)/x;}



ostream& endlu(ostream&  o)
{
    return o<<"\n";
}

ostream& endld(ostream&  o)
{
    return o<<"\r\n";
}


//=============================================================================


Cdiff::Cdiff()
{
ifstream fileSc;

            //test current directory
            fileSc.open("scFact.sft");

            if(!fileSc.is_open()){  //test home directory
            const char *homeDir=std::getenv("NPCLPATH");

                    if(homeDir){
                    const string shomeDir(homeDir);
                    const string fileName(shomeDir+"/scFact.sft");

                        fileSc.open(fileName);
                    }
            }


            if(!fileSc)
                cerr<<"WARNING: file scFact.sft is missing"<<endl;


            fileSc.close();
            scfactors=nullptr;
}
//===================================================================================
bool Cdiff::parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__)
{
const str send("end");

            ptr_uvar=uvar__;
            presetData();

            while(cmd[index]!=send){

                if(DB) cout<<cmd[index]<<endl;


                if( cmd[index].isKey("difftime")){
                    diffTime=true;
                    index++;
                continue;
                }

                if(cmd[index]=="dw"){
                    dw=(cmd[index++][1]);
                continue;
                }


                if(cmd[index]=="extrapolate"){
                    extrapolate=(cmd[index++][1]=="yes");
                continue;
                }

                if(cmd[index]=="fastsinc"){
                    fastsinc=(cmd[index++][1]=="yes");
                continue;
                }

                if(cmd[index]=="Ki"){
                    ki=cmd[index++][1];
                continue;
                }

                if(cmd[index]=="lambda"){
                    lambda=cmd[index][1];
                    index++;
                continue;
                }

                if(cmd[index]=="mode"){
                    mode=(cmd[index++][1]=="laue")? Cdiff::laue  : Cdiff::debyea;
                continue;
                }

                if(cmd[index]=="norm"){

                    if(cmd[index][1].empty())
                        norm=true;
                    else
                    norm=(cmd[index][1]=="yes");

                    index++;
                continue;
                }

                if(cmd[index]=="noiseSQ"){
                    noise=cmd[index][1]+" "+cmd[index][2]+" "+cmd[index][3];
                    index++;
                continue;
                }

                if(cmd[index]=="open"){
                string fileName(cmd[index++][1]);
                    Script::replaceVars(ptr_uvar,fileName);
                    fileNameIn=fileName;
                continue;
                }

                if ( cmd[index].isKey("plot")){
                    plotPrm=cmd[index][1];

                    // if gnuplot statement is empty add a space
                    if(plotPrm.empty()) plotPrm+=" ";

                    index++;
                continue;
                }

                if(cmd[index]=="polarization"){
                    polarization=(cmd[index++][1]=="yes");
                continue;
                }


                if(cmd[index]=="radiation"){

                    if(cmd[index][1]=="xray")
                        radiation=XRAY;
                    else
                        if(cmd[index][1]=="neutron")
                            radiation=NEUTR;
                        else
                            if(cmd[index][1]=="electron")
                                radiation=ELECT;
                            else
                                radiation=OFF;

                    index++;
                continue;
                }


                if(cmd[index]=="range"){
                    range=cmd[index][1]+" "+cmd[index][2]+" "+cmd[index][3];
                    index++;
                continue;
                }


                if(cmd[index]=="save"){
                string fileName(cmd[index++][1]);
                    Script::replaceVars(ptr_uvar,fileName);
                    fileNameOut.push_back(fileName);
                continue;
                }


                if(cmd[index]=="saveopt"){

                    for(int i=1; i<cmd[index].numOfKeyValues(); i+=2){
                    const string key  (cmd[index][i]);
                    const string value(cmd[index][i+1]);

                        if(key=="axis"){
                        const int xyz=(int)(std::tolower(value[0]))-'x';
                            saveopt.axis=(StSaveOpt::EAXIS)(xyz);
                        continue;
                        }
                        if(key=="format"){
                            saveopt.format=(value=="short")? StSaveOpt::SHORT : StSaveOpt::LONG;
                        continue;
                        }
                        warnMsg("unknown save option "+key);
                    }

                    index++;
                continue;
                }


                if(cmd[index]=="filesft"){
                    fileSft=cmd[index++][1];
                     Script::replaceVars(ptr_uvar,fileSft);
                continue;
                }


                if(cmd[index]=="threads"){
                    threads=cmd[index++][1];
                continue;
                }


                if(cmd[index]=="theta"){
                    theta=(cmd[index++][1]=="yes");
                continue;
                }


                if(cmd[index]=="unix2dox"){
                    u2dendl=true;
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


                cerr<<"Error: cpdiff.cpp  unknown command "<< cmd[index]<<" line "<<__LINE__<<endl;
            return false;
            }


            if(fileNameIn.empty()){
                if(range.empty() ){errMsg("ERROR: diffraction range undefined"); return false;}
                if(lambda.empty()){errMsg("ERROR: diffraction lambda undefined"); return false;}
            }


return true;
}
//===================================================================================
void Cdiff::calc()
{

    try{

        if(fileNameIn.empty()){
            if(mode==Cdiff::debyea){

                if(pdh->dataX.empty){
                    cerr<<"ERROR: pdh is empty"<<endl;
                throw Ediffstatus::ERR_EMPTYPDH;
                }

                cout<<" start DIFF computation"<<endl;
                if(printPrm)
                    dispParams();


                if(pdh->grain->atomTypes.size()==1)
                    debDiff_mono();
                else
                    debDiff_multi();

            
            const double comptime=std::difftime(stopTime,startTime);


                if(diffTime) cout<<" computation time (seconds): "<<comptime<<endl;

                if(polarization) polarizationI();

                if(!dw.empty()) debyeaWaller();

                if(!noise.empty()) noiseSQ();

                if(norm)  normIS();

                saveResults();

                if(!plotPrm.empty()) plot();

            }
            else { //laue
                laueDiff();
                saveResults();

            const double comptime=std::difftime(stopTime,startTime);
                if(diffTime) cout<<" computation time (seconds): "<<comptime<<endl;
            }
        }
        else{
            openFile();
            if(!plotPrm.empty()) plot();
        }

    }
    catch(Ediffstatus e){
        cerr<<"ERROR: ";
        switch (e){
        case ERR_EMPTYPDH: cerr<<"empty PDH"; break;
        case ERR_SCATCOEFF: cerr<<"missing scattering factors (scFAct.sft)";break;
        default: cerr<<" diffraction block, code "<<e;break;
        }

        delete [] scfactors;
        scfactors=nullptr;

        cerr<<endl;
        throw e;
    }
    catch (const std::out_of_range& e) {
       std::cerr << "ERROR: Out of Range error: " << e.what() << endl;
     }
    catch(Script::ResultRepVar e){
        throw e;
    }
    catch(...){
        cerr<<"diffraction block"<<endl;
        throw 0;
    }
}
//-----------------------------------------------------------------------------

//===================================================================================
void Cdiff::debDiff_mono()
{           
            if(DB) cout<<" monolattice Diffraction computations "<<endl;            

            if(radiation==OFF)
                ptr_scattfactor =&Cdiff::noScatFact;
            else{
                /// presetData() procedure  deletes scfactors
                scfactors=new Uscattfactors[1];
                if(!loadScattFactors(pdh->grain->atomTypes[0].name,0))
                throw Ediffstatus::ERR_SCATCOEFF;

                  switch(radiation){
                  case XRAY: ptr_scattfactor=&Cdiff::xrayScatEqu; break;
                  case NEUTR: ptr_scattfactor=&Cdiff::neutronCrossSec; break;
                  case ELECT: ptr_scattfactor=&Cdiff::electronScatEqu; break;
                  default: &Cdiff::noScatFact;
                  }
            }
                            

vector<string> ststst(split<string>(range," "));//start step stop
const double start=std::stod(ststst[0]);
const double step =std::stod(ststst[1]);
const double stop =std::stod(ststst[2]);


const int extraSize=(extrapolate) ? start/step : 0;
const int thetaSize=1+std::floor( (stop-start)/step);

            if(thetaSize<2 || extraSize<0){
                cerr<<"ERROR: wrong values of 'range' directive"<<endl;
            throw Ediffstatus::ERR_RANGE;
            }

const size_t totSize=thetaSize+extraSize;
position *binX= pdh->dataX.values;
position *binY= pdh->dataYii.values;
const size_t binSize=pdh->dataSize();
const size_t N=pdh->grain->atomNamesNumber[0];


            allocMem(totSize);

position * const dataX=t.values;
position * const dataQ=q.values;
position * const dataI=I.values;
position * const dataS=S.values;



            //// preliminary steps of diffraction calculation algorithm

                //// 1. searching non zero bins

size_t noneZeroBinsSize=0;
size_t i,j;

            omp_set_num_threads(std::stoi(threads));

            #pragma omp parallel for reduction(+:noneZeroBinsSize) if(binSize>10000)
            for(i=0;i<binSize;i++)
                if(binY[i]>0)
                    noneZeroBinsSize++;

size_t *noneZeroBins=new size_t[noneZeroBinsSize];

                /// save indices of non zero bins
            for(i=0,j=0;i<binSize;i++)
                if(binY[i]>0)
                    noneZeroBins[j++]=i;


                ////// 2. setting initial values
vector<string> lambdaToks{split<string>(lambda," ")};
std::map<std::string, std::string>::const_iterator cathLambda=Elements::xrad.find(lambdaToks[0]);
const position lambda__=(cathLambda!=Elements::xrad.end()) ? std::stod(cathLambda->second) : std::stod(lambdaToks[0]);

const position k1=M_PI/360;
const position invWL=1.0/lambda__;
const position k2=4*M_PI*invWL;
const position startAlpha=(start-extraSize*step)*k1;
const position stepAlpha=step*k1;
position alpha,th,lsin,qx;

                ////// 3. extrapolate

           for(int i=0;i<extraSize;i++){
                alpha=step*i;
                th=alpha*k1;
                lsin=sin(th);
                qx=k2*lsin;
                dataX[i]=2*th;  // radians
                dataQ[i]=qx;     // Q
                dataI[i]=0;
                dataS[i]=0;
            }


            ///// start MAIN TASK /////
position qrij,fi,suma,fiArg;
const position k0=2.0/N;

position (*ptr_sinc)(const position &)= (fastsinc) ? &sincfast : &sinc;
CProgress progress;

			progress.title=std::string(" diffraction ");
			progress.start(thetaSize);


            startTime = std::time(nullptr);

            #pragma omp parallel for private(alpha,th,lsin,qx,fiArg,fi,suma,j,qrij)
            for(i=extraSize;i<totSize; i++){
				
                #pragma omp critical
                {
					progress++;
                }

                th=startAlpha+stepAlpha*i;
                lsin=sin(th);
                qx=k2*lsin;
                fiArg=lsin*invWL;                
                fi=(this->*ptr_scattfactor)(fiArg,0);

                dataX[i]=2*th;
                dataQ[i]=qx;

                suma=0;

                //#pragma omp parallel for reduction(+:suma) private(qrij)
                for(j=0;j<noneZeroBinsSize;j++){
                    qrij=qx*binX[noneZeroBins[j]];
                    suma+=binY[noneZeroBins[j]]*ptr_sinc(qrij);
                }

                dataS[i]=k0*suma+1;//phys. rev. b 73, 184113 2006
                dataI[i]=N*fi*fi*dataS[i]; 
            }

            stopTime = std::time(nullptr);
            delete [] noneZeroBins;

            delete [] scfactors;
            scfactors=nullptr;

}
//===================================================================================
void Cdiff::debDiff_multi()
{

                if(pdh->dataYnn.empty()){
                    cerr<<"ERROR: empty PDH"<<endl;
                    diffstatus=Cdiff::ERR_EMPTYPDH;
                return;
                }

const size_t numOfAtomTypes=pdh->grain->atomTypes.size();
vector<string> ststst(split<string>(range," "));//start step stop
const double start=std::stod(ststst[0]);
const double step =std::stod(ststst[1]);
const double stop =std::stod(ststst[2]);
const int extraSize=(extrapolate) ? start/step : 0;
const int thetaSize=1+std::floor( (stop-start)/step);
const size_t totSize=thetaSize+extraSize;


                if(thetaSize<2 || extraSize<0){
                    cerr<<"ERROR: wrong values of 'range' directive"<<endl;
                throw Ediffstatus::ERR_RANGE;
                }

                allocMem(totSize);



                if(radiation==OFF)
                    ptr_scattfactor =&Cdiff::noScatFact;
                else{
                    /// presetData() procedure  deletes scfactors
                    scfactors=new Uscattfactors[numOfAtomTypes];

                    for(size_t i=0;i<numOfAtomTypes;i++){
                            if(!loadScattFactors(pdh->grain->atomTypes[i].name,i))
                            throw Ediffstatus::ERR_SCATCOEFF;
                   }

                    switch(radiation){
                    case XRAY: ptr_scattfactor=&Cdiff::xrayScatEqu; break;
                    case NEUTR: ptr_scattfactor=&Cdiff::neutronCrossSec; break;
                    case ELECT: ptr_scattfactor=&Cdiff::electronScatEqu; break;
                    default: &Cdiff::noScatFact;
                    }
                }

position * const dataX=t.values;
position * const dataQ=q.values;
position * const dataI=I.values;
position * const dataS=S.values;

std::map<std::string, std::string>::const_iterator cathLambda=Elements::xrad.find(lambda);
const position lambda__=(cathLambda!=Elements::xrad.end()) ? std::stod(cathLambda->second) : std::stod(lambda);

const position k1=M_PI/360;
const position invWL=1.0/lambda__;
const position k2=4*M_PI*invWL;
const position startAlpha=(start-extraSize*step)*k1;
const position stepAlpha=step*k1;

position fi,fj;
position alpha,th,lsin,qx,fiArg;
size_t Nf2,nn;
position sumBins,sumCoh,sumTot,sumNf2;
position (*ptr_sinc)(const position &)=(fastsinc) ? &sincfast : &sinc;


                //// prepare nonzero PDH distributions

                buildNonZeroBinsPdh();

                ////// . extrapolate
                   for(int i=0;i<extraSize;i++){
                        alpha=step*i;
                        th=alpha*k1;
                        lsin=sin(th);
                        qx=k2*lsin;
                        dataX[i]=2*th;  // radians
                        dataQ[i]=qx;     // Q
                        dataI[i]=0;
                        dataS[i]=0;
                    }


                   omp_set_num_threads(std::stoi(threads));
CProgress progress;

					progress.title=(std::string(" diffraction "));
					progress.start(totSize);

                   startTime = std::time(NULL);

                   /// diffraction calculations
                   for(size_t k=extraSize;k<totSize;k++,progress++){

                       th=startAlpha+stepAlpha*k;
                       lsin=sin(th);
                       qx=k2*lsin;
                       //fiArg=sqr(lsin*invWL);
                       fiArg=lsin*invWL;
                       dataX[k]=2*th;
                       dataQ[k]=qx;

                       sumTot=0;
                       sumNf2=0;
                       nn=0;/// select pdh distribution


                       for(size_t i=0;i<numOfAtomTypes;i++){                           
                           fi=(this->*ptr_scattfactor)(fiArg,i);
                            sumCoh=0;

                            for(size_t j=i;j<numOfAtomTypes;j++,nn++){                                
                                fj=(this->*ptr_scattfactor)(fiArg,j);
                                sumBins=0;

                                    #pragma omp parallel
                                    {

                                    StNonZeroBins  *ptr_pdh=(pdhNZB.data()+nn);
                                    const size_t pdhSize=ptr_pdh->dataX.size();
                                    position *ptr_valuesX=ptr_pdh->dataX.data();
                                    position *ptr_valuesY=ptr_pdh->dataY.data();

                                            #pragma omp for reduction(+:sumBins)
                                            for(size_t h=0;h<pdhSize;h++)
                                                sumBins+=ptr_valuesY[h]*ptr_sinc(ptr_valuesX[h]*qx);
                                    }

                                 sumCoh+=fi*fj*sumBins;
                            }

                            Nf2=pdh->grain->atomNamesNumber[i]*fi*fi;

                            sumTot+=( Nf2+ 2*sumCoh) ;
                            sumNf2+=Nf2;
                        }

                        dataI[k]=sumTot;
                        dataS[k]=sumTot/sumNf2;// J. Mol. Str 383 1996 303-308

                       }


                   stopTime = std::time(NULL);

}

//====================================================================================
struct StKivector{
position x,y,z;
enum EKistatus{OK,ERR_SZERO};

    StKivector(const string &Ki,const double lambda)
    {
    vector<string> kisplit(split<string>(Ki," "));

        x=std::stod(kisplit[0]);
        y=std::stod(kisplit[1]);
        z=std::stod(kisplit[2]);

    const double kmod=std::sqrt(sqr(x)+sqr(y)+sqr(z));

        if(kmod<1e-9){
            cerr<<"ERROR: mod(k)<1e-9"<<endl;
        throw ERR_SZERO;
        }

    const double ikmod=1.0/lambda/kmod;
        x*=ikmod;
        y*=ikmod;
        z*=ikmod;
    }



};

//====================================================================================
void Cdiff::laueDiff()
{
const size_t numOfAtomTypes=pdh->grain->atomTypes.size();
vector<string> ststst(split<string>(range," "));//start step stop
const double start=std::stod(ststst[0]);
const double step =std::stod(ststst[1]);
const double stop =std::stod(ststst[2]);

const int KSize=1+std::floor( (stop-start)/step);
//position K;



std::map<std::string, std::string>::const_iterator cathLambda=Elements::xrad.find(lambda);
const position lambda__=(cathLambda!=Elements::xrad.end()) ? std::stod(cathLambda->second) : std::stod(lambda);

                if(KSize<2 ){
                    cerr<<"ERROR: wrong values of 'range' directive"<<endl;
                throw Ediffstatus::ERR_RANGE;
                }

                if(cathLambda==Elements::xrad.end() && lambda=="0"){
                    cerr<<"ERROR: wavelenght shouldn't be zero"<<endl;
                throw Ediffstatus::ERR_RANGE;
                }


                /////////////////////////////////////////////////////////////////////////////////////
                ///  definitions according to:                                                    ///
                ///  Williams D B and Carter C B 2009 .                                           ///
                ///  Diffraction in TEM Transmission Electron Microscopy Part 2: Diffraction,     ///
                /// (New York: Springer) pp 197–20.                                               ///
                ///                                                                               ///
                /// see also:                                                                     ///
                ///      https://docs.lammps.org/compute_xrd.html                                 ///
                ///                                                                               ///
                ///                                                                               ///
                ///   incident   wavevector:    ki=[0,0,1/lambda]                                 ///
                ///                                                                               ///
                /////////////////////////////////////////////////////////////////////////////////////

//const position ilambda=1.0/lambda__;
const StKivector Ki(ki,lambda__);
//const position deg2rad=M_PI/180.0;
const size_t N=pdh->grain->atomNamesNumber[0];
//position fi; //atomic scattering factor


            if(radiation==OFF)
                ptr_scattfactor =&Cdiff::noScatFact;
            else{

                /// presetData() procedure  deletes scfactors
                scfactors=new Uscattfactors[numOfAtomTypes];

                for(size_t i=0;i<numOfAtomTypes;i++){
                        if(!loadScattFactors(pdh->grain->atomTypes[i].name,i))
                        throw Ediffstatus::ERR_SCATCOEFF;
               }

                switch(radiation){
                case XRAY: ptr_scattfactor=&Cdiff::xrayScatEqu; break;
                case NEUTR: ptr_scattfactor=&Cdiff::neutronCrossSec; break;
                case ELECT: ptr_scattfactor=&Cdiff::electronScatEqu; break;
                case OFF: ptr_scattfactor=&Cdiff::noScatFact;
                }
            }



            kpoints.clear();
            kpoints.resize(KSize);

            for (size_t i=0;i<KSize;i++){
                kpoints[i]=(start+i*step);
            }

            laueDiffData.clear();
            laueDiffData.resize(KSize);


            for(auto &s: laueDiffData){
                s.resize(KSize);
                for(auto &ss: s){
                    ss.resize(KSize);
                    for(auto &sss: ss)
                        sss=0;
                }
            }

CProgress progress;

            progress.title=(std::string(" laue diffraction "));
            progress.start(N);

            startTime = std::time(NULL);

            omp_set_num_threads(std::stoi(threads));

            #pragma omp parallel
            {
            cx3dspace local_kSpace;
            position xkx,yky,zkz,arg,fi;
            position kx2,kxy2,kxyz2,sinThLambda;
            int i,j,k;


                    local_kSpace.resize(KSize);
                    for(auto &s: local_kSpace){
                        s.resize(KSize);
                        for(auto &ss: s){
                            ss.resize(KSize);
                            for(auto &sss: ss)
                                sss=0;
                        }
                    }


                    #pragma omp for
                    for(size_t a=0;a<N;a++){
                    auto &atom=pdh->grain->atoms[a];


                        for (i=0;i<KSize;i++){
                            xkx=atom.x*(kpoints[i]-Ki.x);
                            kx2=sqr(kpoints[i]);

                            for (j=0;j<KSize;j++){
                                yky=atom.y*(kpoints[j]-Ki.y);
                                kxy2=kx2+sqr(kpoints[j]);

                                for(k=0;k<KSize;k++){
                                    zkz=atom.z*(kpoints[k]-Ki.z);
                                    kxyz2=kxy2+sqr(kpoints[k]);

                                    sinThLambda=0.5*std::sqrt(kxyz2);
                                    fi=(this->*ptr_scattfactor)(sinThLambda,atom.atype);

                                    arg=2*M_PI*(xkx+yky+zkz);

                                    local_kSpace[i][j][k]+=fi*std::exp(cxposition(0,arg));

                                }
                            }
                        }

                        #pragma omp critical
                        progress++;
                    }



                    #pragma omp critical
                    {
                          for(i=0;i<KSize;i++)
                              for(j=0;j<KSize;j++)
                                  for(k=0;k<KSize;k++)
                                      laueDiffData[i][j][k]+=local_kSpace[i][j][k];                                                                                                       
                    }
            }


            stopTime = std::time(NULL);

}


//===================================================================================
//*** the function is searching for nonzero bins of PDH distributions; it helps to speed up
//    the diffraction computation by removing empty shots


void Cdiff::buildNonZeroBinsPdh()
{
const size_t origPdhNumber=pdh->dataYnn.size();
const size_t origDataSize =pdh->dataX.size;
size_t i,j;
position *ptr_valuesX,*ptr_valuesY;
vector<StNonZeroBins>::iterator iter_pdh;


                pdhNZB.resize(origPdhNumber);
                iter_pdh=pdhNZB.begin();

                for(i=0;i<origPdhNumber;i++,iter_pdh++){
                    iter_pdh->reserveMemory(origDataSize);
                    ptr_valuesX=pdh->dataX.values;
                    ptr_valuesY=(pdh->dataYnn.data()+i)->values;

                    for(j=0;j<origDataSize;j++,ptr_valuesX++,ptr_valuesY++)
                        if(*ptr_valuesY)
                            iter_pdh->push_values(ptr_valuesX,ptr_valuesY);

                    iter_pdh->shrink();
                }


}
//===================================================================================
void Cdiff::saveResults()
{

            for(string file : fileNameOut){

                    if(!createDirsIfDontExist(file)){
                        errMsg("couldn't create nested directories for "+file);
                    throw Cdiff::ERR_FILEOPEN;
                    }

                    if(file.rfind(".dat")!=string::npos){
                        if(mode==debyea) saveDatFile(file);
                        else warnMsg(" .dat format works only for Debye mode");
                    }
                    else
                        if(file.rfind(".diffb")!=string::npos){
                            if(mode==debyea) saveBinFile(file);
                            else warnMsg(" .diffb  format works only for Debye mode");

                        }
						else
                            if(file.rfind(".diff")!=string::npos){
                                if(mode==debyea) saveDiffFile(file);
                                else warnMsg(" .diff format works only for Debye mode");
                            }
							else
                                if(file.rfind(".laueb")!=string::npos){
                                    if(mode==laue) saveLaueFileBin(file);
                                    else warnMsg(" .laueb format works only for Laue mode");
                                }
                                else
                                    if(file.rfind(".laue")!=string::npos){
                                        if(mode==laue) saveLaueFile(file);
                                        else warnMsg(" .laue format works only for Laue mode");
                                    }
                                    else
                                        cerr<<"ERROR: Unknown file format"<<endl;
            }

}
//===================================================================================
void Cdiff::saveDatFile(const string &fileName)
{
fstream file(fileName,ios::out);

            if(!file){
                cerr<<"ERROR: couldn't open a file: "<<fileName<<endl;
            throw Cdiff::ERR_FILEOPEN;
            }

const double k0=180*M_1_PI;

            for(size_t i=0;i<t.size;i++)
                file<<t[i]*k0<<"\t"<<q[i]<<"\t"<<I[i]<<"\t"<<S[i]<<endl;


            file.close();

}
//===================================================================================
void Cdiff::saveDiffFile(const string &fileName)
{
fstream file(fileName,ios::out);

            if(!file){
                cerr<<"ERROR: couldn't open a file: "<<fileName<<endl;
            throw Cdiff::ERR_FILEOPEN;
            }

const double k0=180*M_1_PI;
std::time_t datetime = std::time(nullptr);
const size_t numOfAtoms=pdh->grain->atoms.size();

//ostream& endlw(ostream&  o)

auto (*eol)(ostream &)=(u2dendl)? &endld : &endlu;


            file<<"#ver: 04"<<eol;
            file<<"#title: diffraction data"<<eol;
            file<<"#date: "<<std::asctime(std::localtime(&datetime));
            file<<"#sizeRC: "<<t.size<<"\t4"<<eol;
            file<<"#numOfAtoms: "<<numOfAtoms<<eol;
            file<<(*this);            

csize nprec=10;
csize colwh=12;

            /// wypisuje nazwy column
            file<<"#"<<setw(colwh)<<" Theta"<<setw(colwh+1)<<" Q"<<setw(colwh+1)<<" I"<<setw(colwh+1  )<<" S"<<eol;
            /// koniec naglowka zaznaczony #
            file.fill('#');
            file<<setw(4*(1+colwh))<<'#'<<eol;
            file.fill(' ');

            for(size_t i=0;i<t.size;i++)
                file<<" "<<setprecision(nprec)<<setw(colwh)<<t[i]*k0
                    <<" "<<setprecision(nprec)<<setw(colwh)<<q[i]
                    <<" "<<setprecision(nprec)<<setw(colwh)<<I[i]
                    <<" "<<setprecision(nprec)<<setw(colwh)<<S[i]<<eol;

            file.close();
}
//===================================================================================
void Cdiff::saveBinFile(const string &fileName)
{
fstream file(fileName,ios::out | ios::binary);

            if(!file){
                cerr<<"ERROR: couldn't open a file: "<<fileName<<endl;
            throw Cdiff::ERR_FILEOPEN;
            }

const string version("v01");
const string data("DATA");
const size_t sopt=(size_t) saveopt.dataToSave.to_ulong();
const size_t dataSize=t.size;
const size_t numOfAtoms=pdh->grain->atoms.size();

            file.write( (char *) version.c_str(), version.length()*sizeof(char));
            file.write( (char *) &dataSize,       sizeof(size_t));
            file.write( (char *) &sopt,           sizeof(size_t));
            file.write( (char *) &numOfAtoms,     sizeof(size_t));
            file.write( (char *) data.c_str(),    data.length()*sizeof(char) );


vector<dataDiff * >  dataToSave;

            dataToSave.reserve(4);

            if(saveopt.dataToSave[Cdiff::StSaveOpt::BT])
                dataToSave.push_back(&t);
            if(saveopt.dataToSave[Cdiff::StSaveOpt::BQ])
                dataToSave.push_back(&q);
            if(saveopt.dataToSave[Cdiff::StSaveOpt::BI])
                dataToSave.push_back(&I);
            if(saveopt.dataToSave[Cdiff::StSaveOpt::BS])
                dataToSave.push_back(&S);


            dataToSave.shrink_to_fit();


            for(size_t i=0;i<dataSize;i++){
                for(dataDiff * ptr:dataToSave)
                    file.write( (char * )  & (*ptr)[i], sizeof( position) );
            }



            file.close();

}
//===================================================================================
void Cdiff::saveLaueFile(const string &fileName)
{
fstream file(fileName,ios::out);

            if(!file){
                cerr<<"ERROR: couldn't open a file: "<<fileName<<endl;
            throw Cdiff::ERR_FILEOPEN;
            }

//const double k0=180*M_1_PI;
std::time_t datetime = std::time(nullptr);
const size_t numOfAtoms=pdh->grain->atoms.size();

//ostream& endlw(ostream&  o)

auto (*eol)(ostream &)=(u2dendl)? &endld : &endlu;

const size_t dataSize=laueDiffData.size();
constexpr int nprec=10;
constexpr int colwh=15;
size_t i,j,k;

            file<<"#ver: 01"<<eol;
            file<<"#title: laue diffraction data"<<eol;
            file<<"#date: "<<std::asctime(std::localtime(&datetime));
            file<<"#sizeRC: "<<setprecision(nprec)<<setw(colwh)<<cub(dataSize)<<"\t4"<<eol;
            file<<"#sizeK: "<<dataSize<<eol;
            file<<"#numOfAtoms: "<<numOfAtoms<<eol;
            file<<(*this);
            file<<"#"<<" "<<setw(colwh-2)<<"kx"
                     <<" "<<setw(colwh)<<"ky"
                     <<" "<<setw(colwh)<<"kz"
                     <<" "<<setw(colwh)<<"I"<<endl;


            if(dataSize<10){
                for(k=0;k<dataSize;k++)
                    for(j=0;j<dataSize;j++)
                        for(i=0;i<dataSize;i++){
                            file<<" "<<setprecision(nprec)<<setw(colwh)<<kpoints[i]
                                <<" "<<setprecision(nprec)<<setw(colwh)<<kpoints[j]
                                <<" "<<setprecision(nprec)<<setw(colwh)<<kpoints[k]
                                <<" "<<setprecision(nprec)<<setw(colwh)<<std::abs(laueDiffData[i][j][k])<<endl;
                        }
            }
            else{
            CProgress progress;

                    progress.title=(std::string(" saving data "));
                    progress.start(dataSize);

                    for(k=0;k<dataSize;k++,progress++)
                        for(j=0;j<dataSize;j++)
                            for(i=0;i<dataSize;i++){
                                file<<" "<<setprecision(nprec)<<setw(colwh)<<kpoints[i]
                                    <<" "<<setprecision(nprec)<<setw(colwh)<<kpoints[j]
                                    <<" "<<setprecision(nprec)<<setw(colwh)<<kpoints[k]
                                    <<" "<<setprecision(nprec)<<setw(colwh)<<std::abs(laueDiffData[i][j][k])<<endl;
                            }

                    progress.stop();

            }

            file.close();
}
//===================================================================================
void Cdiff::saveLaueFileBin(const string &fileName)
{
 fstream file(fileName,ios::out | ios::binary);

             if(!file){
                 cerr<<"ERROR: couldn't open a file: "<<fileName<<endl;
             throw Cdiff::ERR_FILEOPEN;
             }

 const string version("v02");
 const string dataTag("DATA");
 const size_t dataSize=cub(laueDiffData.size());
 const size_t ksize=laueDiffData.size();
 vector<string> ststst(split<string>(range," "));//start step stop
 const double start=std::stod(ststst[0]);
 const double step =std::stod(ststst[1]);
 const double stop =std::stod(ststst[2]);
 const size_t numOfAtoms=pdh->grain->atoms.size();
 const position lambda__=std::stod(lambda);
 const char radType=radiationToChar();
 const char axis='X'+(size_t) saveopt.axis;
 const char format=(saveopt.format==StSaveOpt::SHORT) ? 'S' : 'L';

 int i,j,k;
 int *ptr_0, *ptr_1, *ptr_2;
 position I;

 const auto sizeofPosition=sizeof(position);
 CProgress progress;


            if(saveopt.axis==StSaveOpt::X){
                ptr_0=&i;
                ptr_1=&k;
                ptr_2=&j;
            }
            else{
                if(saveopt.axis==StSaveOpt::Y){
                    ptr_0=&j;
                    ptr_1=&i;
                    ptr_2=&k;
                }
                else{
                    ptr_0=&k;
                    ptr_1=&j;
                    ptr_2=&i;
                }
            }


            progress.title=(std::string(" saving data "));
            progress.start(dataSize);


            file.write( (char *) version.c_str(), version.length()*sizeof(char));
            file.write( (char *) &ksize,          sizeof(size_t));
            file.write( (char *) &start,          sizeofPosition);
            file.write( (char *) &step,           sizeofPosition);
            file.write( (char *) &stop,           sizeofPosition);
            file.write( (char *) &dataSize,       sizeof(size_t));
            file.write( (char *) &numOfAtoms,     sizeof(size_t));
            file.write( (char *) &radType,        sizeof(char));
            file.write( (char *) &axis,           sizeof(char));
            file.write( (char *) &format,         sizeof(char));
            file.write( (char *) &lambda__,       sizeof(double));
            file.write( (char *) &saveopt.layer,  sizeof(int));
            file.write( (char *) dataTag.c_str(), dataTag.length()*sizeof(char) );



            if(format=='S'){

                    for(*ptr_0=0;*ptr_0<ksize;(*ptr_0)++,progress++)
                        for(*ptr_1=0;*ptr_1<ksize;(*ptr_1)++)
                            for(*ptr_2=0;*ptr_2<ksize;(*ptr_2)++){


                                I=std::abs(laueDiffData[i][j][k]);
                                file.write( (char *)  &I,  sizeofPosition);
                            }
            }
            else{
                for(*ptr_0=0;*ptr_0<ksize;(*ptr_0)++,progress++)
                    for(*ptr_1=0;*ptr_1<ksize;(*ptr_1)++)
                        for(*ptr_2=0;*ptr_2<ksize;(*ptr_2)++){
                                I=std::abs(laueDiffData[i][j][k]);

                                file.write( (char *)  &kpoints[i], sizeofPosition);
                                file.write( (char *)  &kpoints[j], sizeofPosition);
                                file.write( (char *)  &kpoints[k], sizeofPosition);
                                file.write( (char *)  &I,          sizeofPosition);
                            }
            }


             file.close();
             progress.stop();


}
//===================================================================================
void Cdiff::openFile()
{
        if(fileNameIn.find(".diff")!=string::npos){
            openDiffFile();
        return;
        }


        cerr<<" unknown format "<<endl;
        throw Ediffstatus::ERR_FFORMAT;
}
//===================================================================================
void Cdiff::openDiffFile()
{
fstream fin(fileNameIn,ios::in);

            if(!fin){
                cerr<<"couldn't open file for reading"<<endl;
                fin.close();
            throw Ediffstatus::ERR_FOPEN;
            }

size_t size,headerSize,numOfAtoms;
std::string tagkey,tagvalue,fline;

            fin>>tagkey>>tagvalue;

            if(tagkey!="#ver:") throw Ediffstatus::ERR_FFORMAT;
            if(tagvalue=="04") headerSize=17-1; else throw Ediffstatus::ERR_FFORMAT;

            while(fin.get()!='\n' && !fin.eof()) ;


            ///----------- read header --------------------
            for(size_t hline=0;hline<headerSize;hline++){
                std::getline(fin,fline);

                if(hline>7) continue;

            vector<string> tokens(split<string>(fline," \t"));

                if(tokens[0]=="#sizeRC:") {
                    size=std::stoi(tokens[1]);
                continue;
                }

                if(tokens[0]=="#numOfAtoms:"){
                    numOfAtoms=std::stoi(tokens[1]);
                    pdh->grain->atoms.resize(numOfAtoms);
                continue;
                }

                if(tokens[0]=="#radiation:"){
                    radiation=(tokens[1]=="xray") ? ERADIATION::XRAY :  ( (tokens[1]=="neutron") ? ERADIATION::NEUTR :ERADIATION::ELECT);
                continue;
                }

                if(tokens[0]=="#lambda:"){
                    lambda=tokens[1];
                continue;
                }
            }
            ///-------------------------------------------------
            ///

position ft,fq,fI,fS;

            t.allocMem(size);
            q.allocMem(size);
            I.allocMem(size);
            S.allocMem(size);

            for(size_t i=0;i<size;i++){
                fin>>ft>>fq>>fI>>fS;
                t[i]=ft;
                q[i]=fq;
                I[i]=fI;
                S[i]=fS;
            }

            fin.close();
}
//===================================================================================
void Cdiff::allocMem(const size_t &size)
{

        t.allocMem(size);
        q.allocMem(size);
        I.allocMem(size);
        S.allocMem(size);

}
//===================================================================================
void Cdiff::clearMem()
{
        t.freeMem();
        q.freeMem();
        I.freeMem();
        S.freeMem();
}
//===================================================================================


bool Cdiff::loadScattFactors(const string &aname, const int index)
{
ifstream fileSc;
string elemName('['+aname+']');
string str;


            if(!fileSft.empty())
                fileSc.open(fileSft);

            if(!fileSc.is_open()) //current directory
                fileSc.open("scFact.sft");

            if(!fileSc.is_open()){  //home directory
            const char *homeDir=std::getenv("NPCLPATH");

                    if(!homeDir) {cerr<<"WARNING: empty homeDir"<<endl; return false;}

            const string shomeDir(homeDir);
            const string fileName(shomeDir+"/scFact.sft");

                fileSc.open(fileName);
            }

            if(!fileSc){
                fileSc.close();
                cerr<<"ERROR: file scFact.sft is missing"<<elemName<<endl;
            return false;
            }


            for(auto &c :elemName)
                c=std::toupper(c);

            do{
                fileSc>>str;
            }while(str.find(elemName)==string::npos  && !fileSc.eof() );


            if(fileSc.eof()) {
                fileSc.close();
                cerr<<"ERROR: file scFact.sft, missing element "<<elemName<<endl;
            return false;
            }


vector<string> tokens;
size_t pos,numOfcoeff=0;
const string coeff("a1b1a2b2a3b3a4b4cd1e1d2e2d3e3d4e4d5e5");
bool retval=false,ncs=false;

            try{

                scfactors[index].a1=-1;
                scfactors[index].ncs=-1;

                do{
                    fileSc>>str;
                  tokens=split<string>(str,"=");

                    if(tokens[0]=="ncs"){
                        scfactors[index].sf_array[9]=std::stof(tokens[1]);
                        ncs=true;
                    }
                    else{
                        pos=coeff.find(tokens[0]);

                        if(pos==string::npos) //ignore if not a scattering factor
                            continue;
                        else
                            numOfcoeff++;

                        if(pos<17)
                            pos>>=1; // podziel przez 2
                        else{
                            pos+=3;
                            pos>>=1;
                        }

                        scfactors[index].sf_array[pos]=std::stof(tokens[1]);
                    }

                }while ( !regex_match(str,std::regex("[:s:]*\\[\\w*\\][:s:]*")) &
                         !fileSc.eof() & !(numOfcoeff==19 && ncs));

               // cout<<"sf   "<<elemName<<endl;
               //for (float sf : scfactors[index].sf_array)
                //    cout<<"\t"<<setprecision(9)<<sf;

               //cout<<endl;


                retval=true;
            }
            catch (const std::invalid_argument& ia) {
                std::cerr << "\nERROR: Invalid argument: " << ia.what() <<":  "<<tokens[1]<<endl;
            }

            fileSc.close();

            return retval;
}
//===================================================================================
double Cdiff::xrayScatEqu(const double &q, const size_t index)
{
const double qq=q*q;
const Uscattfactors  * const ptr_scf=&scfactors[index];
const double suma=  ptr_scf->a1*exp(-ptr_scf->b1*qq)+
                    ptr_scf->a2*exp(-ptr_scf->b2*qq)+
                    ptr_scf->a3*exp(-ptr_scf->b3*qq)+
                    ptr_scf->a4*exp(-ptr_scf->b4*qq)+
                    ptr_scf->c;

return suma;
}

//===================================================================================
double Cdiff::electronScatEqu(const double &q, const size_t index)
{
const double qq=q*q;
const Uscattfactors  * const ptr_scf=&scfactors[index];
const double suma=  ptr_scf->d1*exp(-ptr_scf->e1*qq)+
                    ptr_scf->d2*exp(-ptr_scf->e2*qq)+
                    ptr_scf->d3*exp(-ptr_scf->e3*qq)+
                    ptr_scf->d4*exp(-ptr_scf->e4*qq)+
                    ptr_scf->d5*exp(-ptr_scf->e5*qq);

return suma;
}
//===================================================================================
double Cdiff::neutronCrossSec(const double &q, const size_t index)
{
        (void)q;
return scfactors[index].ncs;

}
//===================================================================================
void Cdiff::polarizationI()
{
auto sqrd=[](const double &x){return x*x;};
double pf;

            for(size_t i=0;i<I.size;i++){
                pf=(0.5+0.5*sqrd(std::cos(t.values[i])));
                I.values[i]*=pf;
                S.values[i]*=pf;
            }

}
//===================================================================================
void Cdiff::debyeaWaller()
{
const double u2=std::stod(dw);
constexpr double k3=1.0/3.0;
double fT;

            //phys. rev. b 73, 184113 2006, p. 12
            /// exp(-2M)=exp(-Q^2*u2/3))
            ///

            for(size_t i=0;i<I.size;i++){
                fT=std::exp(-u2*sqrd(q.values[i])*k3);
                I.values[i]*=fT;
                S.values[i]*=fT;
            }

}

//===================================================================================
void Cdiff::noiseSQ()
{
vector<string> params(split<string>(noise," "));

const double a=std::stod(params[1]);
const double b=std::stod(params[2]);


Crandom *fdist=(params[0]=="uniform") ? (Crandom *)  new CrandomUni(a,b) :
              ((params[0]=="norm")    ? (Crandom *)  new CrandomNor(a,b)
                                      : (Crandom *)  new CrandomLogNor(a,b) );

        for(size_t i=0;i<S.size;i++)
            S[i]+=fdist->randNumber();



        delete fdist;
}
//===================================================================================
/*void Cdiff::noise()
{

}*/
//===================================================================================
// normalization gives area =1
void Cdiff::normIS()
{

double sumI=0;
double sumS=0;


            for(size_t i=0;i<I.size;i++){
                sumI+=I[i];
                sumS+=S[i];
            }

            if(sumI>0 && sumS>0){
            const double invSumI=1.0/sumI;
            const double invSumS=1.0/sumS;


                    for(size_t i=0;i<I.size;i++){
                        I.values[i]*=invSumI;
                        S.values[i]*=invSumS;
                    }

            }
            else{
                cerr<<"ERROR: normalization failure"<<endl;
            }

}
//===================================================================================
void Cdiff::dispParams()
{
        cout<<"*** diffraction parameters ***"<<endl;
        cout<<(*this);

        for(string & fn: fileNameOut)
            cout<<"# output file : "<<fn<<endl;
}
//===================================================================================
void Cdiff::plot()
{
std::smatch matchPrm;

            if(DB){cout<<__FILE__<<":"<<__LINE__<<"data plotting"<<endl; }

            if(!std::regex_search(plotPrm,matchPrm,std::regex("[[:s:]]+u[[:s:]]+[12][[:s:]]*:[[:s:]]*[34][[:s:]]+"))){
                cout<<endl;
                errMsg(" 'plot' command doesn't recognize data type to be read; usage, e.g.: plot u 2:4 <other parameters>. ");
            return;
            }


string gplotPrm("u 1:2 "+plotPrm.substr(matchPrm.length()));
vector<string> tokens(split<string>(matchPrm.str(),"u :"));
position *x,*y;
string xlabel("set xlabel "),ylabel("set ylabel ");

                if(tokens[0]=="1"){
                    x=t.values;
                    xlabel+=("\"2*theta\"\n");
                }
                else{
                    x=q.values;
                    xlabel+=("\"Q\"\n");
                }

                if(tokens[1]=="3"){
                    y=I.values;
                    ylabel+="\"Intensity\"\n";
                }
                else{
                    y=S.values;
                    ylabel+="\"S\"\n";
                }


#ifdef __linux__
FILE *plotHandle = NULL;

                if(plotHandle == NULL){
                  plotHandle = popen("gnuplot -persist 2>/dev/null ", "w");
                  assert(plotHandle != NULL);
                }

                fprintf(plotHandle,xlabel.c_str());
                fprintf(plotHandle,ylabel.c_str());
                fprintf(plotHandle,"set grid\n");

                if(gplotPrm.find("title")==string::npos)
                    gplotPrm+=" notitle";


const size_t dataSize=t.size;
const std::string sizeType{(sizeof(position)==sizeof(double))? "double" : "float"};
const std::string plotcmd("plot '-' binary format=\"%%2"+sizeType+"\" record="+
                                                                std::to_string(dataSize)+" endian=little "+gplotPrm+" \n");

                 fprintf(plotHandle,plotcmd.c_str());

size_t i,j;
position *plotData=new position[2*dataSize];

               for(i=0,j=0;i<dataSize;i++,x++,y++){
                   plotData[j++]=(*x);
                   plotData[j++]=(*y);
               }

               fwrite(plotData,sizeof(position),2*dataSize,plotHandle);

               //fflush(plotHandle);
               pclose(plotHandle);

                delete [] plotData;
#endif									  
}
//===================================================================================
char Cdiff::radiationToChar()
{
    switch (radiation){
    case OFF :  return 'O';
    case NEUTR: return 'N';
    case ELECT: return 'E';
    case XRAY : return 'X';
    }
}
//===================================================================================
ostream & operator<<(ostream &o, Cdiff &cdiff)
{
std::map<std::string, std::string>::const_iterator cathLambda=Elements::xrad.find(cdiff.lambda);

        o<<"#range: "<<cdiff.range<<endl;
        o<<"#radiation: "<<((cdiff.radiation==Cdiff::XRAY)?"xray":
                                                       (cdiff.radiation==Cdiff::NEUTR) ? "neutron" :
                                                                                         (cdiff.radiation==Cdiff::ELECT) ?"electron" : "nt")<<endl;

        if(cathLambda!=Elements::xrad.end())
            o<<"#lambda: "<<cathLambda->second<<"\t"<<cathLambda->first<<endl;
        else
            o<<"#lambda: "<<cdiff.lambda<<endl;

        o<<"#extrapolate: "<<cdiff.extrapolate<<endl;
        o<<"#threads: "<<cdiff.threads<<endl;        
        o<<"#fastsinc: "<<cdiff.fastsinc<<endl;
        o<<"#polarization: "<<cdiff.polarization<<endl;
        o<<"#normalization: "<<cdiff.norm<<endl;
        o<<"#numOfatoms: "<<cdiff.pdh->grain->atomNamesNumber[0]<<endl;
        o<<"#dw: "<< ( (cdiff.dw.empty()) ? "0" : cdiff.dw)<<endl;


        if(!cdiff.comment.empty())
            o<<"#"<<cdiff.comment<<endl;

return o;
}


//===================================================================================
/*
void Cdiff::binLattice()
{

            scfactors=new Uscattfactors[2];

            if( !loadScattFactors(pdh->grain->atomNames[0],0)  ||
                !loadScattFactors(pdh->grain->atomNames[1],1)  ){
            throw Ediffstatus::ERR_SCATCOEFF;
            }

            if(pdh->grain->atomNamesNumber.empty()){
                cerr<<"ERROR: number of atoms for each type is not given"<<endl;
            throw Ediffstatus::ERR_ATYPE_NUM;
            }

vector<string> ststst(split<string>(range," "));//start step stop
const double start=std::stod(ststst[0]);
const double step =std::stod(ststst[1]);
const double stop =std::stod(ststst[2]);
const int extraSize=(extrapolate) ? start/step : 0;
const int thetaSize=1+std::floor( (stop-start)/step);

            if(thetaSize<2 || extraSize<0){
                cerr<<"ERROR: wrong values of 'range' directive"<<endl;
            throw Ediffstatus::ERR_RANGE;
            }

const size_t totSize=thetaSize+extraSize;
position *binX= pdh->dataX.values;
position *binYii= pdh->dataYii.values;
position *binYij= pdh->dataYij.values;
position *binYjj= pdh->dataYjj.values;
const size_t binSize=pdh->dataSize();
const size_t N=pdh->grain->atoms.size();
const position N4=N*0.25;

            allocMem(totSize);

position * const dataX=t.values;
position * const dataQ=q.values;
position * const dataI=I.values;
position * const dataS=S.values;

size_t noneZeroBinsSizeII=0;
size_t noneZeroBinsSizeIJ=0;
size_t noneZeroBinsSizeJJ=0;
size_t i,j;

            omp_set_num_threads(std::stoi(threads));

            #pragma omp parallel sections shared(noneZeroBinsSizeII,noneZeroBinsSizeIJ,noneZeroBinsSizeJJ) private(i)
            {

                #pragma omp section
                {
                    for(i=0;i<binSize;i++)
                        if(binYii[i]>0) noneZeroBinsSizeII++;
                }


                #pragma omp section
                {
                    for(i=0;i<binSize;i++)
                        if(binYij[i]>0)noneZeroBinsSizeIJ++;
                }


                #pragma omp section
                {
                    for(i=0;i<binSize;i++)
                        if(binYjj[i]>0) noneZeroBinsSizeJJ++;
                }
            }

size_t *noneZeroBinsII=nullptr;
size_t *noneZeroBinsIJ=nullptr;
size_t *noneZeroBinsJJ=nullptr;


            if(noneZeroBinsSizeII && noneZeroBinsSizeIJ && noneZeroBinsSizeJJ){
                noneZeroBinsII=new size_t[noneZeroBinsSizeII];
                noneZeroBinsIJ=new size_t[noneZeroBinsSizeIJ];
                noneZeroBinsJJ=new size_t[noneZeroBinsSizeJJ];
            }
            else{
                cerr<<" size II ,IJ, JJ =0"<<endl;
            return;
            }

            #pragma omp parallel sections shared(noneZeroBinsII,noneZeroBinsIJ) private(i,j)
            {

                #pragma omp section
                {
                    for(i=0,j=0;i<binSize;i++)
                        if(binYii[i]>0) noneZeroBinsII[j++]=i;
                }

                #pragma omp section
                {

                    for(i=0,j=0;i<binSize;i++)
                        if(binYij[i]>0) noneZeroBinsIJ[j++]=i;
                }

                #pragma omp section
                {

                    for(i=0,j=0;i<binSize;i++)
                        if(binYjj[i]>0) noneZeroBinsJJ[j++]=i;
                }
            }

                ////// 2. setting initial values

const position lambda__=std::stod(lambda);
const position k1=M_PI/360;
const position invWL=1.0/lambda__;
const position k2=4*M_PI*invWL;
const position startAlpha=(start-extraSize*step)*k1;
const position stepAlpha=step*k1;
position alpha,th,lsin,qx;

                ////// 3. extrapolate

                   for(int i=0;i<extraSize;i++){
                        alpha=step*i;
                        th=alpha*k1;
                        lsin=sin(th);
                        qx=k2*lsin;
                        dataX[i]=2*th;  // radians
                        dataQ[i]=qx;     // Q
                        dataI[i]=0;
                        dataS[i]=0;
                    }

                /////// 4. MAIN TASK

position fi,fj, fi2,fj2,fijN4,fiArg;
//auto sqr=[](const position &x){return x*x;};
position sumaII,sumaIJ,sumaJJ,sumaTot;

            runDIFF=true;
            progressDIFFMaxValueInv=100.0/thetaSize;
            progressDIFFValue=0;

            startTime = std::time(NULL);

            for(i=extraSize;i<totSize;i++,progressDIFFValue++){

                th=startAlpha+stepAlpha*i;
                lsin=sin(th);
                qx=k2*lsin;
                //
                fiArg=(lsin*invWL);

                fi=electronScatEqu(fiArg,0);
                fj=electronScatEqu(fiArg,1);
                fi2=fi*fi;
                fj2=fj*fj;

                dataX[i]=2*th;
                dataQ[i]=qx;

                sumaII=0;
                sumaIJ=0;
                sumaJJ=0;



                #pragma omp parallel
                {
                size_t ii;
                position qr_ii,sinc_ii;


                    #pragma omp for reduction(+:sumaII) nowait
                    for(ii=0;ii<noneZeroBinsSizeII;ii++){
                        qr_ii=qx*binX[noneZeroBinsII[ii]];
                        sinc_ii=sin(qr_ii)/qr_ii;
                        sumaII+=binYii[noneZeroBinsII[ii]]*sinc_ii;
                    }

                size_t ij;
                position qr_ij,sinc_ij;

                    #pragma omp for reduction(+:sumaIJ) nowait
                    for(ij=0;ij<noneZeroBinsSizeIJ;ij++){
                        qr_ij=qx*binX[noneZeroBinsIJ[ij]];
                        sinc_ij=sin(qr_ij)/qr_ij;
                        sumaIJ+=binYij[noneZeroBinsIJ[ij]]*sinc_ij;
                    }


                 size_t jj;
                 position qr_jj,sinc_jj;

                     #pragma omp for reduction(+:sumaJJ) nowait
                     for(jj=0;jj<noneZeroBinsSizeJJ;jj++){
                         qr_jj=qx*binX[noneZeroBinsJJ[jj]];
                         sinc_jj=sin(qr_jj)/qr_jj;
                         sumaJJ+=binYjj[noneZeroBinsJJ[jj]]*sinc_jj;
                     }
                }


                sumaTot=(fi2*sumaII+fj2*sumaJJ+fi*fj*sumaIJ );
                fijN4=N4*(fi2+fj2);

                dataS[i]=1+sumaTot/fijN4;    // always true
                dataI[i]=2*(fijN4+sumaTot);

            }

            runDIFF=false;

            stopTime = std::time(NULL);

            delete [] noneZeroBinsII;
            delete [] noneZeroBinsIJ;
            delete [] noneZeroBinsJJ;

}

do{
                    fileSc>>str;
                    tokens=split<string>(str,"=");

                    if(tokens[0]=="ncs"){
                        scfactors[index].sf_array[9]=std::stof(tokens[1]);
                        numOfcoeff++;
                    }
                    else{
                        pos=coeff.find(tokens[0]);

                        if(pos==string::npos)
                            continue;
                        else
                            numOfcoeff++;

                        pos>>=1; // podziel przez 2
                        scfactors[index].sf_array[pos]=std::stof(tokens[1]);
                    }

                }while (numOfcoeff!=10 && !fileSc.eof() );//&& tokens[0][0]!='['


*/


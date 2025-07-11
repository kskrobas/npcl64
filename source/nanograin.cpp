/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * nanograin.cpp
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


#include "nanograin.h"
#include <cmath>
#include <random>
#include <chrono>
#include <iomanip>
#include <string>
#include <functional>
#include <omp.h>
#include <regex>
#include <algorithm>

#include "createdir.h"
#include "crandom.h"
#include "colormsg.h"
#include "parse_expr.h"
#include "cprogress.h"
#include "stsymmgroupgenerator.h"

#ifndef __linux
#define M_PI 3.1415926539
#endif

//#define DEBUG


#ifdef DEBUG
#define DB true
#else
#define DB false
#endif

#define FLINE __FILE__<<__LINE__


//-----------------------------------------------------------------------------
inline position sqr(const position &x){return x*x;}
inline position cube(const position &x){return x*x*x;}


cpos wi=0.5;//sin(30)
cpos wj=sqrt(3)/2;//sin(60)
cpos k2s3=2.0*sqrt(3);


typedef const double cdouble;
//-----------------------------------------------------------------------------

size_t NanoGrain::StAtom::idIterator=0;
list<size_t> NanoGrain::StNanoGrain::savedNumOfAtoms;


//-----------------------------------------------------------------------------
void gotoEOL(fstream &file)
{
        while(file.get()!='\n'  &&  !file.eof()) ;
}

//-----------------------------------------------------------------------------


ostream & NanoGrain::operator<<(ostream &o,const NanoGrain::StAtom &a)
{
    o<<a.id<<"\t"<<a.atype<<"\t"<<a.x<<"\t"<<a.y<<"\t"<<a.z;
return o;
}

//-----------------------------------------------------------------------------
string parseABC(const string &expression)
{
CParseABC parseABC;
string abcexpr;

            if(! parseABC.run(expression,abcexpr))
                throw NanoGrain::Status::ERR_ABCPARSE;

return abcexpr;
}
//-----------------------------------------------------------------------------


void NanoGrain::StNanoGrain::resetPrms()
{
    /*
        if(callbackSetThreads!=nullptr){
            (this->*callbackSetThreads)(threads);
            //((NanoGrain::StNanoGrain*)this)->NanoGrain::StNanoGrain::callbackSetThreads
        }
        else
            threads="2";
            */

        threads= (mthreads.empty()) ? "2" : mthreads;


        radius.clear();
        scaleFactors.clear();
        clp.clear();
        disloc.clear();
        rmatoms.clear();
        uc.clear();
        replicate.clear();
        rename.clear();
        tric.reset();
        fileCIF.clear();

        atomTypes.clear();
        atomNamesNumber.clear();

        shape.clear();
        shapePrm.clear();
        shapePrm2D.clear();
        delete ssShape;
        ssShape=nullptr;

        hcpABC.clear();
        hcpfillup.clear();
        hcpfault.clear();
        hcpcs="circle";
        hcpu="0.75";
        hcpFaultPos.clear();

        coreshell.prm.clear();

        fileNameIn.clear();
        fileNameOut.clear();
        fileNameHeader.clear();

        voidsprm.clear();

        atoms.clear();
        StAtom::resetId();

        atomDisperse.clear();
        comment.clear();
        source=AUTOCREATE;
        lmpstyle="atomic";
        faultmode=OFF;
        atomsRemoved.str()="0";

        center=COFF;
        catomType.clear();

        disperse=false;
        dispParams=false;
        hcpsl=true;
        hcpsurf=EHCPSURF::sB;
        numOfAtomsTest=false;
        numOfAtomsPush.clear();
        saveFileStatusPush.clear();
        margins.clear();
        rotatePrm.clear();

        saveopt.reset();
        vremoveAtomsPrms.clear();
}

//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::sc()
{

cpos fclp=getLP();
cpos fradius=getRadius(fclp);

cpos ucv=fclp*fclp*fclp; // unit cell volume
cpos maxLp=(cpos) floor(fradius/fclp);
cpos pSTART=-maxLp*fclp;
position X,Y,Z;
const size_t atomNameA=0;


            if(shape=="sphere"){
               // eshape=StNanoGrain::SPHERE;
            cpos  gv=4*M_PI*cube(fradius)/3.0; //spherer grain volume
            const size_t   appSize=(size_t)  gv/ucv;
            cpos fradii2=sqr(fradius);
            position radii2,X2,Y2,r2xy;

                atoms.reserve(2*appSize);

                for(X=pSTART;X<=fradius;X+=fclp){
                    X2=X*X;
                    for(Y=pSTART;Y<=fradius;Y+=fclp){
                        Y2=Y*Y;
                        r2xy=X2+Y2;
                        for(Z=pSTART;Z<=fradius;Z+=fclp){
                            radii2=Z*Z+r2xy;

                            if(radii2<=fradii2)
                                atoms.push_back(StAtom(X,Y,Z,atomNameA,radii2));
                        }
                    }
                }
            }
            else{
              //  eshape=StNanoGrain::CUBE;
            cpos csize=fradius*2+fclp;
            const size_t   appSize=(size_t)  cube(csize)/ucv;

                    atoms.reserve(appSize);

                    for(X=pSTART;X<=fradius;X+=fclp){
                        for(Y=pSTART;Y<=fradius;Y+=fclp)
                            for(Z=pSTART;Z<=fradius;Z+=fclp){
                                atoms.push_back(StAtom(X,Y,Z,atomNameA));
                            }
                    }
            }


            atoms.shrink_to_fit();
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::bcc()
{

cpos fclp=getLP();
cpos fradius=getRadius(fclp);
cpos ucv=fclp*fclp*fclp; // unit cell volume
cpos maxLp=(cpos) floor(fradius/fclp);
cpos pSTART=-maxLp*fclp;
position X,Y,Z,Xh,Yh,Zh;
position X2h,Y2h;
const size_t atomNameA=0;


            if(shape=="sphere"){
              //  eshape=StNanoGrain::SPHERE;
            cpos  gv=4*M_PI*cube(fradius)/3.0; //spherer grain volume
            const size_t   appSize=(size_t)  gv/ucv;
            cpos fradii2=sqr(fradius);
            position radii2,X2,Y2,r2xy,r2xyh;

                atoms.reserve(2*appSize);

                for(X=pSTART,Xh=pSTART+fclp*0.5;X<=fradius;X+=fclp,Xh+=fclp){

                    X2=X*X;
                    X2h=Xh*Xh;

                    for(Y=pSTART,Yh=pSTART+fclp*0.5;Y<=fradius;Y+=fclp,Yh+=fclp){

                        Y2=Y*Y;
                        Y2h=Yh*Yh;

                        r2xy=X2+Y2;
                        r2xyh=X2h+Y2h;

                        for(Z=pSTART,Zh=pSTART+fclp*0.5;Z<=fradius;Z+=fclp,Zh+=fclp){
                            radii2=Z*Z+r2xy;

                            if(radii2<=fradii2){
                                atoms.push_back(StAtom(X,Y,Z,atomNameA,radii2));
                             //   atoms.push_back(StAtom(Xh,Yh,Zh,&atomNameA));
                            }

                            radii2=r2xyh+Zh*Zh;
                            if(radii2<=fradii2){
                             //   atoms.push_back(StAtom(X,Y,Z,&atomNameA,radii2));
                                atoms.push_back(StAtom(Xh,Yh,Zh,atomNameA));
                            }


                        }
                    }
                }
            }
            else{
                 //   eshape=StNanoGrain::CUBE;
                cpos csize=fradius*2+fclp;
                const size_t   appSize=(size_t)  cube(csize)/ucv;

                        atoms.reserve(2*appSize);

                        for(X=pSTART,Xh=pSTART+fclp*0.5;X<=fradius;X+=fclp,Xh+=fclp){
                            for(Y=pSTART,Yh=pSTART+fclp*0.5;Y<=fradius;Y+=fclp,Yh+=fclp)
                                for(Z=pSTART,Zh=pSTART+fclp*0.5;Z<=fradius;Z+=fclp,Zh+=fclp){
                                    atoms.push_back(StAtom(X,Y,Z,atomNameA));
                                    atoms.push_back(StAtom(Xh,Yh,Zh,atomNameA));
                                }
                        }
            }

            atoms.shrink_to_fit();
}
//-----------------------------------------------------------------------------
struct acceptRadii{
virtual bool isInside (const double &R) =0;
virtual ~acceptRadii(){ };
};
struct acceptRadiiOut: acceptRadii{
const double rout;
        acceptRadiiOut(cpos &Rout):rout(Rout){ }
        bool isInside (const double &R){ return R<=rout; }

};
struct acceptRadiiInOut: acceptRadii{
const double rout;
const double rin;
        acceptRadiiInOut(cpos &Rout,cpos &Rin):rout(Rout),rin(Rin){ }
        bool isInside (const double &R){ return R<=rout && R>=rin; }

};
//bool operator ==
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::fcc(EFCCTYPE fcctype)
{

cpos fclp=getLP();
cpos fradius=getRadius(fclp);

cpos ucv=fclp*fclp*fclp; // unit cell volume
cpos maxLp=(cpos) (floor(fradius/fclp));
cpos pSTART=-maxLp*fclp;
cpos pSTOP=-pSTART;
position X,Y,Z,Xh,Yh,Zh;
position X2h,Y2h;
const size_t atomNameA=0;

            if(shape=="sphere"){
              //  eshape=StNanoGrain::SPHERE;
            cpos  gv=4*M_PI*cube(fradius)/3.0; //spherer grain volume
            size_t   appSize=(size_t)  gv/ucv;
            const size_t mult=(fcctype) ? 8 : 4;
            cpos fradii2=sqr(fradius);
            position radii2,X2,Y2,r2xy,r2xyh;                                  
            acceptRadii *atomInside;

                    if(cradius.empty())  //cavity radius
                        atomInside=new acceptRadiiOut(fradii2);
                    else{
                    cpos cradii2{sqr(getRadius(fclp,cradius))};
                        atomInside=new acceptRadiiInOut(fradii2,cradii2);
                    }

                    atoms.reserve(mult*appSize);

                    for(X=pSTART,Xh=pSTART+fclp*0.5;X<=pSTOP;X+=fclp,Xh+=fclp){

                        X2=X*X;
                        X2h=Xh*Xh;

                        for(Y=pSTART,Yh=pSTART+fclp*0.5;Y<=pSTOP;Y+=fclp,Yh+=fclp){

                            Y2=Y*Y;
                            Y2h=Yh*Yh;

                            r2xy=X2+Y2;
                            r2xyh=X2h+Y2h;

                            for(Z=pSTART,Zh=pSTART+fclp*0.5;Z<=pSTOP;Z+=fclp,Zh+=fclp){
                                radii2=Z*Z+r2xy;

                                //if(radii2<=fradii2)
                                if( atomInside->isInside(radii2))
                                    atoms.push_back(StAtom(X,Y,Z,atomNameA,radii2));

                                radii2=r2xyh+Z*Z;
                                //if(radii2<=fradii2)
                                if( atomInside->isInside(radii2))
                                    atoms.push_back(StAtom(Xh,Yh,Z,atomNameA));                                


                                radii2=X2h+Y2+Zh*Zh;
                                //if(radii2<=fradii2)
                                if( atomInside->isInside(radii2))
                                    atoms.push_back(StAtom(Xh,Y,Zh,atomNameA));

                                radii2=X2+Y2h+Zh*Zh;
                                //if(radii2<=fradii2)
                                if( atomInside->isInside(radii2))
                                    atoms.push_back(StAtom(X,Yh,Zh,atomNameA));


                            }
                        }
                    }

                    delete atomInside;
            }
            else{    /// cube
                if(shape=="cube"){
                cpos csize=fradius*2+fclp;
                size_t   appSize=(size_t) cube(csize)/ucv;
                const size_t mult=(fcctype) ? 8 : 4;

                        if(disloc.empty()){

                                atoms.reserve(mult*appSize);

                                for(X=pSTART,Xh=pSTART+fclp*0.5;X<=fradius;X+=fclp,Xh+=fclp){
                                    for(Y=pSTART,Yh=pSTART+fclp*0.5;Y<=fradius;Y+=fclp,Yh+=fclp)
                                        for(Z=pSTART,Zh=pSTART+fclp*0.5;Z<=fradius;Z+=fclp,Zh+=fclp){
                                            atoms.push_back(StAtom(X,Y,Z,atomNameA));
                                            atoms.push_back(StAtom(Xh,Yh,Z,atomNameA));
                                            atoms.push_back(StAtom(Xh,Y,Zh,atomNameA));
                                            atoms.push_back(StAtom(X,Yh,Zh,atomNameA));
                                        }
                                }
                        }
                        else { /// add dislocations
                        vector<string> dislocTokens(split<string>(disloc," "));
                        const double from2=sqr(std::stod(dislocTokens[2]));
                        const double to2=  sqr(std::stod(dislocTokens[3]));
                        const double prob=std::stod(dislocTokens[4]);
                        Crandom *distr= new CrandomUni(0,1);
                        double r2;

                            atoms.reserve((4+1)*appSize);

                            if(DB)cout<<__FILE__<<":"<<__LINE__<<"  "<<dislocTokens[1]<<endl;

                            if(dislocTokens[0]=="dumbell"){
                            const string mode((dislocTokens.size()==5)? "z" : dislocTokens[5]);
                            struct Stxyz{double x,y,z; Stxyz(cdouble &x__, cdouble &y__, cdouble &z__ ){x=x__;y=y__;z=z__; }};
                            auto modPos=(mode=="z")
                                                 ? std::function<Stxyz(cdouble &, cdouble &, cdouble &, cdouble )>{ [&fclp] (cdouble &x, cdouble &y, cdouble &z,cdouble op=1){ return Stxyz(x,y,z+op*fclp*0.25);  } }
                                   : (mode=="x") ? std::function<Stxyz(cdouble &, cdouble &, cdouble &, cdouble )>{ [&fclp] (cdouble &x, cdouble &y, cdouble &z,cdouble op=1){ return Stxyz(x+op*fclp*0.25,y,z);  } }
                                    /*mode=="y"*/: std::function<Stxyz(cdouble &, cdouble &, cdouble &, cdouble )>{ [&fclp] (cdouble &x, cdouble &y, cdouble &z,cdouble op=1){ return Stxyz(x,y+op*fclp*0.25,z);  } };

                                for(X=pSTART,Xh=pSTART+fclp*0.5;X<=fradius;X+=fclp,Xh+=fclp){
                                    for(Y=pSTART,Yh=pSTART+fclp*0.5;Y<=fradius;Y+=fclp,Yh+=fclp)
                                        for(Z=pSTART,Zh=pSTART+fclp*0.5;Z<=fradius;Z+=fclp,Zh+=fclp){
                                            atoms.push_back(StAtom(X,Y,Z,atomNameA));
                                            atoms.push_back(StAtom(Xh,Y,Zh,atomNameA));
                                            atoms.push_back(StAtom(X,Yh,Zh,atomNameA));

                                            r2=sqr(Xh)+sqr(Yh)+sqr(Z);

                                            if(from2<=r2 && r2<=to2){
                                               if(distr->randNumber()<=prob){
                                               Stxyz xyz(modPos(Xh,Yh,Z,1));

                                                   atoms.emplace_back(StAtom(xyz.x,xyz.y,xyz.z,atomNameA));
                                                   atoms.back().mark=true;

                                                   xyz=modPos(Xh,Yh,Z,-1);
                                                   //atoms.push_back(StAtom(Xh,Yh, Z-fclp*0.25,atomNameA));
                                                   atoms.emplace_back(StAtom(xyz.x,xyz.y,xyz.z,atomNameA));
                                                   atoms.back().mark=true;
                                                  // cout<<" yes \n";
                                               }
                                               else
                                                   atoms.push_back(StAtom(Xh,Yh,Z,atomNameA));
                                            }
                                            else {
                                                atoms.push_back(StAtom(Xh,Yh,Z,atomNameA));
                                            }
                                      }
                                }
                            }
                            else {

                                    for(X=pSTART,Xh=pSTART+fclp*0.5;X<=fradius;X+=fclp,Xh+=fclp){
                                        for(Y=pSTART,Yh=pSTART+fclp*0.5;Y<=fradius;Y+=fclp,Yh+=fclp)
                                            for(Z=pSTART,Zh=pSTART+fclp*0.5;Z<=fradius;Z+=fclp,Zh+=fclp){
                                                atoms.push_back(StAtom(X,Y,Z,atomNameA));
                                                atoms.push_back(StAtom(Xh,Y,Zh,atomNameA));
                                                atoms.push_back(StAtom(X,Yh,Zh,atomNameA));
                                                atoms.push_back(StAtom(Xh,Yh,Z,atomNameA));

                                                r2=sqr(Xh)+sqr(Yh)+sqr(Z);

                                                if(from2<=r2 && r2<=to2){
                                                    if(distr->randNumber()<=prob){
                                                    cdouble aveX=(X+Xh)*0.5;
                                                    cdouble aveY=(Y+Yh)*0.5;
                                                    cdouble aveZ=(Z+Zh)*0.5;

                                                        atoms.push_back(StAtom(aveX,aveY,aveZ,atomNameA));
                                                        atoms.back().mark=true;
                                                    }
                                                }
                                          }
                                    }
                            }


                            delete distr;
                        }
                }
                else{
                    if(shape=="cone"){
                     vector<string> coneParams(split<string>(shapePrm," "));//scale factors
                            if(coneParams[1]=="0"){
                                errMsg("cone height must be greater than 0");
                            throw Status::ERR_HEIGHT;
                            }

                     cpos minRadius=std::stod(coneParams[0]);
                     cpos coneHeight=std::stod(coneParams[1]);
                     cpos slope=(minRadius-fradius)/coneHeight;
                     cpos zshift=fradius+slope*coneHeight*0.5;
                     cpos  gv=1.0*M_PI*sqr(fradius)*coneHeight/3.0; //cone grain volume
                     size_t   appSize=(size_t)  gv/ucv;
                     const size_t mult=(fcctype) ? 8 : 4;
                     position radii2,X2,Y2,r2xy,fradii2,fradii2h;
                     cpos heightHalf=0.5*coneHeight;
                     cpos cmaxLp=(cpos) (floor( (fradius>minRadius) ?  fradius/fclp : minRadius/fclp));
                     cpos cSTART=-cmaxLp*fclp;
                     cpos cSTOP=-cSTART;

                            atoms.reserve(mult*appSize);

                            for(X=cSTART,Xh=cSTART+fclp*0.5;X<=cSTOP+fclp;X+=fclp,Xh+=fclp){

                                X2=X*X;
                                X2h=Xh*Xh;

                                for(Y=cSTART,Yh=cSTART+fclp*0.5;Y<=cSTOP+fclp;Y+=fclp,Yh+=fclp){

                                    Y2=Y*Y;
                                    Y2h=Yh*Yh;

                                    r2xy=X2+Y2;

                                    for(Z=-heightHalf,Zh=-heightHalf+fclp*0.5;Z<=heightHalf;Z+=fclp,Zh+=fclp){
                                        radii2=r2xy;
                                        fradii2= sqr(slope*Z + zshift );
                                        //fradii2h=sqr(slope*Zh+ zshift );

                                        if(radii2<=fradii2){
                                            atoms.push_back(StAtom(X,Y,Z,atomNameA,radii2));
                                            //if(r2xyh<=fradii2)
                                            atoms.push_back(StAtom(Xh,Yh,Z,atomNameA));
                                            //if(Xh<cSTOP)
                                            atoms.push_back(StAtom(Xh,Y,Zh,atomNameA));
                                            //if(Yh<cSTOP)
                                            atoms.push_back(StAtom(X,Yh,Zh,atomNameA));
                                        }
                                    }
                                }
                            }
                    }
                    else{
                        errMsg(" unknown shape type ");
                    throw Status::ERR_UNKSHAPE;
                    }
                }
            }

            if(fcctype==EFCCTYPE::ZB){//zinc blende structure
            cpos fclpQ=fclp*0.25;
            const size_t csize=atoms.size();

            StAtom atomB;

                   // atomB.atype=(atomNameB.empty()) ? atomNameA  : &atomNameB;
                    atomB.atype=(atomTypes.size()==1) ? atomNameA  : 1;

                    for(size_t i=0;i<csize;i++){

                        atomB.x=atoms[i].x+fclpQ;
                        atomB.y=atoms[i].y+fclpQ;
                        atomB.z=atoms[i].z+fclpQ;
                        atomB.calcR2();

                        atoms.push_back(atomB);

                    }
            }


            if(fcctype==EFCCTYPE::UO2){
            cpos fclpQ=fclp*0.5;
            StAtom atomB;

                   // atomB.atype=(atomNameB.empty()) ? atomNameA  : &atomNameB;
                    atomB.atype=(atomTypes.size()==1) ? atomNameA  : 1;

                    for(X=pSTART,Xh=pSTART+fclp*0.25;X<=fradius;X+=fclp,Xh+=fclp){
                        for(Y=pSTART,Yh=pSTART+fclp*0.25;Y<=fradius;Y+=fclp,Yh+=fclp)
                            for(Z=pSTART,Zh=pSTART+fclp*0.25;Z<=fradius;Z+=fclp,Zh+=fclp){

                                atoms.push_back(StAtom(Xh,Yh,Zh,  atomB.atype));
                                atoms.push_back(StAtom(Xh,Yh,fclpQ+Zh,atomB.atype));

                                atoms.push_back(StAtom(Xh,fclpQ+Yh,Zh,  atomB.atype));
                                atoms.push_back(StAtom(Xh,fclpQ+Yh,fclpQ+Zh,atomB.atype));

                                atoms.push_back(StAtom(fclpQ+Xh,fclpQ+Yh,Zh,  atomB.atype));
                                atoms.push_back(StAtom(fclpQ+Xh,fclpQ+Yh,fclpQ+Zh,atomB.atype));

                                atoms.push_back(StAtom(fclpQ+Xh,Yh,Zh,  atomB.atype));
                                atoms.push_back(StAtom(fclpQ+Xh,Yh,fclpQ+Zh,atomB.atype));

                            }
                    }

            }

            atoms.shrink_to_fit();
}
//-----------------------------------------------------------------------------

double shapePrmValue(const string  &shapePrm)
{
            if(shapePrm.find("rand")!=string::npos){
            vector<string> randPrms(split<string>(shapePrm,"(,)"));
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator (seed);
            std::uniform_real_distribution<double> distribution(std::stod(randPrms[1]),std::stod(randPrms[2]));

            return distribution(generator);
            }
            else {
                return std::stod(shapePrm);
            }
}

//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::hcp()
{             
//size_t posLP=clp.find("lp");

cpos fclp=getLP();
cpos fradius=getRadius(fclp);

            hcpRadius=fradius;
            hcpRadius2=sqr(hcpRadius);
            hcpa=fclp/std::sqrt(2);
            hcpc=fclp/std::sqrt(3);

vector<string> hcpParams(split<string>(hcpABC," "));

            //if(DB) cout<<__FILE__<<":"<<__LINE__<<" set hcp parameters "<<endl;

            if(hcpParams.size()==4)
                buildRandomABC(hcpParams);

            if(hcpABC.empty() && hcpfillup.empty())
            throw NanoGrain::ERR_UNKSEQ;

            if(!hcpfillup.empty())
                buildHcpABC();

            if(faultmode!=StNanoGrain::OFF)
                insertFaults();

            ///cout<<"shapePrm :"<<shapePrm<<endl;


            if(shapePrm.empty()){  /// non supersphere shapes
                    if(shape=="cylinder"){
                        fbuildHcpShape=&StNanoGrain::buildHcpCylinder;
                        hcpCylinderHeight=hcpc*hcpABC.length();
                    }
                    else
                        fbuildHcpShape=&StNanoGrain::buildHcpSphere;
            }
            else{
            vector<string> shapeParams(split<string>(shapePrm," "));
            const double p=shapePrmValue(shapeParams[0]);
            const double R=0.75*hcpRadius;

                if(shape=="cubic")
                    if(shapeParams.size()==1) ssShape=new CCubic(p,R);
                    else ssShape=new CCubic(p,R,std::stod(shapeParams[1]),std::stod(shapeParams[2]),std::stod(shapeParams[3]));
                else{
                    if(shape=="oct")
                        ssShape=new COctahedral(p,R);
                    else{
                        if(shape=="dod")
                            ssShape=new CDodecahedral(p,R);
                        else {//poly
                        const double a=shapePrmValue(shapeParams[1]);
                        const double b=shapePrmValue(shapeParams[2]);

                            if(DB){
                                cout<<"shapeParams : ";
                                for (auto &value:shapeParams)
                                    cout<<"  "<<value;

                                cout<<endl;
                            }

                            if(shapeParams.size()==1) ssShape=new CPolyhedral(p,R,a,b);
                            else ssShape=new CPolyhedral(p,R,a,b,std::stod(shapeParams[3]),std::stod(shapeParams[4]),std::stod(shapeParams[5]));
                        }
                    }
                }

                fbuildHcpShape=&StNanoGrain::buildHcpSuperSphere;
            }


             if(hcpcs=="rhomb")
                buildHcpBaseRhm();
             else
                buildHcpBaseHex();

             buildHcpLayers();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::zb110()
{
cpos fclp=getLP();
cpos fradius=getRadius(fclp);

cpos a=fclp*sqrt(2);
cpos b=a*0.5;
cpos c=fclp;

position X,Y,Z,X2,X4,X34,Z4;
const size_t atomNameA=0;

// atomB.atype=(atomTypes.size()==1) ? atomNameA  : 1;


            if(shape=="cube"){
            const int aMax=static_cast<int>(fradius/a);
            const int bMax=static_cast<int>(fradius/b);
            const int cMax=static_cast<int>(fradius/c);
            const int numOfCells=aMax*bMax*cMax;


                    atoms.reserve(static_cast<size_t>(numOfCells*8));

                    for(int i=-aMax;i<aMax;i++){
                        X=i*a;
                        X2=X+a*0.5;
                        X34=X+a*0.75;
                        X4=X+a*0.25;

                        for(int j=-bMax;j<bMax;j++){
                            Y=j*b;

                            for(int k=-cMax;k<cMax;k++){
                                Z=k*c;
                                Z4=Z+c*0.25;
                                atoms.push_back(StAtom(Z,X,Y,atomNameA));
                                atoms.push_back(StAtom(Z,X2,Y,atomNameA));
                                atoms.push_back(StAtom(Z4,X4,Y,atomNameA));
                                atoms.push_back(StAtom(Z4,X34,Y,atomNameA));
                            }
                        }
                    }

                const size_t asize=atoms.size();
                StAtom atomB;

                        atomB.atype=(atomTypes.size()==1) ? atomNameA  : 1;

                    for(size_t i=0;i<asize;i++){
                        atomB.x=atoms[i].x+c*0.5;
                        atomB.y=atoms[i].y+a*0.25;
                        atomB.z=atoms[i].z+b*0.5;

                        atoms.push_back(atomB);
                    }


            }
            else {
            #pragma message (" !!! zb110  not implemented for multitypes structure")
            }


}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::feni()
{
cpos fclp=getLP();
cpos fradius=getRadius(fclp);

cpos ucv=fclp*fclp*fclp; // unit cell volume
cpos maxLp=(cpos) (floor(fradius/fclp));
cpos pSTART=-maxLp*fclp;
position X,Y,Z,Xh,Yh,Zh;
const size_t atomNameA=0;
cpos csize=fradius*2+fclp;
size_t   appSize=(size_t) cube(csize)/ucv;
const size_t mult=4;/// atoms per unit cell for FCC
StAtom atomB;

        atomB.atype=(atomTypes.size()==1) ? atomNameA  : 1;

        ////cube

        atoms.reserve(mult*appSize);

        for(X=pSTART,Xh=pSTART+fclp*0.5;X<=fradius;X+=fclp,Xh+=fclp){
            for(Y=pSTART,Yh=pSTART+fclp*0.5;Y<=fradius;Y+=fclp,Yh+=fclp)
                for(Z=pSTART,Zh=pSTART+fclp*0.5;Z<=fradius;Z+=fclp,Zh+=fclp){
                    atoms.push_back(StAtom(X,Y,Z,atomNameA));
                    atoms.push_back(StAtom(Xh,Yh,Z,atomNameA));
                    atoms.push_back(StAtom(Xh,Y,Zh,atomB.atype));
                    atoms.push_back(StAtom(X,Yh,Zh,atomB.atype));
                }
        }

        atoms.shrink_to_fit();

        
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildFromUC()
{
            if(DB){ cout<<" building by UC replication"<<endl;}

int is,js,ks;
int pm;
const int repX=std::stoi(replicate[0]);
const int repY=std::stoi(replicate[1]);
const int repZ=std::stoi(replicate[2]);

            if(replicate.size()==4)  // +/- option
                { is=-repX;    js=-repY;    ks=-repZ;     pm=2;}
            else
                { is=js=ks=0;    pm=1;     }


const size_t sizeBase=uc.atoms.size();

size_t n;
StVector vi,vj,vk;
position X,Y,Z;
CrandomUni randUniDistr(0,1);
position prob;
//new CrandomUni(a,b) :

            atoms.reserve(pm*repX*pm*repY*pm*repZ*sizeBase);

            //-----------------------------------------

            atomTypes.reserve(uc.atoms.size());

            for(auto &baseAtom: uc.atoms){
            auto iter=std::find(atomTypes.begin(),atomTypes.end(),NanoGrain::StAtomType(baseAtom.name));

                if(iter==atomTypes.end()){
                    atomTypes.push_back(baseAtom.name);                    
                    baseAtom.id=atomTypes.size()-1;
                }
                else
                    baseAtom.id=(size_t) std::distance(atomTypes.begin(),iter);

            }

            atomTypes.shrink_to_fit();
            //-------------------------------------------

            for(int k=ks;k<repZ;k++){
                vk=uc.vtrans[2]*k;

                for(int j=js;j<repY;j++){
                    vj=uc.vtrans[1]*j+vk;

                    for(int i=is;i<repX;i++){
                        vi=uc.vtrans[0]*i+vj;

                        for(n=0;n<sizeBase;n++){
                            X=vi.x+uc.atoms[n].x;
                            Y=vi.y+uc.atoms[n].y;
                            Z=vi.z+uc.atoms[n].z;

                            prob=randUniDistr.randNumber();
                            if(prob>uc.atoms[n].rmProb)
                                atoms.push_back(StAtom(X,Y,Z,uc.atoms[n].id));
                        }
                    }
                }
            }

            for(n=0;n<sizeBase;n++)
                atoms[n].rdh=uc.atoms[n].rdh;

            atoms.shrink_to_fit();

            if(atomTypes.size()>1)
                sortAtomsByName();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildFromTric()
{

StVector &a=tric.a;
StVector &b=tric.b;
StVector &c=tric.c;

cpos deg2rad=M_PI/180;
cpos alpha=std::stod(tric.alpha)*deg2rad;
cpos beta =std::stod(tric.beta)*deg2rad;
cpos gamma=std::stod(tric.gamma)*deg2rad;

cpos A=std::stod(tric.lpa);
cpos B=std::stod(tric.lpb);
cpos C=std::stod(tric.lpc);


                ///
                ///// https://docs.lammps.org/Howto_triclinic.html
                ///

                a.x=A;
                a.y=0;
                a.z=0;

                b.x=B*cos(gamma);
                b.y=B*sin(gamma);
                b.z=0;

                c.x=C*cos(beta);
                c.y=C*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
                c.z=std::sqrt(C*C-sqr(c.x)-sqr(c.y));

                tric.xy=(b.x);
                tric.xz=(c.x);
                tric.yz=(c.y);


                /////////////////////////////////////////////////////


                uc.vtrans.resize(3);
                uc.vtrans[0]=std::move(a);
                uc.vtrans[1]=std::move(b);
                uc.vtrans[2]=std::move(c);

                uc.rdhAtoms=tric.rdhAtoms;

                /////////////////////////////////////////////////////


                /// conversion from triclinic to Cartesian system
                /// see:
                /// https://www.ucl.ac.uk/~rmhajc0/frorth.pdf
                if(tric.csys==StTric::TRIC){
                const size_t nOfatoms=tric.atoms.size();
                vector<StUcAtom> ucAtoms(nOfatoms);
                const auto &tricAtoms=tric.atoms;
                StMatrix mconv;

                        mconv.m11=A; mconv.m12=b.x; mconv.m13=c.x;
                        mconv.m21=0; mconv.m22=b.y; mconv.m23=-C*sin(beta)*cos(beta);
                        mconv.m31=0; mconv.m32=0;   mconv.m33= C*sin(beta)*sin(alpha);

                        for(size_t i=0;i<nOfatoms;i++){
                        const auto &atomt=tricAtoms[i];
                                ucAtoms[i]=atomt;
                                ucAtoms[i]=mconv*StVector(atomt.x,atomt.y,atomt.z);
                        }

                        uc.atoms=std::move(ucAtoms);
                }
                else  /// tric.atoms in Cartesian system
                    uc.atoms=std::move(tric.atoms);

                buildFromUC();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildFromCIF()
{
fstream fin(fileCIF,ios::in);

        if(!fin){
            fin.close();
            errMsg(" file *.CIF is missing");
        throw NanoGrain::Status::ERR_FOPEN;
        }


        ///https://www.iucr.org/__data/assets/pdf_file/0019/22618/cifguide.pdf

std::string fline;
StTric tmpTric;
auto ifbra=[](const char &c){return c=='('|| c==')';};
size_t dataCompletness=0;
vector<string> sgOperators;

        while(!fin.eof() && dataCompletness!=7 ){
            std::getline(fin,fline);
            trim(fline);

            if(fline.empty() || fline[0]=='#' || fline[0]==';') continue;

                                            /// it checks both formats:  1.234567  & 1.2345(67)
            if(regex_match(fline,std::regex("_cell_length_[abc][[:s:]]+[0-9]+[.]?[0-9]*([\\(][0-9]*[\\)])?"))){
            vector<string> tokens(split<string>(fline," \t"));
            string value;

                    value.resize(tokens[1].size());
                    for(size_t i=0,j=0;i<tokens[1].size();i++) {
                        if(ifbra(tokens[1][i]))  continue;

                        value[j++]=tokens[1][i];
                    }

                   switch (tokens[0].back()){
                   case 'a': tmpTric.lpa=value;break;
                   case 'b': tmpTric.lpb=value;break;
                   case 'c': tmpTric.lpc=value;break;
                   }

                   dataCompletness++;

            continue;
            }


            if(regex_match(fline,std::regex("_cell_angle_(alpha|beta|gamma)[[:s:]]+[0-9]+[.]?[0-9]*([\\(][0-9]*[\\)])?"))){
            vector<string> tokens(split<string>(fline," \t"));
            string value;

                    value.resize(tokens[1].size());
                    for(size_t i=0,j=0;i<tokens[1].size();i++) {
                        if(ifbra(tokens[1][i]))  continue;

                        value[j++]=tokens[1][i];
                    }

                   switch (tokens[0].at(12)){
                   case 'a': tmpTric.alpha=tokens[1];break;
                   case 'b': tmpTric.beta=tokens[1];break;
                   case 'g': tmpTric.gamma=tokens[1];break;
                   }

                   dataCompletness++;

            continue;
            }


            if(regex_match(fline,std::regex("loop_"))){

                std::getline(fin,fline);
                trim(fline);

                if(fline.find("_atom_site_")!=string::npos){
                    if(fline.find("_aniso_")!=string::npos) continue;

                int aName,aLabel,aX,aY,aZ;
                int tagPos=0;

                    aName=aX=aY=aZ=-1;

                    do{
                        if(regex_match(fline,std::regex("_atom_site_type_symbol")))
                            aName=tagPos;
                        else
                            if(regex_match(fline,std::regex("_atom_site_label")))
                                aLabel=tagPos;
                            else{
                                if(regex_match(fline,std::regex("_atom_site_fract_[xyz]"))){
                                    switch (fline.back()){
                                    case 'x': aX=tagPos;break;
                                    case 'y': aY=tagPos;break;
                                    case 'z': aZ=tagPos;break;
                                    }
                                }
                        }

                        std::getline(fin,fline);
                        trim(fline);
                        tagPos++;
                    }while(fline[0]=='_' && !fin.eof());

                    if(aName<0 && aLabel>=0){
                        aName=aLabel;
                        warnMsg("_atom_site_type_symbol not given, replaced by _atom_site_label value");
                    }



                    if(! fin.good() || aName<0 || aX<0 || aY<0 || aZ<0 || tagPos<4){
                        errMsg(" the cif file format is corrupted, line: "+std::to_string(__LINE__));
                    throw NanoGrain::Status::ERR_FFORMAT;
                    }


                    do{
                    vector<string> tokens(split<string>(fline," \t"));
                    const size_t tsize=tokens.size();

                        if( (tsize!=static_cast<int>(tagPos))) break;

                    string N(tokens[aName]);
                    string X,Y,Z;//([\\(][0-9]*[\\)])?

                    string value;
                    size_t j=0;
                            value.resize(tokens[aX].size());
                            for(size_t i=0;i<tokens[aX].size();i++) {
                                if(ifbra(tokens[aX][i]))  continue;
                                value[j++]=tokens[aX][i];
                            }
                            X=value.substr(0,j);

                            value.resize(tokens[aY].size());
                            j=0;
                            for(size_t i=0;i<tokens[aY].size();i++) {
                                if(ifbra(tokens[aY][i]))  continue;
                                value[j++]=tokens[aY][i];
                            }
                            Y=value.substr(0,j);
                            j=0;

                            value.resize(tokens[aZ].size());
                            for(size_t i=0;i<tokens[aZ].size();i++) {
                                if(ifbra(tokens[aZ][i]))  continue;
                                value[j++]=tokens[aZ][i];
                            }
                            Z=value.substr(0,j);



                    const string numberFormat("[+-]?([01]([.][0-9]*)?|[.][0-9]+)");
                    auto testX=regex_match(X,std::regex(numberFormat));
                    auto testY=regex_match(Y,std::regex(numberFormat));
                    auto testZ=regex_match(Z,std::regex(numberFormat));

                        if(testX && testY && testZ){

                            tmpTric.atoms.push_back(StUcAtom());
                            tmpTric.atoms.back().name=N;
                            tmpTric.atoms.back().x=std::stod(X);
                            tmpTric.atoms.back().y=std::stod(Y);
                            tmpTric.atoms.back().z=std::stod(Z);

                        }
                        else{
                            errMsg(" couldn't recognize fractonial position : "+fline);
                        throw NanoGrain::Status::ERR_FFORMAT;
                        }

                        std::getline(fin,fline);
                        trim(fline);

                    }while( !fin.eof() &&
                            fline.find("loop")==string::npos && fline.find("data")==string::npos &&
                            fline[0]!='_' && fline[0]!='#' || fline[0]!=';');

                    dataCompletness++;
                }

                else
                    if(fline.find("_space_group_symop_operation_xyz")!=string::npos ||
                            fline.find("_symmetry_equiv_pos_as_xyz")!=string::npos){
                    const string cmdName(fline);
                    std::streampos streamPos;

                        sgOperators.reserve(200);

                        while(!fin.eof()){
                            streamPos=fin.tellg();
                            std::getline(fin,fline);
                            trim(fline);

                            if(fline[0]=='_' || fline[0]==';' || fline[0]=='#' ||
                               fline.empty() ||
                                fline.find("loop")!=string::npos){  break; }

                            try{
                            auto itr=std::find_if(fline.begin(),fline.end(),[](const char &ch) {return ch==',';});
                                if(itr!=fline.end()){
                                        itr=std::find_if(itr,fline.end(),[](const char &ch) {return ch==',';});
                                        if(itr!=fline.end())
                                            sgOperators.push_back(fline);
                                        else
                                            throw 1;
                                }
                                else throw 0;
                            }
                            catch(int i)
                            {
                                errMsg("unrecognized/wrong format of " +cmdName+" section within CIF file");
                                fin.close();
                            throw NanoGrain::Status::ERR_FFORMAT;
                            }
                        }
                        fin.seekg(streamPos);
                    }/////

                }
            }/// while

        fin.close();

        if(dataCompletness!=7){
            errMsg(" dataCompletness != 7");
        throw NanoGrain::Status::ERR_FFORMAT;
        }

        if(!tmpTric.atoms.empty() && !sgOperators.empty()){
        vector<StAtomFracPostion> acceptedAtoms;

                    for(auto &ap:  tmpTric.atoms){                        
                        for(auto &sgFormula: sgOperators){
                        StSymmGroupGenerator stg(buildGenerator(sgFormula));
                        StAtomFracPostion    ucAtomPos(stg.ucAtomPos(StAtomFracPostion(ap.x,ap.y,ap.z)));
                              ucAtomPos.reduct();

                        auto itr=std::find_if(acceptedAtoms.begin(),acceptedAtoms.end(),ucAtomPos);
                            if(itr==acceptedAtoms.end()){
                                acceptedAtoms.emplace_back(ucAtomPos);
                                acceptedAtoms.back().name=ap.name;
                            }



                        }
                    }

                    tric=std::move(tmpTric);
                    tric.atoms.clear();
                    tric.atoms.reserve(acceptedAtoms.size());

                    for(auto &a:acceptedAtoms){
                        tric.atoms.emplace_back(StUcAtom(a.x,a.y,a.z,a.name));
                    }


        }
        else
            tric=std::move(tmpTric);


        buildFromTric();
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildHcpABC()
{
size_t laynum;

            if(shape=="cylinder") {
                    if(hcplaynum.empty())
                    throw NanoGrain::ERR_UNKSEQNUM;
               laynum=std::stoi(hcplaynum);
            }
            else
                if(shape=="sphere")
                    laynum=(hcplaynum.empty()) ? (size_t) std::ceil(2*hcpRadius/hcpc) : (size_t) std::stoi(hcplaynum);
                else
                throw NanoGrain::ERR_WRSHAPE;

size_t i,j;
const size_t seqlen=hcpfillup.length();

            hcpABC.resize(laynum);

            for(i=0,j=0;i<laynum;i++,j++){
                j%=seqlen;
                hcpABC[i]=hcpfillup[j];
            }

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildRandomABC(vector<string> &params)
{
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);

float a,b;
unsigned seqLength;

        a=std::stof(params[1]);
        b=std::stof(params[2]);


        if(params[0]=="uniform"){
        std::uniform_int_distribution<unsigned> distribution((int)a, (int) b);
            seqLength=distribution(generator);
        }
        else
            if(params[0]=="normal"){
            std::normal_distribution<float>    distribution(a,b);
                seqLength=(unsigned) std::round(distribution(generator));
            }
            else{
            std::lognormal_distribution<float>    distribution(a,b);
                seqLength=(unsigned) std::round(distribution(generator));

            }


        if(DB) cout<<"seqLength: "<<seqLength<<endl;


        hcpABC.clear();

         if(params[3]=="abc"){
             for(size_t i=0;i<seqLength;i++)
                 hcpABC+="ABC";
         }
         else
             if(params[3]=="random"){
             std::uniform_int_distribution<int> distributionIni(0,2);
             int nr=distributionIni(generator);

                    switch(nr){
                    case 0 : hcpABC="a"; break;
                    case 1 : hcpABC="b";break;
                    case 2 : hcpABC="c";break;
                    }

             std::uniform_int_distribution<int> distribution(1,2);

                    for(size_t i=1;i<seqLength;i++){
                        nr+=distribution(generator);
                        nr%=3;
                        switch(nr){
                        case 0 : hcpABC+="a"; break;
                        case 1 : hcpABC+="b";break;
                        case 2 : hcpABC+="c";break;
                        }
                    }
             }

         if(DB) {cout<<"hcpABC: "<<hcpABC<<endl;}

}
//-----------------------------------------------------------------------------
void pushAtomCScircle(NanoGrain::vatoms &atoms,const NanoGrain::StAtom & atom,cpos & radii2, NanoGrain::StNanoGrain *ptr)
{
        (void) ptr;
        if(atom.r2>radii2)
            return;

        atoms.push_back(atom);

}
//-----------------------------------------------------------------------------
void pushAtomCSHex(NanoGrain::vatoms &atoms,const NanoGrain::StAtom & atom,cpos & radii2, NanoGrain::StNanoGrain *ptr)
{
        (void)radii2;
        (void)ptr;
        atoms.push_back(atom);
}
//-----------------------------------------------------------------------------
void pushAtomCSPoly(NanoGrain::vatoms &atoms,const NanoGrain::StAtom & atom,cpos & radii2, NanoGrain::StNanoGrain *ptr)
{

    if(ptr->ssShape->isValid(atom.x,atom.y,0))
        atoms.push_back(atom);
}


//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildHcpBaseHex()
{
            a.x=0;         a.y= hcpa;     a.z=0;
            b.x= hcpa*wj;  b.y=-hcpa*wi;  b.z=0;
            c=a+b;

            baseAtomsHcp.clear();

            //@INFO
            // maxAtoms= 0.5 * (area of base /area of unit cell)
            // area of unit cell is equal sides (hcpa) triangle
            // area of base is equal area of hexagon with side equal radius
            // 0.5 - each unit cell has 0.5 atom


const size_t maxAtoms=  0.5* 6 * sqr(hcpRadius/hcpa) * atomTypes.size();

            baseAtomsHcp.reserve(maxAtoms);


void (*pushAtom)(NanoGrain::vatoms &,const NanoGrain::StAtom & ,cpos &, NanoGrain::StNanoGrain *);

           // pushAtom=(hcpcs=="circle") ? &pushAtomCScircle :
                 //    (hcpcs=="hex"   ) ? &pushAtomCSHex    : &pushAtomCSPoly;


            if(hcpcs=="circle") { pushAtom=&pushAtomCScircle;}
            else
                if(hcpcs=="hex") {pushAtom=&pushAtomCSHex;   }
                else{
                vector<string> tokens(split<string>(shapePrm2D," "));
                    if(tokens.size()==3){
                    const double xa=std::stod(tokens[0]); //elongation paramter
                    const double yb=std::stod(tokens[1]); //elongation paramter
                    const double p=std::stod(tokens[2]);   //polyhedrality parameter

                        ssShape=new CPolyhedral2D(p,hcpRadius*0.75);
                        ssShape->xa=xa;
                        ssShape->yb=yb;

                        pushAtom=&pushAtomCSPoly;
                    }
                    else { //tokens.size==5
                    cdouble p=std::stod(tokens[0]);  // polyhedrality parameter
                    cdouble a=std::stod(tokens[1]);  // select shape
                    cdouble b=std::stod(tokens[2]);  // select shape
                    cdouble xa=std::stod(tokens[3]); // elongation parameter
                    cdouble yb=std::stod(tokens[4]); // elongation paramter

                            ssShape=new CPolyhedral2D_HOD(p,hcpRadius*0.75,a,b);

                            ssShape->xa=xa;
                            ssShape->yb=yb;

                            pushAtom=&pushAtomCSPoly;
                    }
                }


const int maxN= (int) std::ceil(hcpRadius/hcpa)+1;
int l;
StVector va,vb,vd;
StAtom atom;

             // 1. step  - string of atoms along diameter
             for(l=-maxN;l<=maxN;l++){
                 va=a*l;

                 atom.x=va.x;
                 atom.y=va.y;
                 atom.z=0;
                 atom.calcR2();

                 pushAtom(baseAtomsHcp,atom,hcpRadius2,this);
             }

             vd=a*maxN;

             // 2. step - fill up the cross section

             for(l=1;l<=maxN;l++){
                 va=vd*(-1.0)+c*l;
                 vb=vd-c*l;

                 for(int p=2*maxN-l+1;p>0;p--){

                    atom.x=va.x;  // atom on the left side of cross section
                    atom.y=va.y;
                    atom.z=0;
                    atom.calcR2();

                    pushAtom(baseAtomsHcp,atom,hcpRadius2,this);

                    atom.x=vb.x;  // atom on the right side of cross section
                    atom.y=vb.y;
                    atom.z=0;
                    atom.calcR2();

                    pushAtom(baseAtomsHcp,atom,hcpRadius2,this);

                    va=va+a;
                    vb=vb-a;
                  }
             }


            baseAtomsHcp.shrink_to_fit();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildHcpBaseRhm()
{

            //a.x=0;         a.y=-hcpa;     a.z=0;
            //b.x= hcpa*wj;  b.y=hcpa*wi;  b.z=0;
    
            b.y=0;         b.x= -hcpa;     b.z=0;
            a.y= hcpa*wj;  a.x= hcpa*wi;  a.z=0;                                                
            c=a+b;

            baseAtomsHcp.clear();

            //@INFO
            // maxAtoms= 0.5 * (area of base /area of unit cell)
            // area of unit cell is equal sides (hcpa) triangle
            // area of base is equal area of hexagon with side equal radius
            // 0.5 - each unit cell has 0.5 atom

const int maxN= (int) std::ceil(hcpRadius/hcpa)+1;
StVector va,vb;
StAtom atom;
//const size_t maxAtoms=  0.5* 6 * sqr(hcpRadius/hcpa) * atomTypes.size();

            baseAtomsHcp.reserve(maxN*maxN*4);


            for(int i=-maxN;i<maxN;i++){
                va=a*i;

                for(int j=-maxN;j<maxN;j++){
                    vb=va+ b*j;

                    atom.x=vb.x;
                    atom.y=vb.y;
                    atom.z=0;
                    atom.calcR2();

                    baseAtomsHcp.push_back(atom);

                }
            }

            baseAtomsHcp.shrink_to_fit();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::insertFaults()
{
const int hcpFaultNum=std::stoi(hcpfault);


            if(hcpFaultNum==0){
                cerr<<"WARNING: number of faults is equal 0; operation ignored"<<endl;
            return;
            }

            if(hcpFaultNum<0){
                cerr<<"WARNING: number of faults is less than 0; operation ignored"<<endl;
            return;
            }



            if(DB) cout<<" ABC original: \n"<<hcpABC<<" "<<hcpABC.length()<<endl;

            //decrease number of ini ABC layers by number of fault layers
            hcpABC=hcpABC.substr(0,hcpABC.length()-hcpFaultNum);

const size_t laynum=hcpABC.size();
const size_t frac=laynum/hcpFaultNum;


            if(frac<3){
                cerr<<"WARNING: number of faults is too big"<<endl;
            return;
            }


            if(DB) cout<<" ABC reducted: \n"<<hcpABC<<" "<<hcpABC.length()<<endl;


            hcpFaultPos.reserve(hcpFaultNum);

            if(faultmode==StNanoGrain::RANDOM){
            std::default_random_engine generator (std::chrono::system_clock::now().time_since_epoch().count());
            size_t from,to;
            int i;


                    //// UWAGA: algorytm wstawia sf tylko WEWNATRZ wybranej sekwencji, nigdy po bokach,
                    /// dlatego distribution jest w zakresie (from,to-1);
                    ///
                    /// WARNING: algorithm inserts sf inside of selected sequence, never adds sf on boundaries
                    /// so range of distribution is (from,to-1);




                    for(i=0,from=0,to=frac-1;i<hcpFaultNum;i++){
                    std::uniform_int_distribution<unsigned> distribution(from,to-1);
                    const unsigned randPos=distribution(generator);

                            if(DB){
                                cout<<"from, to, randPos: "<<from<<", "<<to<<", "<<randPos<<endl;

                                if(randPos>hcpABC.length()){
                                    cerr<<"ERROR: randPos out of range"<<endl;
                                continue;
                                }
                            }

                            hcpFaultPos.push_back(randPos+1);
                    string subABC(hcpABC.substr(randPos,2));

                            subABC[0]=std::tolower(subABC[0]);
                            subABC[1]=std::tolower(subABC[1]);


                            try{

                                switch (subABC[0]){
                                case 'a':
                                    hcpABC.insert(randPos+1, (subABC[1]=='b') ? "c": "b"); break;

                                case 'b':
                                    hcpABC.insert(randPos+1, (subABC[1]=='c') ? "a": "c"); break;

                                case 'c':
                                    hcpABC.insert(randPos+1, (subABC[1]=='a') ? "b": "a"); break;
                                }

                            }
                            catch(const std::out_of_range& e){
                                cerr<<"WARNING: the insert fault exception has been thrown"<<endl;
                            }

                            if(DB){cout<<"hcpABC.length(): "<<hcpABC.length()<<endl;}

                            from+=frac;
                            to=from+frac;
                            from+=1;

                    }

            }

            if(DB) {
                cout<<" ABC after update: "<<endl;
                //cout<<hcpABC<<" "<<hcpABC.length()<<endl;

                for(size_t i=0;i<hcpABC.length();i++){

                    if(i%(1+frac)==0) cout<<" ";

                    cout<<hcpABC[i];
                }

                cout<<endl;

                for(size_t i=0;i<hcpABC.length()-1;i++){
                    cout<<" ";
                    if(std::tolower(hcpABC[i])==std::tolower(hcpABC[i+1]))
                        cout<<"^"<<i;
                }
            }


}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildHcpLayers()
{
cpos ht=hcpa*sqrt(3)/2; //triangule height
StVector b,c,vd,vres;

             b.x=ht/3; // shift 1/3
             b.y=hcpa/2.0;
             b.z=0;

             c.x=ht*2.0/3.0; //shift 2/3
             c.y=0;
             c.z=0;

             d.x=0;	d.y=0;	d.z=hcpc;


const size_t seqSize=hcpABC.length();
size_t i,layer0Atoms=0,layerPrevAtoms=0,layerLastAtoms=0;
int ii;
const size_t atomNameA=0;


            atoms.reserve(seqSize*baseAtomsHcp.size()*atomTypes.size());

            for(i=0,ii=-seqSize/2;i<seqSize;i++,ii++){

                vd=d*ii;

                switch(hcpABC[i]){
                case 'a':
                case 'A': vres=  vd;break;
                case 'b':
                case 'B': vres=b+vd;break;
                case 'c':
                case 'C': vres=c+vd;break;
                default:  cerr<<"ERROR: unknown sequence type: "<<hcpABC[i]<<endl; return;
                }

                for(StAtom &atom : baseAtomsHcp)
                    (this->*fbuildHcpShape)(atom,vres);

                if(!layer0Atoms)
                    layer0Atoms=atoms.size();


                if(atoms.size()>layerPrevAtoms){
                    layerLastAtoms=atoms.size()-layerPrevAtoms;
                    layerPrevAtoms=atoms.size();
                }

             }

            if(atomTypes.size()==2 || hcpsl ){
            StAtom atomB;

            size_t i;
            cpos prc=std::stod(hcpu);
            cpos dc=d.z*prc;

                    //atomB.p_name=(atomNameB.empty()) ? atomNameA  : &atomNameB;

                        atomB.atype=(atomTypes.size()==1) ? atomNameA  : 1;


                    /// dwa typy powierz
                    if(hcpsurf==EHCPSURF::sB || hcpsurf==EHCPSURF::sH){
                    const size_t asize=atoms.size();

                            for(i=0;i<asize;i++){
                                atomB.x=atoms[i].x;
                                atomB.y=atoms[i].y;
                                atomB.z=atoms[i].z+dc;                                
                                atomB.calcR2();

                                atoms.push_back(atomB);
                            }


                            if(hcpsurf==EHCPSURF::sH){

                                atomTypes.emplace_back(StAtomType("H"));
                            const size_t atypeHpos=atomTypes.size()-1;


                                for(i=0;i<layer0Atoms;i++){
                                StAtom atomH(atoms[i]);

                                    atomH.z-=1.1;
                                    atomH.atype=atypeHpos;
                                    atomH.calcR2();
                                    atoms.emplace_back(atomH);
                                }


                                for(i=2*asize-layerLastAtoms;i<2*asize;i++){
                                StAtom atomH(atoms[i]);

                                    atomH.z+=1.1;
                                    atomH.atype=atypeHpos;
                                    atomH.calcR2();
                                    atoms.emplace_back(atomH);

                                }
                            }
                    }
                    else   {                           
                    const size_t asize=(hcpsurf==EHCPSURF::sA) ? atoms.size()-layerLastAtoms
                                                               : atoms.size(); // surf AB

                    //const size_t asize=atoms.size()-atomsExcLastLayer;

                            for(i=0;i<asize;i++){
                                atomB.x=atoms[i].x;
                                atomB.y=atoms[i].y;
                                atomB.z=atoms[i].z+dc;

                                atoms.push_back(atomB);
                            }

                            atoms.erase(atoms.begin(),atoms.begin()+layer0Atoms);
                    }
            }


            atoms.shrink_to_fit();
            baseAtomsHcp.clear(); // clear memory

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildHcpCylinder(const StAtom &atom, const StVector &vres)
{
StAtom tmpatom;

                tmpatom.x=atom.x+vres.x;
                tmpatom.y=atom.y+vres.y;
                tmpatom.z=atom.z+vres.z;
                tmpatom.calcR2();
                tmpatom.atype=0;
                atoms.emplace_back(tmpatom);
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildHcpSphere(const NanoGrain::StAtom &atom, const StVector &vres)
{
StAtom tmpatom;

                tmpatom.x=atom.x+vres.x;
                tmpatom.y=atom.y+vres.y;
                tmpatom.z=atom.z+vres.z;
                tmpatom.calcR2();
                tmpatom.atype=0;

                if(tmpatom.r2<=hcpRadius2)
                    atoms.emplace_back(tmpatom);
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildHcpSuperSphere(const NanoGrain::StAtom &atom, const StVector &vres)
{
StAtom tmpatom;

                    tmpatom.x=atom.x+vres.x;
                    tmpatom.y=atom.y+vres.y;
                    tmpatom.z=atom.z+vres.z;

                    if(ssShape->isValid(tmpatom.x,tmpatom.y,tmpatom.z)){
                        tmpatom.calcR2();
                        tmpatom.atype=0;
                        atoms.emplace_back(tmpatom);

                    }
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::buildSuperSphere()
{
StMinMax minmax;

            grainCentering();
            minmax.searchForMinMax(atoms);

vector<string> shapeParams(split<string>(shapePrm," "));
const double p=shapePrmValue(shapeParams[0]);
const double R=minmax.getMinPos()*1.05;

            if(shape=="cubic")
                if(shapeParams.size()==1) ssShape=new CCubic(p,R);
                else ssShape=new CCubic(p,R,std::stod(shapeParams[1]),std::stod(shapeParams[2]),std::stod(shapeParams[3]));
            else
                if(shape=="oct")
                    ssShape=new COctahedral(p,R);
                else
                    if(shape=="dod")
                        ssShape=new CDodecahedral(p,R);
                    else {//poly
                    const double a=shapePrmValue(shapeParams[1]);
                    const double b=shapePrmValue(shapeParams[2]);


                        if(shapeParams.size()==1) ssShape=new CPolyhedral(p,R,a,b);
                        else ssShape=new CPolyhedral(p,R,a,b,std::stod(shapeParams[3]),std::stod(shapeParams[4]),std::stod(shapeParams[5]));
                    }


vatoms tmpAtoms(std::move(atoms));

                atoms.clear();
                atoms.reserve(tmpAtoms.size());

                for(auto & atom : tmpAtoms){
                    if(ssShape->isValid(atom.x,atom.y,atom.z)){
                        atoms.emplace_back(atom);
                    }
                }

                atoms.shrink_to_fit();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::grainCentering()
{

        omp_set_num_threads(stoi(threads));

position Xm,Ym,Zm;
const size_t numOfAtoms=atoms.size();

        Xm=Ym=Zm=0;
//if(numOfAtoms>10000)
        #pragma  omp parallel for reduction(+:Xm,Ym,Zm)
        for(size_t i=0; i<numOfAtoms; i++){
            Xm+=atoms[i].x;
            Ym+=atoms[i].y;
            Zm+=atoms[i].z;
        }


const position tX=Xm/numOfAtoms;
const position tY=Ym/numOfAtoms;
const position tZ=Zm/numOfAtoms;

        #pragma omp parallel for
        for(size_t i=0;i<numOfAtoms;i++){
            atoms[i].x-=tX;
            atoms[i].y-=tY;
            atoms[i].z-=tZ;
        }


}
//-----------------------------------------------------------------------------
bool noATypeCompare(const NanoGrain::StAtom &a,const position &r2, int &atype)
{
    return a.r2<r2;
}

//-----------------------------------------------------------------------------
bool whATypeCompare(const NanoGrain::StAtom &a,const position &r2, int &atype )
{
    return (a.atype==static_cast<size_t>(atype)) && (a.r2<r2);
    //if (a.atype==atype && a.r2<r2)
     //   return true;

///return false;
}

//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::catomCentering()
{
position Xm,Ym,Zm;
const size_t numOfAtoms=atoms.size();

        /// searching for geometrical center
        grainCentering();

        omp_set_num_threads(stoi(threads));


size_t id;
position r2=INFINITY;

int atype;
bool (*ptr_cmp)(const NanoGrain::StAtom &a,const position &r2, int &atype);

        if(catomType.empty())
            ptr_cmp=&noATypeCompare;
        else{
            if( (atype=findAtomName(catomType)) <0){
               cerr<<"WARNING: atom type "<<catomType<<" unknown; centering by catom disabled"<<endl;
            return;
            }
            else
               ptr_cmp=&whATypeCompare;
        }


        #pragma omp parallel
        {
        position r2_local=INFINITY;
        size_t   id_local;


            #pragma omp for
            for (size_t i=0;i<numOfAtoms;i++){
                atoms[i].calcR2();

                //if(atoms[i].r2<r2_local){
                if(ptr_cmp(atoms[i],r2_local,atype)){
                    r2_local=atoms[i].r2;
                    id_local=i;
                }
            }

            #pragma omp critical
            {
                if(r2_local<r2){
                    r2=r2_local;
                    id=id_local;
                }
            }

            //#pragma omp barier


            #pragma omp single
            {
                Xm=atoms[id].x;
                Ym=atoms[id].y;
                Zm=atoms[id].z;
            }

            //#pragma omp barier
            #pragma omp for
            for (size_t i=0;i<numOfAtoms;i++){
                atoms[i].x-=Xm;
                atoms[i].y-=Ym;
                atoms[i].z-=Zm;
            }
        }


}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::idCentering()
{
const size_t aid=std::stoi(catomType);
const size_t numOfAtoms=atoms.size();


            if(aid>=numOfAtoms){
                cerr<<"WARNING: atom ID "<<aid<<" is out of range; centering by atom's ID disabled"<<endl;
             return;
            }

position Xm,Ym,Zm;

            Xm=atoms[aid].x;
            Ym=atoms[aid].y;
            Zm=atoms[aid].z;


            #pragma omp parallel for
            for (size_t i=0;i<numOfAtoms;i++){
                atoms[i].x-=Xm;
                atoms[i].y-=Ym;
                atoms[i].z-=Zm;
            }


}

//---------------------------------------------------------------------------



//---------------------------------------------------------------------------

void NanoGrain::StNanoGrain::insertDisloc()
{
vector<string> toks{split<string>(dislocPlane," ")};

//plane axis params
cpos  A{std::stod(toks[0])};
cpos  B{std::stod(toks[1])};
cpos  C{std::stod(toks[2])};
cpos  D{std::stod(toks[3])};
cpos  rmin{std::stod(toks[4])};
cpos  rmax{std::stod(toks[5])};
const int N{std::stoi(toks[6])};

// rotation axis =   [0,0,1] x [A,B,C] = [-B,A,0]
StAxis axis(-B,A,0);
cpos mian=std::sqrt(A*A+B*B+C*C);
cpos sa=std::sqrt(A*A+B*B)/mian;
cpos ca=C/mian;
StRotationMatrix rotMat(axis,sa,ca);

vatoms planeAtoms;

std::default_random_engine generatorRadius (std::chrono::system_clock::now().time_since_epoch().count());
std::uniform_real_distribution<double> rdistr(rmin*rmin,rmax*rmax);

std::default_random_engine generatorAngle (std::chrono::system_clock::now().time_since_epoch().count());
std::uniform_real_distribution<double> adistr(0,2*M_PI);

double sqR,Ang;
std::string aname("C");
size_t atype;
StVector v;
double dTh=2*M_PI/(N-1);

            atomTypes.push_back(StAtomType(aname));
            atype=atomTypes.size()-1;
            planeAtoms.reserve(N);
            Ang=0;

            for(size_t i=0;i<N;i++){
                sqR=std::sqrt(rdistr(generatorRadius));
                Ang=dTh*i;
                //Ang=adistr(generatorAngle);

                v.x=sqR*std::cos(Ang);
                v.y=sqR*std::sin(Ang);
                v.z=D;

                v=rotMat*v;

                planeAtoms.emplace_back(StAtom(v.x,v.y,v.z,atype));
            }


            atoms.insert(atoms.begin(),planeAtoms.begin(),planeAtoms.end());

}
//-----------------------------------------------------------------------------
bool NanoGrain::StNanoGrain::testSavedNumOfAtoms(const size_t numOfAtoms)
{

auto & sna=NanoGrain::StNanoGrain::savedNumOfAtoms;

auto isEqual=[&](const size_t v){return numOfAtoms==v;};
auto it=std::find_if(std::begin(sna),std::end(sna),isEqual);

return it!=std::end(sna);

}
//-----------------------------------------------------------------------------

void NanoGrain::StNanoGrain::build()
{

    if(structure.find("uc")!=str::npos)
        buildFromUC();
    else
        if(structure.find("tric")!=str::npos)
            buildFromTric();
        else{
            if(structure.find("cif")!=str::npos)
                buildFromCIF();
            else
                if(structure.find("sc")!=str::npos)
                            sc();
                        else
                            if(structure.find("bcc")!=str::npos)
                                bcc();
                            else
                                if(structure.find("fcc")!=str::npos)
                                    fcc();
                                else
                                    if(structure.find("zb110")!=str::npos)
                                        zb110();
                                    else
                                        if(structure.find("uo2")!=str::npos)
                                            fcc(EFCCTYPE::UO2);
                                        else
                                            if(structure.find("feni")!=str::npos)
                                                feni();
                                            else
                                                if(structure.find("zb")!=str::npos)
                                                    fcc(EFCCTYPE::ZB);
                                                else{
                                                    hcp();
                                                return;
                                                }
        }


        if(!shapePrm.empty())
            if(shape!="cone")
                buildSuperSphere();

}

//-----------------------------------------------------------------------------

void NanoGrain::StNanoGrain::voids()
{
vatoms atomsTmp;
std::default_random_engine generator (std::chrono::system_clock::now().time_since_epoch().count());
std::uniform_real_distribution<float> distribution(0,1);
vector<string> tokens(split<string>(voidsprm," "));//scale factors
float prob;

                if(tokens.size()==1)
                    prob=std::stof(tokens[0]);
                else{
                std::uniform_real_distribution<float> distributionProb( stof(tokens[0]), stof(tokens[1]));
                    prob=distributionProb(generator);
                }

                atomsTmp.reserve(atoms.size());

                for(auto & atom:atoms){
                    if(prob<=distribution(generator))
                        atomsTmp.emplace_back(atom);
                }

                atomsTmp.shrink_to_fit();

                if(dispParams &&  !atoms.empty() ){
                const size_t arem=atoms.size()-atomsTmp.size();
                    atomsRemoved.str("");
                    atomsRemoved<<arem<<" "<<setw(5)<<setprecision(3)<<100.0*arem/atoms.size()<<"%";
                }

                atoms=std::move(atomsTmp);
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::disperseN()
{
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
position rx,ry,rz;
vector<position> vdf;
size_t i;
                vdf.reserve(atomTypes.size());

                if(atomDisperse[0].name=="*"){
                    for(size_t i=0;i<atomTypes.size();i++)
                        vdf.push_back(std::stod(atomDisperse[0].dispfact));
                }
                else{
                    for(auto &atype:atomTypes){
                    const std::string anames(atype.name);
                            for(i=0;i<atomDisperse.size();i++)
                                if(atomDisperse[i].name==anames){
                                    vdf.push_back(std::stod(atomDisperse[i].dispfact));
                                    break;
                                }
                    }


                    if(vdf.size()!=atomTypes.size()){
                        cerr<<"ERROR: unrecognized atom name(s)"<<endl;
                    throw NanoGrain::ERR_UNKAN;
                    }
                }


                for(StAtom &atom:atoms){
                std::normal_distribution<position> distribution(0,vdf[atom.atype]);

                            rx=distribution(generator);
                            ry=distribution(generator);
                            rz=distribution(generator);

                            atom.x+=rx;
                            atom.y+=ry;
                            atom.z+=rz;
                            atom.calcR2();
                }


}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::displaceAtoms()
{


position currRadius,newRadius,corr;

            for(StAtom & atom : atoms){

                if(atom.r2<1e-9) continue;

                currRadius=std::sqrt(atom.r2);
                newRadius=coreshell.newR(currRadius);
                corr=newRadius/currRadius;

                atom*=corr;

            }

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::rescale()
{
vector<string> sfxyz(split<string>(scaleFactors," "));//scale factors

            if(sfxyz.size()==1){
            position sfx;

                sfx=std::stod(sfxyz[0]);

                for(StAtom &atom: atoms)
                    atom.rescale(sfx);
            }
            else{
            position sfx,sfy,sfz;

                sfz=std::stod(sfxyz[2]);
                sfy=std::stod(sfxyz[1]);
                sfx=std::stod(sfxyz[0]);

                for(StAtom &atom: atoms)
                    atom.rescale(sfx,sfy,sfz);
            }
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::renameAtoms()
{    
            for(auto & fromToRename: rename){
                renameAtomsFromTo(fromToRename);
            }

            sortAtomsByName();
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::renameAtomsFromTo(const string fromToProb)
{
vector<string> toks(split<string>(fromToProb," "));
const size_t toksNumber=toks.size();
constexpr size_t argSize=3;
const size_t pairsNumber=toksNumber/argSize;

struct StTypeFromTo{
                    string fname,tname;
                    size_t from,to;
                    double prob;

                    bool operator == (const size_t &k) { return k==from; }
                    };

vector<StTypeFromTo> vfromTo(pairsNumber);

                    for(size_t i=0,k;i<pairsNumber;i++){
                        k=i*argSize;
                        vfromTo[i].fname=toks[k];
                        vfromTo[i].tname=toks[k+1];
                        vfromTo[i].prob=std::stod(toks[k+2]);
                    }



vector<StAtomType> tmpAtomTypes;

                    tmpAtomTypes.resize(atomTypes.size());
                    std::copy(atomTypes.begin(),atomTypes.end(),tmpAtomTypes.begin());

                    for(auto & at: vfromTo){
                    auto itr=std::find(tmpAtomTypes.begin(),tmpAtomTypes.end(),at.fname);

                        if(itr!=tmpAtomTypes.end())
                            at.from=std::distance(tmpAtomTypes.begin(),itr);
                        else
                            continue;


                        itr=std::find(tmpAtomTypes.begin(),tmpAtomTypes.end(),at.tname);

                        if(itr==tmpAtomTypes.end())
                            atomTypes.emplace_back(at.tname);

                        at.to=std::distance(tmpAtomTypes.begin(),itr);

                    }


unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
std::uniform_real_distribution<double> distribution(0,1);
double rvalue;

                    for(auto &atom: atoms){
                    auto itr=std::find(vfromTo.begin(),vfromTo.end(),atom.atype);

                        if(itr==vfromTo.end())
                            continue;

                        rvalue=distribution(generator);

                        if(rvalue<itr->prob)
                            atom.atype=itr->to;

                    }

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::removeAtoms()
{
        if(!rmatoms.empty()) removeRandomAtoms();

        for(auto &rPrms: vremoveAtomsPrms)
            removeAtomsPlane(rPrms);

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::removeRandomAtoms()
{
size_t i,j;
const size_t numOfAtoms=atoms.size();
vector<string> tokens(split<string>(rmatoms," \t"));
const unsigned int minNumOfBonds=std::stoi (tokens[0]);
const double bondLenght=std::stod(tokens[1]);
const double bL2=bondLenght*bondLenght;
auto sqr=[](const double &x, const double &y, const double &z)  { return x*x+y*y+z*z;};

        //cout<<" threads "<<threads<<endl;

        omp_set_num_threads(std::stoi(threads));


        #pragma omp parallel for if(numOfAtoms>10000)
        for(i=0;i<numOfAtoms;i++)
            atoms[i].neighID.reserve(4);


        for(i=0;i<numOfAtoms-1;i++){

             #pragma omp parallel if(numOfAtoms>10000)
             {
             double dx,dy,dz;
             StAtom * const ptr_fatom=&atoms[i];

                    #pragma omp for
                    for(j=i+1;j<numOfAtoms;j++){
                        dx=atoms[j].x-ptr_fatom->x;
                        dy=atoms[j].y-ptr_fatom->y;
                        dz=atoms[j].z-ptr_fatom->z;

                        if(sqr(dx,dy,dz)<=bL2){
                            atoms[j].neighID.push_back(i);

                            #pragma omp critical
                            {
                                ptr_fatom->neighID.push_back(j);
                            }

                        }
                     }
             }
        }

vatoms tmpAtoms(std::move(atoms));

            atoms.clear();
            atoms.reserve(numOfAtoms);

            if(tokens.size()==3){
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator (seed);
            std::uniform_real_distribution<double> distribution(0,100);
            const double prob=std::stod(tokens[2]);


                    for(i=0;i<numOfAtoms;i++){
                        if(tmpAtoms[i].neighID.size()>minNumOfBonds)
                                atoms.emplace_back(tmpAtoms[i]);
                        else {
                            if(distribution(generator)>prob)
                                atoms.emplace_back(tmpAtoms[i]);
                        }
                    }

            }
            else {
                    for(i=0;i<numOfAtoms;i++){
                        if(tmpAtoms[i].neighID.size()>minNumOfBonds)
                            atoms.emplace_back(tmpAtoms[i]);
                    }
            }

            atoms.shrink_to_fit();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::removeAtomsPlane(const string &prms)
{
vector<string> toks{split<string>(prms," ")};
const bool out=(toks[1]=="out");
const position A=std::stod(toks[2]);
const position B=std::stod(toks[3]);
const position C=std::stod(toks[4]);
const position D=std::stod(toks[5]);
const size_t numOfAtoms=atoms.size();
position v;
vatoms tmpAtoms(std::move(atoms));

            atoms.clear();
            atoms.reserve(numOfAtoms);

            for (auto &atom: tmpAtoms){
                v=A*atom.x+B*atom.y+C*atom.z+D;

                if( (v<=0) ^ !out)
                    atoms.push_back(atom);
            }

            atoms.shrink_to_fit();
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::rotate()
{
vector<string> toks{split<string>(rotatePrm," ")};
const size_t nOftoks{toks.size()};

cpos deg2rad=M_PI/180;
cpos alpha{std::stod(toks[0])*deg2rad};
cpos sa(std::sin(alpha));
cpos ca(std::cos(alpha));

cpos A(std::stod(toks[1]));
cpos B(std::stod(toks[2]));
cpos C(std::stod(toks[3]));

cpos  px{  (nOftoks>4) ? std::stod(toks[4]) : 0};
cpos  py{  (nOftoks>4) ? std::stod(toks[5]) : 0};
cpos  pz{  (nOftoks>4) ? std::stod(toks[6]) : 0};

StVector axisPos(px,py,pz);
StAxis rotAxis(A,B,C);
StRotationMatrix rotMat(rotAxis,sa,ca);
StVector v_xyz,vrot;


           for (auto & atom: atoms){
               v_xyz=atom.Pos()-axisPos;
               vrot=rotMat*v_xyz;

               atom.x=vrot.x;
               atom.y=vrot.y;
               atom.z=vrot.z;

           }

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::findMaxR()
{

            maxR=0;
            maxX=0;
            maxY=0;
            maxZ=0;

const size_t size=atoms.size();

            omp_set_num_threads(std::stoi(threads));

            #pragma omp parallel for reduction(max: maxX,maxY,maxZ,maxR)  if(size>99)
            for(size_t i=0;i<size;i++){
            auto &atom=atoms[i];

                if(atom.x>maxX)  maxX=atom.x;
                if(atom.y>maxY)  maxY=atom.y;
                if(atom.z>maxZ)  maxZ=atom.z;
                if(atom.r2>maxR) maxR=atom.r2;
            }

            maxR=std::sqrt(maxR);
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::findAtomNamesNumber()
{
size_t counter;
size_t atype,j,k;
size_t atomsSize=atoms.size();

        atomNamesNumber.resize(atomTypes.size());

        for(atype=0;atype<atomTypes.size();atype++){

            /// find first atom in the buffer with a given type
            for(j=0;j<atomsSize;j++)
                if(atoms[j].atype==atype)
                    break;

            /// count until same type
            for(k=j+1,counter=1;k<atomsSize;k++,counter++){
                if(atoms[k].atype!=atype)
                    break;
            }

            atomNamesNumber[atype]=counter;
        }

}
//-----------------------------------------------------------------------------
double getRandomValue(const vector<string> &tokens)
{
std::default_random_engine generator (std::chrono::system_clock::now().time_since_epoch().count());
const position p0=std::stod(tokens[1]);
const position p1=std::stod(tokens[2]);

            if(tokens[0]=="uniform"){
            std::uniform_real_distribution<position> distribution(p0,p1);

            return distribution(generator);
            }

            if(tokens[0]=="normal"){
            std::normal_distribution<position> distribution(p0,p1);
            const double val=distribution(generator);
            const double tol=std::abs(val-p0);
            const double randValue=(tol<3*p1)? val: ((val-p0) ? p0+3*p1 : p0-3*p1 );

            return randValue;
            }


            if(tokens[0]=="lognormal"){
            std::lognormal_distribution<position> distribution(p0,p1);
            const double val=distribution(generator);

            return val;
            }

return 1.0;
}

//-----------------------------------------------------------------------------
position NanoGrain::StNanoGrain::getLP()
{
vector<string> tokens(split<string>(clp," "));

            if(tokens.size()==1) return std::stod(clp);
            else {
            const double randValue=getRandomValue(tokens);

                if(randValue<=0) throw ERR_LPLTZERO;

                return randValue;
            } //ERR_LPLTZERO,ERR_RLTZERO
}
//-----------------------------------------------------------------------------
position NanoGrain::StNanoGrain::getRadius(cpos &lp,const string &sradius)
{
const string R((sradius.empty()) ? this->radius : sradius);
vector<string> tokens(split<string>(R," "));

        if(tokens.size()==1){
        const size_t posLP=radius.find("lp");
        return (posLP==string::npos)? std::stod(R) : std::stod(R.substr(0,posLP))*lp;
        }
        else {
        const double randValue=getRandomValue(tokens);

                if(randValue<=0) throw ERR_RLTZERO;

        return randValue;
        }

}
//-----------------------------------------------------------------------------
int NanoGrain::StNanoGrain::findAtomName(const string &aname__)
{
int apos=0;

            for(auto &atype :atomTypes){
            std::string & aname(atype.name);
                if( aname==aname__ )
                return apos;

                apos++;
            }

return -1;
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::addAtomName(const string &aname__, bool checkMass)
{
const bool anamedup= (findAtomName(aname__)>=0);

        if(anamedup)
            cerr<<"WARNING: atom name duplicate has been found: "<<aname__<<endl;

        if(checkMass){
            if(Elements::mass.find(aname__)==Elements::mass.end())
                cerr<<"WARNING: unknown mass of an atom: "<<aname__<<endl;
        }

        atomTypes.emplace_back(StAtomType(aname__));
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::saveToFile()
{
const size_t nOfatoms=atoms.size();

            if(nOfatoms<saveopt.min){
                if(DB) cout<<__FILE__<<":"<<__LINE__<<" number of atoms < min "<<endl;
                saveopt.fileSaved=false;
            return;
            }

            if(nOfatoms>saveopt.max){
                if(DB) cout<<__FILE__<<":"<<__LINE__<<" number of atoms > max "<<endl;
                //testSNA=true;
                saveopt.fileSaved=false;
            return;
            }


            //cout<<"saveopt.fileSaved:  "<<saveopt.fileSaved<<endl;
            for(auto & ifc : saveopt.lwh){
                if(ifc=="H>W" && maxX>maxZ) return;
                if(ifc=="H>L" && maxY>maxZ) return;
            }



            if(DB) cout<<__FILE__<<":"<<__LINE__<<" saving atoms "<<endl;


            for (auto & fileName:fileNameOut){

                if(!createDirsIfDontExist(fileName)){
                    errMsg("couldn't create nested directories for "+fileName);
                throw Status::ERR_FOPEN;
                }

                if(fileName.rfind(".dat")!=string::npos){
                    saveDatFile(fileName);
                continue;
                }

                if(fileName.rfind(".nxyz")!=string::npos){
                    saveNXYZFile(fileName);
                continue;
                }

                if(fileName.rfind(".mxyz")!=string::npos){
                    saveMXYZFile(fileName);
                continue;
                }


                if(fileName.rfind(".xyz")!=string::npos){
                    saveXYZFile(fileName);
                    saveopt.fileSaved=true;
                continue;
                }



                if(fileName.rfind(".ndl")!=string::npos){
                    saveNDLFile(fileName);
                    saveopt.fileSaved=true;
                continue;
                }

                if(fileName.rfind(".lmp")!=string::npos){
                    saveLammpsFile(fileName);
                    saveopt.fileSaved=true;
                continue;
                }


                cerr<<" unknown file format "<<fileName<<endl;
            throw Status::ERR_FFORMAT;
            }


}


//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::saveDatFile(const string &fileName)
{


fstream fout(fileName,ios::out);

        if(!fout){
            cerr<<"couldn't open file for saving, atom positions will be lost"<<endl;
            fout.close();
        throw Status::ERR_FOPEN;
        }

        for(StAtom &atom: atoms)
            fout<<atom.x<<"\t"<<atom.y<<"\t"<<atom.z<<endl;

        fout.close();
}


//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::saveXYZFile(const string &fileName)
{
fstream fout(fileName,ios::out);

        if(!fout){
            cerr<<"couldn't open file for saving, atom positions will be lost"<<endl;
            fout.close();
        throw Status::ERR_FOPEN;
        }


        fout<<atoms.size()<<endl;
        fout<<"ver.01"<<endl;

        for(StAtom &atom: atoms)
             fout<<atomTypes[atom.atype].name<<"\t"<<atom.x<<"\t"<<atom.y<<"\t"<<atom.z<<endl;///<<"\t"<<( (atom.rdh) ? "*" : "" )

        fout.close();
}

//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::saveNXYZFile(const string &fileName)
{
fstream fout(fileName,ios::out);

        if(!fout){
            cerr<<"couldn't open file for saving, atom positions will be lost"<<endl;
            fout.close();
        throw Status::ERR_FOPEN;
        }


        fout<<atoms.size()<<endl;
        fout<<"ver.01"<<endl;

        //cout<<"atom.atype "<<atoms[0].atype<<endl;

        for(StAtom &atom: atoms)
             fout<<atom.id<<"\t"<<atomTypes[atom.atype].name<<"\t"<<atom.x<<"\t"<<atom.y<<"\t"<<atom.z<<endl;
            //fout<<atomNames[atom.atype]<<"\t"<<atom.x<<"\t"<<atom.y<<"\t"<<atom.z<<endl;

        fout.close();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::saveMXYZFile(const string &fileName)
{
string nxyzFileName(fileName);
        nxyzFileName[fileName.length()-4]='n';



fstream fout(nxyzFileName,ios::out);
fstream foutM(fileName,ios::out);


        if(!fout){
            cerr<<"couldn't open file for saving, atom positions will be lost"<<endl;
            fout.close();
        throw Status::ERR_FOPEN;
        }


        fout<<atoms.size()<<endl;
        fout<<"ver.01"<<endl;

        foutM<<"         "<<endl;
        foutM<<"selected atoms of "<<nxyzFileName<<endl;


size_t numOfMarkedAtoms=0;

        for(StAtom &atom: atoms){
             fout<<atom.id<<"\t"<<atomTypes[atom.atype].name<<"\t"<<atom.x<<"\t"<<atom.y<<"\t"<<atom.z<<endl;
             if(atom.mark){
                 foutM<<atom.id<<"\t"<<atomTypes[atom.atype].name<<"\t"<<atom.x<<"\t"<<atom.y<<"\t"<<atom.z<<endl;
                 numOfMarkedAtoms++;
             }
        }


        foutM.seekp(0,ios::beg);
        foutM<<numOfMarkedAtoms;



        fout.close();
        foutM.close();
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::saveNDLFile(const string &fileName)
{

fstream fout(fileName,ios::out);

        if(!fout){
            cerr<<"couldn't open file for saving, atom positions will be lost"<<endl;
            fout.close();
        throw Status::ERR_FOPEN;
        }

        fout<<"#:ver: 0"<<endl;
        fout<<(*this);

        for(StAtom &atom: atoms)
            fout<<atomTypes[atom.atype].name<<"\t"<<atom.x<<"\t"<<atom.y<<"\t"<<atom.z<<endl;

        fout.close();
}
//-----------------------------------------------------------------------------
inline position minv(position a,position b,position c)
{
position minv=0;

            if(a<minv) minv=a;
            if(b<minv) minv=b;
            if(c<minv) minv=c;

return minv;
}
//-----------------------------------------------------------------------------
inline position maxv(position a,position b,position c)
{
position maxv=0;

            if(a>maxv) maxv=a;
            if(b>maxv) maxv=b;
            if(c>maxv) maxv=c;

return maxv;
}
//-----------------------------------------------------------------------------

void NanoGrain::StNanoGrain::saveLammpsFile(const string &fileName)
{

fstream file(fileName,ios::out);

            if(!file){
                cerr<<"couldn't open file for saving, atom positions will be lost"<<endl;
                file.close();
            throw Status::ERR_FOPEN;
            }


//cdouble clpH=std::stod(clp)*0.5;

double mX,mY,mZ;
StMinMax minmax;


            if(margins.empty()){
                if(clp.empty()){
                    if(uc.vtrans.empty()){
                        mX=mY=mZ=0;
                    }
                    else{
                        mX=0.5*(uc.vtrans[0].x+uc.vtrans[1].x+uc.vtrans[2].x);
                        mY=0.5*(uc.vtrans[0].y+uc.vtrans[1].y+uc.vtrans[2].y);
                        mZ=0.5*(uc.vtrans[0].z+uc.vtrans[1].z+uc.vtrans[2].z);
                    }
                }
                else
                    mX=mY=mZ=0.5*std::stod(clp);
            }
            else{
            vector<string> toks(split<string>(margins," "));
                if(toks.size()==1)
                    mX=mY=mZ=std::stod(toks[0]);
                else{
                    mX=std::stod(toks[0]);
                    mY=std::stod(toks[1]);
                    mZ=std::stod(toks[2]);
                }
            }



            minmax.searchForMinMax(atoms);

            file<<"# created by npcl"<<endl<<endl;

            //------------------------
            file<<"\t"<<atoms.size()<<"\tatoms"<<endl;
            file<<"\t"<<atomTypes.size()<<"\tatom types"<<endl;



            if(saveopt.lmpTric){
                if(structure.find("tric")!=str::npos || structure.find("cif")!=str::npos){
                 #pragma message (" !!! definition of triclinic box is not implemented for xz, yz planes")
                const int repX=std::stoi(replicate[0]);
                const int repY=std::stoi(replicate[1]);
                const int repZ=std::stoi(replicate[2]);

                cpos XY=repY*tric.b.x;

                cpos XZ=repZ*tric.c.x;
                cpos YZ=repZ*tric.c.y;

                cpos XYZ=XY+XZ;


                cpos xlo=0-minv(XY,XZ,XYZ);
                cpos xhi=repX*tric.a.x-maxv(XY,XZ,XYZ);

                cpos ylo=0-minv(0,0,YZ);
                cpos yhi=repY*tric.b.y-maxv(0,0,YZ);

                    file<<"    "<<0 <<"    "<<repX*tric.a.x*1.25<<"    xlo xhi"<<endl;
                    file<<"    "<<-0.125* repY*tric.b.y  <<"    "<<repY*tric.b.y*1.25<<"    ylo yhi"<<endl;
                    file<<"    "<<0 <<"    "<<repZ*tric.c.z*1.25<<"    zlo zhi"<<endl;

                    file<<"    "<<XY<<"    "<<XZ<<"    "<<YZ<<"   xy xz  yz"<<endl;
                }            
                else{
                    file<<"    "<<minmax.xmin-mX<<"    "<<minmax.xmax+mX<<"    xlo xhi"<<endl;
                    file<<"    "<<minmax.ymin-mY<<"    "<<minmax.ymax+mY<<"    ylo yhi"<<endl;
                    file<<"    "<<minmax.zmin-mZ<<"    "<<minmax.zmax+mZ<<"    zlo zhi"<<endl;
                    file<<"    0  0  0  xy  xz  yz"<<endl;
                }
            }
            else{
                file<<"    "<<minmax.xmin-mX<<"    "<<minmax.xmax+mX<<"    xlo xhi"<<endl;
                file<<"    "<<minmax.ymin-mY<<"    "<<minmax.ymax+mY<<"    ylo yhi"<<endl;
                file<<"    "<<minmax.zmin-mZ<<"    "<<minmax.zmax+mZ<<"    zlo zhi"<<endl;
            }



            //------------------------
            file<<"Masses"<<endl<<endl;

mapConstIter massIter;

            for( size_t i=0;i<atomTypes.size();i++){

                file<<"    "<<i+1<<"  ";

                if(atomTypes[i].mass.empty()){
                    if( (massIter=Elements::mass.find(atomTypes[i].name))==Elements::mass.end()){
                        file<<"?";
                        warnMsg("unknown mass for "+atomTypes[i].name);
                    }
                    else
                        file<<massIter->second;
                }
                else
                    file<<atomTypes[i].mass;

                file<<"  #  "<<atomTypes[i].name<<endl;

            }
            //------------------------

size_t i;
const size_t numOfatomsTot=atoms.size();
const size_t at=atomTypes.size();
const StAtom * ptr_atom=atoms.data();

            if(lmpstyle=="atomcharge"){

                file<<" Atoms # charge"<<endl<<endl;

                for(i=1;i<=numOfatomsTot;i++,ptr_atom++){
                    file<<setw(6)<<i<<"  "<<( ptr_atom->atype+1 )<<"  "<<  atomTypes[ptr_atom->atype].charge<<"  "
                        <<ptr_atom->x<<" "<<ptr_atom->y<<" "<<ptr_atom->z<<endl;
                }
            }
            else{

                file<<" Atoms"<<endl<<endl;

                for(i=1;i<=numOfatomsTot;i++,ptr_atom++)
                    file<<setw(6)<<i<<"  "<<  ptr_atom->atype+1 <<"  "<<ptr_atom->x<<" "<<ptr_atom->y<<" "<<ptr_atom->z<<endl;
            }


            file.close();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::saveHeader()
{
            if(!createDirsIfDontExist(fileNameHeader)){
                errMsg("couldn't create nested directories for "+fileNameHeader);
            throw Status::ERR_FOPEN;
            }

fstream fout(fileNameHeader,ios::out);

            if(!fout){
                cerr<<"couldn't open file for saving, header will be lost: "<<fileNameHeader<<endl;
                fout.close();
            throw Status::ERR_FOPEN;
            }

            fout<<"#:ver: 1"<<endl;
            fout<<(*this);



            fout.close();
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::openFile()
{

            if(fileNameIn.find(".xyz")!=string::npos){
                openXYZFile();
                sortAtomsByName();
            return;
            }

            if(fileNameIn.find(".lmp")!=string::npos){
                openLMPFile();
                sortAtomsByName();
            return;
            }

            if(fileNameIn.find(".ndl")!=string::npos){
                openNDLFile();
                sortAtomsByName();
            return;
            }


            cerr<<" unknown format "<<endl;
            throw Status::ERR_FFORMAT;
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::openXYZFile()
{
fstream fin(fileNameIn,ios::in);

            if(!fin){
                cerr<<"couldn't open file for reading"<<endl;
                fin.close();
                throw Status::ERR_FOPEN;
            }

int row=0,arows=-1;

            try{
                    fin.exceptions(ifstream::failbit | ifstream::badbit | ifstream::eofbit);

                    fin>>arows;
                    while(fin.get()!='\n' ) ;


                    if(arows<0){
                        cerr<<" rows <0"<<endl;
                    throw Status::ERR_FFORMAT;
                    }


                    while(fin.get()!='\n' ) ;  //ignore a comment row

                    atoms.clear();
                    atoms.reserve(arows);

            position x,y,z;
            string aname;
            int atype;

                    fin.exceptions(ifstream::failbit | ifstream::badbit);

                    for(  ;row<arows;row++){
                        fin>>aname>>x>>y>>z;

                        if( (atype=findAtomName(aname)) <0){
                            //atomNames.push_back(aname);
                            atomTypes.push_back(StAtomType(aname));
                            atype=atomTypes.size()-1;
                        }

                        atoms.emplace_back(StAtom(x,y,z,atype));
                    }

                    atoms.shrink_to_fit();
            }
            catch(std::ifstream::failure e){

                cerr<<" exceptions during file procesing, row: "<<row<<", e.what(): "<<e.what()<<endl;
                fin.close();
                throw Status::ERR_FIN_EXC;
            }
            catch(Status e){

                fin.close();
                throw e;
            }


fin.close();
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::openLMPFile()
{
fstream fin(fileNameIn,ios::in);

            if(!fin){
                cerr<<"couldn't open file for reading"<<endl;
                fin.close();
                throw Status::ERR_FOPEN;
            }

int nOfatoms,nOftypes;
string sline;

            try{
                    fin.exceptions(ifstream::failbit | ifstream::badbit | ifstream::eofbit);

                    //ignore 2 lines
                    gotoEOL(fin);
                    gotoEOL(fin);


                    //read number of atoms
                    std::getline(fin,sline);
             vector<string> tokens(split<string>(sline," \t"));

                    if(tokens.size()<2 && tokens[1]!="atoms"){
                        if(DB) cerr<<__FILE__<<":"<<__LINE__<<"  ";
                        cerr<<" unknown number of atoms";
                    throw Status::ERR_FFORMAT;
                    }

                    nOfatoms=std::stoi(tokens[0]);


                    //read number of types
                    std::getline(fin,sline);
                    tokens=split<string>(sline," \t");

                    nOftypes=std::stoi(tokens[0]);

                    atomTypes.resize(nOftypes);
                    for (size_t i=1;i<=nOftypes;i++)
                        atomTypes[i-1].name=std::to_string(i);

                    // find tag 'Masses'
                    do{
                        std::getline(fin,sline);
                    }while(sline.find("Masses")==string::npos );

                    int countTypes=0;
                    do{
                        std::getline(fin,sline);
                        tokens=split<string>(sline," \t");
                        if(tokens.empty()) continue;

                        countTypes++; 

                    }while(countTypes<nOftypes);


                    // find tag 'Atoms'
                    do{
                        std::getline(fin,sline);
                    }while(sline.find("Atoms")==string::npos );

                    // ignore empty line
                    gotoEOL(fin);


                    // READ atom positions
                    atoms.clear();
                    atoms.reserve(nOfatoms);

            position x,y,z;
            string nofa;
            int id;

                    fin.exceptions(ifstream::failbit | ifstream::badbit);


                    if(nOfatoms<1e6){
                        for(int row=0 ;row<nOfatoms;row++){
                            fin>>nofa>>id>>x>>y>>z;
                            gotoEOL(fin);

                            atoms.emplace_back(StAtom(x,y,z,id-1));
                        }
                    }
                    else{
                    CProgress progress;
                            progress.title=(string(" lmp file reading "));
                            progress.start(nOfatoms);

                            for(int row=0 ;row<nOfatoms;row++,progress++){
                                fin>>nofa>>id>>x>>y>>z;
                                gotoEOL(fin);

                                atoms.emplace_back(StAtom(x,y,z,id-1));
                            }

                            progress.stop();
                            cerr<<endl;
                    }

                    atoms.shrink_to_fit();

            }
            catch(std::ifstream::failure &e){
                cerr<<" exception during file procesing, e.what(): "<<e.what()<<endl;
                fin.close();
                throw Status::ERR_FIN_EXC;
            }
            catch(Status &e){
                fin.close();
                throw e;
            }

            fin.close();
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::openNDLFile()
{
fstream fin(fileNameIn,ios::in);

            if(!fin){
                cerr<<"couldn't open file for reading"<<endl;
                fin.close();
            throw Status::ERR_FOPEN;
            }

int nOfatoms,nOftypes;
string sline;

            try{
                    fin.exceptions(ifstream::failbit | ifstream::badbit | ifstream::eofbit);

                    while(fin.peek()=='#'){

                        std::getline(fin,sline);
                    vector<string> toks{split<string>(sline," \t\r")};

                        if(toks[0].find("#sizeRC")!=string::npos){
                            nOfatoms=std::stoi(toks[1]);
                        continue;
                        }

                        if(toks[0].find("#atoms")!=string::npos){
                            for(size_t i=1;i<toks.size();i++){
                                if(findAtomName(toks[i])<0)
                                   atomTypes.push_back(StAtomType(toks[i]));
                            }                        
                        }

                    }

                    atoms.clear();
                    atoms.reserve(nOfatoms);

            position x,y,z;
            string aname;
            int atype;

                    fin.exceptions(ifstream::failbit | ifstream::badbit);

                    for(int i=0;i<nOfatoms;i++){
                        fin>>aname>>x>>y>>z;

                        if( (atype=findAtomName(aname)) <0){
                            atomTypes.push_back(StAtomType(aname));
                            atype=atomTypes.size()-1;
                            warnMsg("possible wrong file format");
                        }

                        atoms.emplace_back(StAtom(x,y,z,atype));
                    }

                    atoms.shrink_to_fit();

            }
            catch(std::ifstream::failure &e){
                cerr<<" exception during file procesing, e.what(): "<<e.what()<<endl;
                fin.close();
                throw Status::ERR_FIN_EXC;
            }
            catch(Status &e){
                fin.close();
                throw e;
            }

            fin.close();

}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::sortAtomsByName()
{
auto fsort=[](StAtom &a, StAtom&b) { return a.atype<b.atype; } ;
           std::sort(atoms.begin(),atoms.end(),fsort);

}
//-----------------------------------------------------------------------------
bool NanoGrain::StNanoGrain::parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__)
{
const str send("end");

        try{
        resetPrms();
        ptr_uvar=uvar__;

        while(cmd[index]!=send){

                if(DB) cout<<cmd[index]<<endl;


                if(cmd[index]=="atom"){
                    addAtomName(cmd[index++][1]);
                continue;
                }

                if(cmd[index]=="atoms"){
                const int kvsize=cmd[index].numOfKeyValues();

                    for(int i=1;i<kvsize;i++)
                        addAtomName(cmd[index][i]);

                    index++;
                continue;
                }

                if(cmd[index]=="center"){
                    if(cmd[index][1]=="geom")
                        center=GEOM;
                    else
                        if(cmd[index][1]=="catom")
                            center=CATOM;
                        else
                            center=ID;

                    if(cmd[index].numOfKeyValues()==3)
                        catomType=cmd[index][2];

                    index++;
                continue;
                }

                if(cmd[index]=="charge"){
                const string aname(cmd[index][1]);
                const string value(cmd[index][2]);
                const int anamedup= (findAtomName(aname));

                        lmpstyle="atomcharge";

                        if(anamedup>=0)
                            atomTypes[anamedup].charge=value;
                        else{
                            addAtomName(aname,false);
                            atomTypes.back().charge=value;
                        }

                    index++;
                continue;
                }


                if(cmd[index]=="csh"){

                    do{
                        index++;

                        if(cmd[index]=="end")
                        break;

                    vector<string> tokens(split<string>(cmd[index][0]," "));//scale factors

                        if(tokens.size()==2)
                            coreshell.prm.emplace_back(CShells::shell(tokens[0],tokens[1]));
                        else {
                        CrandomUni dev(std::stod(tokens[1]),std::stod(tokens[2]));
                        CrandomUni radius(std::stod(tokens[4]),std::stod(tokens[5]));

                            coreshell.prm.emplace_back(CShells::shell(dev.randNumber(),radius.randNumber()));
                        }

                    }while(true);


                    /// additional shell to simplify coreShell algorithm  placed in newR
                    coreshell.prm.emplace_back(CShells::shell("0","inf"));

                    index++;

                continue;
                }

                if(cmd[index]=="disloc"){

                    if(cmd[index][1]=="plane")
                        dislocPlane=cmd[index][2];
                    else
                       disloc=cmd[index][1];

                    index++;
                continue;
                }

                if(cmd[index]=="disperse"){
                string name= cmd [index][1];
                string value=cmd [index][2];

                    atomDisperse.emplace_back(StAtomDispFactor(name,value));

                    disperse=true;
                    index++;
                continue;
                }


                if(cmd[index]=="geometry"){
                    shape=cmd[index][1];

                string varValue;
                    for(int i=2;i<cmd[index].numOfKeyValues();i++){
                        varValue=cmd[index][i];
                        Script::replaceVars(ptr_uvar,varValue);
                        shapePrm+=varValue;
                    }

                    if(DB){ cout<<__FILE__<<":"<<__LINE__<<"  "<<shapePrm<<endl;}

                    index++;
                continue;
                }


                if(cmd[index]=="lp" ){
                    clp=cmd[index++][1];
                    Script::replaceVars(ptr_uvar,clp);
                continue;
                }


                if(cmd[index]=="lmpstyle"){
                    lmpstyle=cmd[index++][1];
                continue;
                }

                
                if(cmd[index]=="lmpmargin"){
                    margins=cmd[index++][1];
                continue;
                }
                
                
                if(cmd[index]=="mass"){
                const string aname(cmd[index][1]);
                const string value(cmd[index][2]);
                const int anamedup= (findAtomName(aname));

                        if(anamedup>=0)
                            atomTypes[anamedup].mass=value;
                        else{
                            addAtomName(aname,false);
                            atomTypes.back().mass=value;
                        }

                        index++;
                continue;
                }


                if(cmd[index]=="hcpabc"){

                    if(cmd[index][1]=="abc"){
                    string varValue(cmd[index][2]);
                        Script::replaceVars(ptr_uvar,varValue);
                    const string abcexpr(parseABC(varValue));
                        if(DB) cout<<abcexpr<<endl;
                        hcpABC=abcexpr;
                    }
                    else{
                        for(int i=1;i<cmd[index].numOfKeyValues();i++)
                            hcpABC+=" "+cmd[index][i];
                    }

                    index++;
                continue;
                }


                if(cmd[index]=="hcpu"){
                    hcpu=cmd[index++][1];
                continue;
                }


                if(cmd[index]=="hcpcs"){
                    hcpcs=cmd[index][1];

                    if(hcpcs=="poly"){
                    string varValue;

                        for(int i=2;i<cmd[index].numOfKeyValues();i++){
                            varValue=cmd[index][i];
                            Script::replaceVars(ptr_uvar,varValue);
                            shapePrm2D+=" "+varValue;
                        }
                    }

                    index++;
                continue;
                }


                if(cmd[index]=="hcpfillup"){
                    hcpfillup=cmd[index][1];
                    if(cmd[index].numOfKeyValues()==3)
                        hcplaynum=cmd[index][2];
                    else
                        hcplaynum.clear();

                    index++;
                continue;
                }


                if(cmd[index]=="hcpsl"){
                    if(cmd[index].numOfKeyValues()==1)
                        hcpsl=true;
                    else
                        hcpsl=(cmd[index][1]=="yes");

                    index++;
                continue;
                }

                if(cmd[index]=="hcpsurfA"){
                auto val=cmd[index][1];
                const vector<string> opt{"yes","no","AB","H"};
                auto it=std::find_if(opt.begin(),opt.end(),[&val](const string &x){return x==val;});

                     hcpsurf=(EHCPSURF)(std::distance(opt.begin(),it));

                    index++;
                continue;
                }


                if(cmd[index]=="hcpsurftype"){
                auto val=cmd[index][1];

                     Script::replaceVars(ptr_uvar,val);

                     switch ( val[0]) {
                     case '0': hcpsurf=EHCPSURF::sA; break;
                     case '1': hcpsurf=EHCPSURF::sAB; break;
                     case '2': hcpsurf=EHCPSURF::sB; break;
                     default: throw Status::ERR_ST;
                     }

                    index++;
                continue;
                }


                if(cmd[index]=="insfault"){
                ClKeyValues keyValues(cmd[index]);

                    if(keyValues.getValue(2).find("$")!=string::npos)
                        Script::replaceVars(ptr_uvar,keyValues.getValue(2));

                    faultmode= (keyValues.getValue(1)=="random") ? StNanoGrain::RANDOM : StNanoGrain::CUSTOM;
                    hcpfault = keyValues.getValue(2);

                    index++;
                continue;
                }


                if(cmd[index]=="cradius"){
                    cradius=cmd[index++][1];
                    Script::replaceVars(ptr_uvar,cradius);
                continue;
                }


                if(cmd[index]=="numOfatomsTest"){
                    numOfAtomsTest=(cmd[index++][1]=="yes");
                continue;
                }


                if(cmd[index]=="open"){
                string fileName(cmd[index++][1]);

                    Script::replaceVars(ptr_uvar,fileName);
                    if(fileName.rfind(".cif")!=string::npos)
                        fileCIF=fileName;
                    else
                        fileNameIn=fileName;


                continue;
                }


                if(cmd[index]=="push"){
                    if(cmd[index].getValue(1)=="numAtoms"){
                        numOfAtomsPush=cmd[index].getValue(2);
                        index++;
                    continue;
                    }

                    if(cmd[index].getValue(1)=="saveStatus"){
                        saveFileStatusPush=cmd[index].getValue(2);
                    }

                    index++;
                continue;
                }

                if(cmd[index]=="printPrm" || cmd[index]=="printprm"){
                    dispParams=true;

                    index++;
                continue;
                }


                if(cmd[index]=="radius"){
                    radius=cmd[index++][1];
                    Script::replaceVars(ptr_uvar,radius);                    
                continue;
                }


                if(cmd[index]=="replicate"){
                ClKeyValues keyValues(cmd[index]);

                    Script::replaceVars(ptr_uvar,keyValues.getValue(1));
                    Script::replaceVars(ptr_uvar,keyValues.getValue(2));
                    Script::replaceVars(ptr_uvar,keyValues.getValue(3));

                    replicate.resize(keyValues.numOfKeyValues()-1);
                    replicate[0]=keyValues.getValue(1);
                    replicate[1]=keyValues.getValue(2);
                    replicate[2]=keyValues.getValue(3);                    

                    index++;
                continue;
                }

                if(cmd[index]=="rename"){
                std::string val{cmd[index].getValue()};

                    Script::replaceVars(ptr_uvar,val);
                    rename.emplace_back(val);

                    index++;
                continue;
                }

                if(cmd[index]=="rescale"){
                    scaleFactors=cmd[index][1];
                    if(cmd[index].numOfKeyValues()==(1+3) )
                        scaleFactors+=" "+cmd[index][2]+" "+cmd[index][3];

                    index++;
                continue;
                }


                if(cmd[index]=="remove"){
                    if(cmd[index].numOfKeyValues()==2){
                    std::string rparams{cmd[index][1]};

                        Script::replaceVars(ptr_uvar,rparams);
                        vremoveAtomsPrms.emplace_back(rparams);
                    }
                    else{
                        rmatoms=cmd[index][1]+" "+cmd[index][2];
                        if(cmd[index].numOfKeyValues()==(1+3))
                            rmatoms+=" "+cmd[index][3];
                    }

                    index++;
                continue;
                }

                if(cmd[index]=="rotate"){
                string varValue;

                    for(int i=1;i<cmd[index].numOfKeyValues();i++){
                        varValue=cmd[index][i];
                        Script::replaceVars(ptr_uvar,varValue);
                        rotatePrm+=" "+varValue;
                    }

                    index++;
                continue;
                }


                if(cmd[index]=="save"){
                std::string fileName =cmd[index++][1];
                    Script::replaceVars(ptr_uvar,fileName);

                    fileNameOut.emplace_back(fileName);
                continue;
                }


                if(cmd[index]=="saveopt"){

                    for(int i=1; i<cmd[index].numOfKeyValues(); i+=2){
                    const string key  (cmd[index][i]);
                    const string value(cmd[index][i+1]);

                        if(key=="min"){
                            saveopt.min=std::stod(value);
                        continue;
                        }

                        //StRotationMatrix rotMat(mainAxis,sa,ca);



                        if(key=="max"){
                            saveopt.max=std::stod(value);
                        continue;
                        }

                        if(key=="if"){
                            saveopt.lwh.push_back(value);
                        continue;
                        }

                        if(key=="lmp"){
                            saveopt.lmpTric=true;
                        continue;
                        }

                        warnMsg("unknown save option "+key);
                    }

                    index++;
                continue;
                }


                if(cmd[index]=="saveHeader"){
                    fileNameHeader =cmd[index++][1];
                    Script::replaceVars(ptr_uvar,fileNameHeader);
                continue;
                }


               /* if(cmd[index]=="side"){
                   // fileNameOut=cmd[index++][1];
                    //Script::replaceVars(ptr_uvar,fileNameOut);
                continue;
                }*/

                if(cmd[index]=="struct"){
                    structure=cmd[index++][1];
                continue;
                }

                if(cmd[index]=="threads"){
                    threads=cmd[index++][1];
                continue;
                }

                if(cmd[index]=="voids"){
                    voidsprm=cmd[index][1];
                    if(cmd[index].numOfKeyValues()==3)
                        voidsprm+=" "+cmd[index][2];

                    index++;
                continue;
                }

                if(cmd[index]=="tricp"){

                    tric.atoms.reserve(1);

                    do{
                        index++;

                        if(cmd[index]=="end")
                        break;

                        if(cmd[index]=="csys"){
                            if(cmd[index].getValue()=="tric")
                                tric.csys=StTric::TRIC;
                            else
                                tric.csys=StTric::CART;

                        continue;
                        }

                    ClKeyValues keyValues(cmd[index]);

                        Script::replaceVars(ptr_uvar,keyValues.getValue(1));

                        if( ( std::string("lpa lpb lpc").find(cmd[index].getKey())) !=std::string::npos ){
                        string lpabc(cmd[index][1]);

                           Script::replaceVars(ptr_uvar,lpabc);
                           switch(cmd[index].getKey()[2]){
                           case 'a': tric.lpa=lpabc;break;
                           case 'b': tric.lpb=lpabc;break;
                           case 'c': tric.lpc=lpabc;break;
                           }

                        continue;
                        }


                        if( ( std::string("alpha beta gamma").find(cmd[index].getKey())) !=std::string::npos  ){
                        string abg(cmd[index][1]);

                            Script::replaceVars(ptr_uvar,abg);
                            switch(cmd[index].getKey()[0]){
                            case 'a': tric.alpha=abg;break;
                            case 'b': tric.beta=abg;break;
                            case 'g': tric.gamma=abg;break;
                            }

                        continue;
                        }

                        Script::replaceVars(ptr_uvar,keyValues.getValue(2));
                        Script::replaceVars(ptr_uvar,keyValues.getValue(3));

                        tric.atoms.push_back(StUcAtom());
                        tric.atoms.back().name=cmd[index].getKey();
                        tric.atoms.back().x=std::stod(keyValues.getValue(1));
                        tric.atoms.back().y=std::stod(keyValues.getValue(2));
                        tric.atoms.back().z=std::stod(keyValues.getValue(3));

                        // rdh/random mode atoms
                        auto rdhKV=(keyValues.numOfKeyValues()-4);
                        if(rdhKV<2){ // * atoms
                            tric.atoms.back().rdh=(bool)rdhKV;
                            tric.rdhAtoms+=rdhKV; }
                        else{ // rm atoms
                            tric.atoms.back().rmProb=std::stod(keyValues.getValue(5));
                        }
                        //

                    }while(true);

                    index++;

                continue;
                }


                if(cmd[index]=="ucp"){

                    uc.vtrans.resize(3);
                    uc.atoms.reserve(1);

                    do{
                        index++;

                        if(cmd[index]=="end")
                        break;

                    ClKeyValues keyValues(cmd[index]);

                        Script::replaceVars(ptr_uvar,keyValues.getValue(1));
                        Script::replaceVars(ptr_uvar,keyValues.getValue(2));
                        Script::replaceVars(ptr_uvar,keyValues.getValue(3));

                        if(cmd[index].getKey()[0]=='v'){
                        auto  ptr=uc.vtrans.begin()+( cmd[index].getKey()[1]-'x');

                            ptr->x=std::stod(keyValues.getValue(1));
                            ptr->y=std::stod(keyValues.getValue(2));
                            ptr->z=std::stod(keyValues.getValue(3));

                        continue;
                        }

                        uc.atoms.push_back(StUcAtom());
                        uc.atoms.back().name=cmd[index].getKey();
                        uc.atoms.back().x=std::stod(keyValues.getValue(1));
                        uc.atoms.back().y=std::stod(keyValues.getValue(2));
                        uc.atoms.back().z=std::stod(keyValues.getValue(3));

                        // rdh mode atoms
                        auto rdhKV=(keyValues.numOfKeyValues()-4);
                        if(rdhKV<2){ // * atoms
                            uc.atoms.back().rdh=(bool)rdhKV;
                            uc.rdhAtoms+=rdhKV; }
                        else{ // rm atoms
                            uc.atoms.back().rmProb=std::stod(keyValues.getValue(5));
                        }
                        //

                    }while(true);

                    index++;
                continue;
                }


                if(cmd[index]=="comment"){
                    comment=cmd[index++][1];
                continue;
                }

                cerr<<"Error: nanograin.cpp  unknown command "<< cmd[index]<<" line "<<__LINE__<<endl;
                return false;
        }



            if(fileNameIn.empty()){
            const bool vttriccif=uc.vtrans.empty() &
                                    tric.empty() &
                                    fileCIF.empty();

                if(clp.empty() && vttriccif)                throw Status::ERR_LP;
                if(radius.empty() && replicate.empty() )    throw Status::ERR_RADII;
                if(shape.empty()     && vttriccif)          throw Status::ERR_GEOM;
                if(atomTypes.empty() && vttriccif)          throw Status::ERR_ATYPES;

                build();

                if(!coreshell.prm.empty()) displaceAtoms();

            }
            else
                openFile();


            /// -------------- common operations --------------------

            if(atoms.empty()) throw Status::ERR_ANUM_ZERO;

            if(!voidsprm.empty()) voids();
            if(!rmatoms.empty() || !vremoveAtomsPrms.empty() ) removeAtoms();

            if(!scaleFactors.empty()) rescale();
            if(disperse) disperseN();
            if(!rename.empty()) renameAtoms();

            findAtomNamesNumber();

            switch(center){
            case GEOM:  grainCentering();break;
            case CATOM: catomCentering();break;
            case ID:    idCentering();break;
            case COFF:break;
            }

            if(!rotatePrm.empty()) rotate();

            findMaxR();

            if(!numOfAtomsPush.empty()){
            const string varName(numOfAtomsPush);
            auto iterVar=std::find(ptr_uvar->begin(),ptr_uvar->end(),varName);

                    if(iterVar==ptr_uvar->end())
                        ptr_uvar->emplace_back(strpair(varName,std::to_string(atoms.size())));
                    else
                        iterVar->getValue()=std::to_string(atoms.size());
            }


            if(!dislocPlane.empty())
                insertDisloc();


            if(numOfAtomsTest){
                if(DB) infoMsg(" testing number of atoms");

                // if true  ->  grain with the same number of atoms already exists                               
                if(testSavedNumOfAtoms(atoms.size())){
                    if(DB){
                        const string warn(" grains with the same number of atoms are ignored "+fileNameOut[0]+" "+std::to_string(atoms.size()));
                        warnMsg(warn);
                    }
                    saveopt.fileSaved=false;
                   
                }
                else{
                    savedNumOfAtoms.push_back(atoms.size());
                    if(!fileNameOut.empty()) {
                        saveToFile();
                    }
                }
            }
            else
                if(!fileNameOut.empty()) saveToFile();


            if(!saveFileStatusPush.empty()){
            const string varName(saveFileStatusPush);
            auto iterVar=std::find(ptr_uvar->begin(),ptr_uvar->end(),varName);
            auto sval(std::to_string(saveopt.fileSaved));

                    if(iterVar==ptr_uvar->end())
                        ptr_uvar->emplace_back(strpair(varName,sval));
                    else
                        iterVar->getValue()=sval;
            }

            if(!fileNameHeader.empty() && saveopt.fileSaved ) saveHeader();

            if(dispParams) dispParameters();

         return true;
         }
         catch(Status e){

             cerr<<"ERROR ";
             switch(e){
             case Status::ERR_FFORMAT :  cerr<<" wrong format";               break;
             case Status::ERR_FOPEN :    cerr<<" file couldn't be opened";    break;
             case Status::ERR_LP :       cerr<<" unknown lattice parameter";  break;
             case Status::ERR_RADII :    cerr<<" unknown radius  ";           break;
             case Status::ERR_GEOM  :    cerr<<" unknown geometry ";          break;
             case Status::ERR_ATYPES :   cerr<<" unknown atom types ";        break;
             case Status::ERR_ANUM_ZERO: cerr<<" number of atoms is 0";       break;
             case NanoGrain::ERR_UNKSEQ: cerr<<" empty hcp string";           break;
             case Status::ERR_LPLTZERO:  cerr<<" random value of latt. param. is less than zero"; break;
             case Status::ERR_RLTZERO:   cerr<<" random value of radius is less than zero"; break;
             case Status::ERR_ABCPARSE:  cerr<<" hcpabc failure";               break;
             case Status::ERR_ST:        cerr<<" unrecognized surface type"; break;
             default: cout<<"unknown error (Status type), code "<<e<<" ";
             }

             cout<<endl;

         }
         catch (const std::out_of_range& e) {
           std::cerr << "ERROR: Out of Range error: " << e.what() << endl;
         }
         catch(Script::ResultRepVar e){
            throw e;
        }
         catch(...){
             cerr<<"ERROR: type unknown (grain building)"<<endl;
         }

return false;
}
//-----------------------------------------------------------------------------
void NanoGrain::StNanoGrain::dispParameters()
{
            cout<<"*** nanograin parameters ***"<<endl;
            cout<<(*this);
}
//-----------------------------------------------------------------------------
ostream &NanoGrain::operator<<(ostream &o, NanoGrain::StNanoGrain &grain)
{

        o<<"#clp: "<<grain.clp<<endl;
        o<<"#maxSize: "<<grain.maxR<<endl;
        o<<"#shape: "<<grain.shape<<endl;
        o<<"#structure: "<<grain.structure<<endl;

        if(grain.structure=="hcp"){
            o<<"#layers: "<<grain.hcpABC.size(); if(grain.hcpABC.size()<=120) o<<" "<<grain.hcpABC;
            o<<endl;

            o<<"#stackingFaults: ";
            switch(grain.faultmode){
            case NanoGrain::StNanoGrain::OFF    : o<<"off";    break;
            case NanoGrain::StNanoGrain::RANDOM : o<<"random"; for(size_t &s: grain.hcpFaultPos) o<<" "<<s;  break;
            case NanoGrain::StNanoGrain::CUSTOM : o<<"custom"; break;
            }
            o<<endl;

            o<<"#hcpU: "<<grain.hcpu<<endl;
            o<<"#hcpcs: "<<grain.hcpcs<<endl;
            o<<"#hcpSubLattice: "<<grain.hcpsl<<endl;
            o<<"#shapeParams: "<< ( (grain.ssShape)? grain.ssShape->getInfo() : "<empty set>") <<endl;
        }

        o<<"#inputFile: "<<grain.fileNameIn<<endl;
        o<<"#outputFile: "<<grain.fileNameOut.size();

            for(auto &fileName: grain.fileNameOut)
                o<<"    "<<fileName;
        o<<endl;


        o<<"#numberAtomsForEachType:"; for(size_t i=0;i<grain.atomTypes.size();i++) o<<" "<<grain.atomTypes[i].name<<" "<<grain.atomNamesNumber[i];
        o<<endl;


        o<<"#atomsRemoved: "<<grain.atomsRemoved.str()<<endl;

        o<<"#totalNumberOfAtoms: "<<grain.atoms.size()<<endl;

        if(!grain.comment.empty())
            o<<"#"<<grain.comment<<endl;


        if(grain.coreshell.prm.empty())
            o<<"#coreshell: 0"<<endl;
        else {
            o<<"#coreshell: "<<grain.coreshell.prm.size()-1<<endl;
            for(size_t i=1;i<grain.coreshell.prm.size();i++)
                o<<"#           "<<(grain.coreshell.prm[i-1].dev-1)*100<<"   "<<grain.coreshell.prm[i-1].radius<<endl;
        }



return o;
}
//---------------------------------------------------------------------------
void NanoGrain::StMinMax::searchForMinMax(const NanoGrain::vatoms &atoms)
{
        xmin=ymin=zmin=std::stod("inf");
        xmax=ymax=zmax=-xmin;

        for(auto & atom : atoms){
            if(atom.x>xmax) xmax=atom.x;
            if(atom.y>ymax) ymax=atom.y;
            if(atom.z>zmax) zmax=atom.z;

            if(atom.x<xmin) xmin=atom.x;
            if(atom.y<ymin) ymin=atom.y;
            if(atom.z<zmin) zmin=atom.z;
        }
}
//---------------------------------------------------------------------------
position NanoGrain::StMinMax::getMaxAbs()
{
auto tmp=(xmax>ymax) ? xmax : ymax;
return (tmp>zmax) ? tmp : zmax;
}
//---------------------------------------------------------------------------
position NanoGrain::StMinMax::getMinAbs()
{
auto tmp=(xmin<ymin) ? xmin : ymin;
return (tmp<zmin) ? tmp : zmin;
}
//---------------------------------------------------------------------------
position NanoGrain::StMinMax::getMinPos()
{
auto tmp=(xmax<ymax) ? xmax : ymax;
return (tmp<zmax) ? tmp : zmax;
}
//---------------------------------------------------------------------------
position NanoGrain::StNanoGrain::CShells::newR(const position &r)
{
size_t prmI;
position dr=0;
position rshellP=0;
position rshellN=0;
const size_t sizeInner=prm.size();

        for(prmI=0;prmI<sizeInner;prmI++){

            if(r<prm[prmI].radius){
                dr=(r-rshellP)*prm[prmI].dev;
            break;
            }

            rshellN+=(prm[prmI].radius-rshellP)*prm[prmI].dev;
            rshellP =prm[prmI].radius;
        }

return rshellN+dr;
}
//---------------------------------------------------------------------------

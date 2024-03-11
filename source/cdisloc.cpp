/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * cdisloc.cpp
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
#include "cdisloc.h"
#include "crandom.h"
#include "affinemat.h"
#include "colormsg.h"

#ifndef __linux
#define M_PI 3.1415926539
#endif

//#define DEBUG


#ifdef DEBUG
#define DB true
#else
#define DB false
#endif



//-----------------------------------------------------------------------------
Cdisloc::Cdisloc()
{

}

//-----------------------------------------------------------------------------
int Cdisloc::findAtomName(cstring &aname__)
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
void Cdisloc::addAtomName(cstring &aname__)
{
const bool anamedup= (findAtomName(aname__)>=0);

        if(anamedup)
            cerr<<"WARNING: atom name duplicate has been found: "<<aname__<<endl;

        if(Elements::mass.find(aname__)==Elements::mass.end())
            cerr<<"WARNING: unknown mass of an atom: "<<aname__<<endl;

        atomTypes.emplace_back(NanoGrain::StAtomType(aname__));
}

//-----------------------------------------------------------------------------
void Cdisloc::clearData()
{
    atomTypes.clear();
    axis.clear();
    axispos.clear();
    rangeR.clear();
    rangeA.clear();
    rangeRoll.clear();
    mindist.clear();
    scatter.clear();
    mode.clear();
    vrpy.clear();

    fileNameIn.clear();
    fileNameOut.clear();

    vrpy.reserve(3);

}

//-----------------------------------------------------------------------------


bool Cdisloc::parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__)
{
const str send("end");

           clearData();
           ptr_uvar=uvar__;

           while(cmd[index]!=send){


               if(DB) cout<<cmd[index]<<endl;

               if(cmd[index]=="angle"){
                   angle=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,angle);
               continue;
               }

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

               if(cmd[index]=="axis"){
                   axis=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,axis);
               continue;
               }

               if(cmd[index]=="mindist"){
                   //axis=cmd[index][1];
                   //Script::replaceVars(ptr_uvar,fileNameIn);
                   mindist=cmd[index++][1];
               continue;
               }

               if(cmd[index]=="mode"){
                   mode=cmd[index++][1];
               continue;
               }


               if(cmd[index]=="position"){
                   axispos=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,axispos);
               continue;
               }

               if(cmd[index]=="projh"){
                   projh=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,projh);
               continue;
               }


               if(cmd[index]=="rangeA"){
                   rangeA=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,rangeA);
               continue;
               }

               if(cmd[index]=="rangeRoll"){
                   rangeRoll=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,rangeA);
               continue;
               }


               if(cmd[index]=="roll" || cmd[index]=="pitch" || cmd[index]=="yaw"){
               std::string keyvalue{cmd[index][0]+" "+cmd[index][1]};
                    Script::replaceVars(ptr_uvar,keyvalue);
                    vrpy.push_back(keyvalue);

                    index++;
                continue;
               }


               if(cmd[index]=="rangeR"){
                   rangeR=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,rangeR);
               continue;
               }


               if(cmd[index]=="save"){
               std::string fileName =cmd[index++][1];
                   Script::replaceVars(ptr_uvar,fileName);
                   fileNameOut.emplace_back(fileName);
               continue;
               }


               if(cmd[index]=="saveopt"){
                   saveopt=cmd[index][1];
                   Script::replaceVars(ptr_uvar,saveopt);
                   index++;
               continue;
               }


               if(cmd[index]=="scatter"){
                   scatter=cmd[index++][1];
                   Script::replaceVars(ptr_uvar,scatter);
               continue;
               }


               cerr<<"Error: unknown command"<<__FILE__<< cmd[index]<<" line "<<__LINE__<<endl;
           return false;
           }



return true;
}
//-----------------------------------------------------------------------------
void Cdisloc::calc()
{
        if(mode=="loop")
            insertLoop();
        else{
            if(mode=="rot"){
                if(axis.empty() || axispos.empty() || rangeR.empty() || angle.empty() || projh.empty()){
                    errMsg(" some instructions not defined");
                throw ERR_NOFPARAMS;
                }
                rotateLoop();
            }            
            if(mode=="rpy"){
                if(axis.empty() || axispos.empty() || rangeR.empty() ){
                    errMsg(" some instructions not defined");
                throw ERR_NOFPARAMS;
                }
                rpyLoop();
            }
            if(mode=="cyl"){
                if(axis.empty() || axispos.empty() || rangeR.empty() ){
                    errMsg(" some instructions not defined");
                throw ERR_NOFPARAMS;
                }
                cylLoop();
            }
        }



        if(!saveopt.empty()){
        auto otoks{split<string>(saveopt,"if")};
        auto ptoks{split<string>(otoks[0],"=")};

        const int  lop=std::stoi(ptoks[0]);
        const int  rop=std::stoi(ptoks[1]);

            if (lop==rop){
                if(!this->fileNameOut.empty()){
                    grain->fileNameOut.clear();
                    grain->fileNameOut.push_back(this->fileNameOut[0]);
                    //std::copy(grain->fileNameOut.begin(),this->fileNameOut.begin(),this->fileNameOut.end());
                    grain->saveToFile();
                }
            }

        }
        else{
            if(!this->fileNameOut.empty()){
                grain->fileNameOut.clear();
                grain->fileNameOut.push_back(this->fileNameOut[0]);
                //std::copy(grain->fileNameOut.begin(),this->fileNameOut.begin(),this->fileNameOut.end());
                grain->saveToFile();
            }
        }

}

//-----------------------------------------------------------------------------


void Cdisloc::insertLoop()
{
//axis params
vector<string> taxis{split<string>(axis," ")};
cpos  A{std::stod(taxis[0])};
cpos  B{std::stod(taxis[1])};
cpos  C{std::stod(taxis[2])};


//position params
vector<string> tpos{split<string>(axispos," ")};
cpos  px{std::stod(tpos[0])};
cpos  py{std::stod(tpos[1])};
cpos  pz{std::stod(tpos[2])};


//radius ranging
vector<string> tradius{split<string>(rangeR," ")};
            if(tradius.size()!=3) {errMsg("rangeR: wrong number of params"); throw Edisstatus::ERR_NOFPARAMS;}

cpos  rmin {std::stod(tradius[0])};
cpos  rstep{std::stod(tradius[1])};
cpos  rmax {std::stod(tradius[2])};


//angle ranging
vector<string> tangle{split<string>(rangeA," ")};
cpos  amin {std::stod(tangle[0])};
cpos  astep{std::stod(tangle[1])};
cpos  amax {std::stod(tangle[2])};


            //
            // main direction is assumed along Z-axis [0,0,1]
            // so:
            // rotation axis =   [0,0,1] x [A,B,C] = [-B,A,0]
            //

StAxis axis(-B,A,0);
cpos mian=std::sqrt(A*A+B*B+C*C);
cpos sa=std::sqrt(A*A+B*B)/mian;
cpos ca=C/mian;
StRotationMatrix rotMat(axis,sa,ca);

string aname=atomTypes[0].name;
int atype;
StVector v;
cpos theta2Rad=M_PI/180.0;
position Ang;
NanoGrain::vatoms dislocAtoms;



        if( (atype=grain->findAtomName(aname)) <0){
            grain->atomTypes.push_back(NanoGrain::StAtomType(aname));
            atype=grain->atomTypes.size()-1;
        }


typedef const size_t csize;
csize nOfrsteps=static_cast<size_t>( (rmax-rmin)/rstep+1 );
csize nOfasteps=static_cast<size_t>( (amax-amin)/astep+1 );
csize disSize=nOfrsteps*nOfasteps;

        dislocAtoms.reserve(disSize);


        for (position r=rmin;r<rmax;r+=rstep){
            for(position th=amin;th<amax;th+=astep){

                Ang=theta2Rad*th;
                v.x=px+r*std::cos(Ang);
                v.y=py+r*std::sin(Ang);
                v.z=pz;

                v=rotMat*v;

                dislocAtoms.emplace_back(NanoGrain::StAtom(v.x,v.y,v.z,atype));
            }
        }


        if(!scatter.empty()){
        std::default_random_engine generator (std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<double> udistr(-1,1);

        //scatter params
        vector<string> tscatt{split<string>(scatter," ")};
        cpos  tx {std::stod(tscatt[1])};
        cpos  ty {std::stod(tscatt[2])};
        cpos  tz {std::stod(tscatt[3])};

                //std::sqrt(rdistr(generatorRadius));

                for(auto & atom: dislocAtoms){
                    atom.x+=tx*udistr(generator);
                    atom.y+=ty*udistr(generator);
                    atom.z+=tz*udistr(generator);
                }


        }


        if(!mindist.empty())
            validateDisloc(dislocAtoms);

        grain->atoms.insert(grain->atoms.begin(),dislocAtoms.begin(),dislocAtoms.end());

}
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
void Cdisloc::rotateLoop()
{
//axis params
vector<string> taxis{split<string>(axis," ")};
cpos  A{std::stod(taxis[0])};
cpos  B{std::stod(taxis[1])};
cpos  C{std::stod(taxis[2])};


//position params
vector<string> tpos{split<string>(axispos," ")};
cpos  px{std::stod(tpos[0])};
cpos  py{std::stod(tpos[1])};
cpos  pz{std::stod(tpos[2])};


//radius ranging
vector<string> tradius{split<string>(rangeR," ")};
            if(tradius.size()!=2) {errMsg("rangeR: wrong number of params (should be 2)"); throw Edisstatus::ERR_NOFPARAMS;}

cpos  rmin {std::stod(tradius[0])};
cpos  rmax {std::stod(tradius[1])};

cpos angle_=std::stod(angle);
cpos projh_=std::stod(projh);

/*
int atype=grain->atomTypes.size()-1;

                if(  !atomTypes.empty()){
                const string aname=atomTypes[0].name;

                    if( grain->findAtomName(aname) < 0 ){
                        grain->atomTypes.push_back(NanoGrain::StAtomType(aname));
                        atype=grain->atomTypes.size()-1;
                    }
                }*/

StVector point;
const StAxis axis(A,B,C,px,py,pz);
position d,ph;

cpos sa=std::sin(angle_*M_PI/180);
cpos ca=std::cos(angle_*M_PI/180);
StRotationMatrix rotMat(axis,sa,ca);

#ifdef DB
size_t nOfrot=0;
#endif

                if(scatter.empty()){
                    for(auto &atom: grain->atoms){
                        point.x=atom.x;
                        point.y=atom.y;
                        point.z=atom.z;

                        d=pointPlaneDistance(axis,point);

                        if(rmin<d && d<rmax){
                            ph=std::fabs(projLength(axis,point));
                            if(ph<projh_){
                                point.x-=axis.xo;
                                point.y-=axis.yo;
                                point.z-=axis.zo;
                                point=rotMat*point;
                                atom.x=point.x+axis.xo;
                                atom.y=point.y+axis.yo;
                                atom.z=point.z+axis.zo;
                                //atom.atype=atype;

                                #ifdef DB
                                nOfrot++;
                                #endif
                            }
                        }
                    }
                }
                else{
                std::default_random_engine generator (std::chrono::system_clock::now().time_since_epoch().count());
                std::uniform_real_distribution<double> udistr(-1,1);

                //scatter params
                vector<string> tscatt{split<string>(scatter," ")};
                cpos  tx {std::stod(tscatt[1])};
                cpos  ty {std::stod(tscatt[2])};
                cpos  tz {std::stod(tscatt[3])};

                        //std::sqrt(rdistr(generatorRadius));
                        for(auto &atom: grain->atoms){
                            point.x=atom.x;
                            point.y=atom.y;
                            point.z=atom.z;

                            d=pointPlaneDistance(axis,point);

                            if(rmin<d && d<rmax){
                                ph=std::fabs(projLength(axis,point));
                                if(ph<projh_){
                                    point.x-=axis.xo;
                                    point.y-=axis.yo;
                                    point.z-=axis.zo;
                                    point=rotMat*point;
                                    atom.x=point.x+axis.xo+tx*udistr(generator);
                                    atom.y=point.y+axis.yo+ty*udistr(generator);
                                    atom.z=point.z+axis.zo+tz*udistr(generator);

                                    #ifdef DB
                                    nOfrot++;
                                    #endif

                                    //atom.atype=atype;
                                }
                            }
                        }
                }


                #ifdef DB
                infoMsg("number of rotated atoms: "+std::to_string(nOfrot));
                #endif


}
//-----------------------------------------------------------------------------
StAxis froll(const StVector &b, const StVector &c)
{
return StAxis(crossProduct(b,c));
}

StAxis fpitch(const StVector &b, const StVector &c)
{
return StAxis(c);
}

StAxis fyaw(const StVector &b, const StVector &c)
{
return StAxis(b);
}
///         axes/vectos orientation, all vectors are coplanar  vec_c=vec_bx(vec_a x vec_b) ; vec_b is defined by user
///
///             ^ vec_b
///             |
///             |
///             |    _ o _
///             |    /| |\
///             |   /     \ vec_dr
///             |  /vec_a  \
///             | /         \
///             |------------*------>vec_c
///
///
///
///
//-----------------------------------------------------------------------------
class acceptAtom
{
public:
    virtual ~acceptAtom() {  }
    virtual bool accept(StVector & r, StVector &a, StVector &b)=0;

    position rmin, rmax;



};


class testRadii: public acceptAtom
{
public:
    bool accept(StVector & r, StVector &a, StVector &b ){
        return r.getModule()<rmax;
    }


};


class testRadiiAngle: public acceptAtom
{

public:
    position angleMin,angleMax,pn=1;


    bool accept(StVector & r, StVector &a, StVector &b){

            if (r.getModule()>rmax) return false;

    cpos  ca{cosa(r,a)};
    cpos  cb{pn*cosa(r,b)};
    cpos angle{std::atan2(ca,cb)*180/M_PI};

    return ( angleMin<angle ) && (angle<angleMax);
    //return true;
    }

};

//-----------------------------------------------------------------------------
//SVector sumVec vec_dr=rotMat*vec_dr;

void rotAtomNoRnd(StRotationMatrix &rm, StVector &a, StVector &rndShift)
{
        a=rm*a;
}

void rotAtomRndShift(StRotationMatrix &rm, StVector &a, StVector &rndShift)
{
static std::default_random_engine generator (std::chrono::system_clock::now().time_since_epoch().count());
static std::uniform_real_distribution<double> udistr(-1,1);

        a=rm*a;
        a.x+=rndShift.x*udistr(generator);
        a.y+=rndShift.y*udistr(generator);
        a.z+=rndShift.z*udistr(generator);
}


//-----------------------------------------------------------------------------
void Cdisloc::rpyLoop()
{
//axis params
vector<string> taxis{split<string>(axis," ")};
cpos  A{std::stod(taxis[0])};
cpos  B{std::stod(taxis[1])};
cpos  C{std::stod(taxis[2])};

//position params
vector<string> tpos{split<string>(axispos," ")};
cpos  px{std::stod(tpos[0])};
cpos  py{std::stod(tpos[1])};
cpos  pz{std::stod(tpos[2])};
StVector axisPos(px,py,pz);

//radius ranging
vector<string> tradius{split<string>(rangeR," ")};
            if(tradius.size()!=2) {errMsg("rangeR: wrong number of params (should be 2)"); throw Edisstatus::ERR_NOFPARAMS;}

cpos  rmin {std::stod(tradius[0])};
cpos  rmax {std::stod(tradius[1])};
cpos  rave {0.5*(rmin+rmax)};
cpos  dr   {rmax-rave};


//
StVector vec_a,vec_c,vec_dr;
StVector vec_b(A,B,C);
StVector point;
StAxis (*frpy)(const StVector &, const StVector &);
position imodC,modA;
position cos_vavb;

acceptAtom *testAtom;

                ///////////////////////////////////////

                if(rangeRoll.empty())
                    testAtom=new testRadii();
                else{
                    testAtom=new testRadiiAngle();

                testRadiiAngle *tra{dynamic_cast<testRadiiAngle *>(testAtom)};
                vector<string> ratoks{split<string>(rangeRoll," \t")};

                    tra->angleMin=std::stod(ratoks[0]);
                    tra->angleMax=std::stod(ratoks[1]);

                    if(tra->angleMin<-180 || tra->angleMax>180){
                        errMsg(" wrong values of  roll angle range");
                        delete testAtom;
                    throw Edisstatus::ERR_AROLLRANGE;
                    }

                    if(ratoks.size()==3){ //pn option
                        tra->pn=(ratoks[2]=="p") ? 1 : -1;
                    }
                }


                testAtom->rmax=dr;

                ///////////////////////////////////////

                if(  !atomTypes.empty()){
                const string aname=atomTypes[0].name;
                    if( grain->findAtomName(aname) < 0 )
                        grain->atomTypes.push_back(NanoGrain::StAtomType(aname));
                }


const size_t atype=grain->atomTypes.size()-1;

                ///////////////////////////////////////

void  (*fRotRndShift)(StRotationMatrix &rm, StVector &a, StVector &rndShift);
StVector rndShift;

                if(scatter.empty()){
                    fRotRndShift=&rotAtomNoRnd;
                    rndShift.x=rndShift.y=rndShift.z=0;
                }
                else{
                vector<string> tscatt{split<string>(scatter," ")};
                cpos  tx {std::stod(tscatt[1])};
                cpos  ty {std::stod(tscatt[2])};
                cpos  tz {std::stod(tscatt[3])};

                    fRotRndShift=&rotAtomRndShift;
                    rndShift=StVector(tx,ty,tz);
                }


                ///////////////////////////////////////
                ///
                ///

position ca,cb;

                for(auto &rpy: vrpy){
                vector<string> rpyToks(split<string>(rpy," "));
                cpos sa=std::sin(std::stod(rpyToks[1])*M_PI/180);
                cpos ca=std::cos(std::stod(rpyToks[1])*M_PI/180);

                        if(rpyToks[0]=="roll")
                            frpy=&froll;
                        else{
                            if(rpyToks[0]=="pitch")
                                frpy=&fpitch;
                            else
                                frpy=&fyaw;
                        }

                        for(auto &atom: grain->atoms){                                
                                vec_a=atom.Pos()-axisPos;
                                cos_vavb=std::fabs(cosa(vec_a,vec_b));
                                modA=vec_a.getModule();

                                if(  modA>rmax || modA<rmin ||
                                     cos_vavb>1-1e-6) { continue; }


                                vec_c=crossProductTriple(vec_b,vec_a,vec_b);
                                imodC=rave/vec_c.getModule();

                                vec_c*=imodC;
                                vec_dr=vec_a-vec_c;

                                //ca=cosa(vec_dr,vec_c);
                                //cb=cosa(vec_dr,vec_b);

                                if( testAtom->accept(vec_dr,vec_b,vec_c) ){
                                //if(vec_dr.getModule()<dr  ){
                                StRotationMatrix rotMat(frpy(vec_b,vec_c),sa,ca);

                                        //vec_dr=rotMat*vec_dr;
                                        fRotRndShift(rotMat,vec_dr,rndShift);

                                        atom.x=vec_c.x+vec_dr.x+px;
                                        atom.y=vec_c.y+vec_dr.y+py;
                                        atom.z=vec_c.z+vec_dr.z+pz;

                                        atom.atype=atype;
                                }
                        }
                }///end for


                delete testAtom;
}
//-----------------------------------------------------------------------------
void Cdisloc::cylLoop()
{
//axis params
vector<string> taxis{split<string>(axis," ")};
cpos  A{std::stod(taxis[0])};
cpos  B{std::stod(taxis[1])};
cpos  C{std::stod(taxis[2])};

//position params
vector<string> tpos{split<string>(axispos," ")};
cpos  px{std::stod(tpos[0])};
cpos  py{std::stod(tpos[1])};
cpos  pz{std::stod(tpos[2])};
StVector axisPos(px,py,pz);

//radius ranging
vector<string> tradius{split<string>(rangeR," ")};
            if(tradius.size()!=2) {errMsg("rangeR: wrong number of params (should be 2)"); throw Edisstatus::ERR_NOFPARAMS;}

cpos  rmin {std::stod(tradius[0])};
cpos  rmax {std::stod(tradius[1])};
cpos  rave {0.5*(rmin+rmax)};
cpos  dr   {rmax-rave};

//
StVector vec_a,vec_c,vec_dr,vec_r;
StVector vec_b(A,B,C);

StAxis mainAxis(A,B,C,px,py,pz);
StAxis (*frpy)(const StVector &, const StVector &);

//cpos angle_=std::stod(angle);
cpos projh_=std::stod(projh);


//angle ranging
vector<string> tangle{split<string>(rangeA," ")};
cpos  amin {std::stod(tangle[0])};
//cpos  astep{std::stod(tangle[1])};
cpos  amax {std::stod(tangle[1])};

cpos ao=(amax-amin)/(2*projh_);
cpos bo=ao*2*projh_-amin;


void  (*fRotRndShift)(StRotationMatrix &rm, StVector &a, StVector &rndShift);
StVector rndShift;

                if(scatter.empty()){
                    fRotRndShift=&rotAtomNoRnd;
                    rndShift.x=rndShift.y=rndShift.z=0;
                }
                else{
                vector<string> tscatt{split<string>(scatter," ")};
                cpos  tx {std::stod(tscatt[1])};
                cpos  ty {std::stod(tscatt[2])};
                cpos  tz {std::stod(tscatt[3])};

                    fRotRndShift=&rotAtomRndShift;
                    rndShift=StVector(tx,ty,tz);
                }

                ///////////////////////////////////////

                if(  !atomTypes.empty()){
                const string aname=atomTypes[0].name;
                    if( grain->findAtomName(aname) < 0 )
                        grain->atomTypes.push_back(NanoGrain::StAtomType(aname));
                }


const size_t atype=grain->atomTypes.size()-1;

                ///////////////////////////////////////
position d,ph,alpha,sa,ca;
cpos deg2rad=M_PI/180;


                for(auto &atom: grain->atoms){                
                    vec_r=atom.Pos()-axisPos;
                    d=pointAxisDistance(mainAxis,vec_r);

                    if(rmin<d && d<rmax){                                                
                        ph=projLength(mainAxis,atom.Pos());

                        if(std::fabs(ph)<projh_){
                        alpha=deg2rad*(ao*ph+bo);
                        sa=std::sin(alpha);
                        ca=std::cos(alpha);
                        StRotationMatrix rotMat(mainAxis,sa,ca);

                            fRotRndShift(rotMat,vec_r,rndShift);

                            atom.x=vec_r.x+px;
                            atom.y=vec_r.y+py;
                            atom.z=vec_r.z+pz;

                            atom.atype=atype;

                        }
                    }
                }


}

//-----------------------------------------------------------------------------

inline double sqrd (const position &x)  {return x*x;}


void Cdisloc::validateDisloc(NanoGrain::vatoms &dislocAtoms)
{
auto dtoks{split<string>(mindist," ")} ;
cpos minDist2=sqrd(std::stod(dtoks[0]));
NanoGrain::vatoms tmpAtoms(std::move(dislocAtoms));
const size_t nOfdislocAtoms=tmpAtoms.size();
const size_t nOfgrainAtoms=grain->atoms.size();
double r2;
auto calcR2=[](NanoGrain::StAtom *A, NanoGrain::StAtom *B){
        return sqrd(A->x-B->x) + sqrd(A->y-B->y) + sqrd(A->z-B->z) ; };


            dislocAtoms.clear();
            dislocAtoms.reserve(tmpAtoms.size());


NanoGrain::StAtom *ptrA=tmpAtoms.data();
NanoGrain::StAtom *ptrB;
bool valid;

             for(size_t i=0;i<nOfdislocAtoms;i++,ptrA++){
                valid=true;
                ptrB=grain->atoms.data();

                for(size_t j=0;j<nOfgrainAtoms;j++,ptrB++){
                    r2=calcR2(ptrA,ptrB);
                    if(r2<minDist2){
                        valid=false;
                        break;
                    }
                }

                if(valid)
                    dislocAtoms.emplace_back(*ptrA);

             }
}



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
    mindist.clear();
    scatter.clear();
    mode="loop";

    fileNameIn.clear();
    fileNameOut.clear();

}




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
            if(mode=="spin"){
                spinLoop();

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
StRotationMatrix rotMat(axis,sa);

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
StRotationMatrix rotMat(axis,sa);

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
                            ph=projHeight(axis,point);
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
                                ph=projHeight(axis,point);
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
void Cdisloc::spinLoop()
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
cpos  rave {0.5*(rmin+rmax)};
cpos  dr   {rmax-rave};

cpos angle_{std::stod(angle)};
//cpos projh_{std::stod(projh)};

StVector vec_a;
//const StAxis axis(A,B,C,px,py,pz);
StVector vec_b(A,B,C);
StVector point;
position d,ph;

cpos sa=std::sin(angle_*M_PI/180);
//StRotationMatrix rotMat(axis,sa);


int atype=grain->atomTypes.size()-1;

                if(  !atomTypes.empty()){
                const string aname=atomTypes[0].name;

                    if( grain->findAtomName(aname) < 0 ){
                        grain->atomTypes.push_back(NanoGrain::StAtomType(aname));
                        atype=grain->atomTypes.size()-1;
                    }
                }


//int ia=0;
//std::default_random_engine generator (std::chrono::system_clock::now().time_since_epoch().count());
//std::uniform_real_distribution<double> udistr(-0.1,0.1);


                for(auto &atom: grain->atoms){
                        //ia++;
                StVector vec_a(atom.x-px,atom.y-py,atom.z-pz);
                StVector vec_c{crossProductTriple(vec_b,vec_a,vec_b)};
                cpos imodC=rave/vec_c.getModule();
                          vec_c*=imodC;

                StVector vec_dr=vec_a-vec_c;

                        if(vec_dr.getModule()<dr){

                        StAxis spinAxis( crossProduct(vec_b,vec_c));
                        //position sa=udistr(generator);
                        StRotationMatrix rotMat(spinAxis,sa);

                                //if(ia>160 && ia<170) {
                                 //   rotMat.showMatrix(); cout<<endl;}

                                vec_dr=rotMat*vec_dr;
                                atom.x=vec_c.x+vec_dr.x+px;
                                atom.y=vec_c.y+vec_dr.y+py;
                                atom.z=vec_c.z+vec_dr.z+pz;
                                atom.atype=1;

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



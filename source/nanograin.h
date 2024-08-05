
/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * nanograins.h
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



#ifndef NANOGRAIN_H
#define NANOGRAIN_H


#include<vector>
#include<iostream>
#include<string>
#include<fstream>
#include<cmath>
#include<sstream>
#include<list>

#include "scriptanalyser.h"
#include "elements.h"
#include "affinemat.h"

#include <random>
#include <chrono>

using namespace std;


//
/// Susumu Onaka
/// Nanomaterials 2016, 6, 27; doi:10.3390/nano6020027
//


//-----------------------------------------------------------------------------
typedef double position;
typedef const position cpos;
typedef const double cdouble;
typedef std::string str;




//----------------------------------------------------------------------------------------------------------------------------------------------------------
//
// CSuperSphere class types definitions are placed in:  csupersphere.cpp
//
//----------------------------------------------------------------------------------------------------------------------------------------------------------

class CSuperSphere
{
protected:
        double  Rp;
        string info;
        const std::string name;

public:
        const double p;
        const double R;
        double xa,yb,zc;

        CSuperSphere():name("UNK"),p(0),R(1) { }
        CSuperSphere(const std::string name__,const double p__, const double R__)
            :name(name__),p(p__),R(R__),xa(1),yb(1),zc(1)
        {
            Rp=std::pow(R,p);  info=name+" (R,p)=("+std::to_string(R)+","+std::to_string(p)+")";
        }

        CSuperSphere(const std::string name__,cdouble p__, cdouble R__, cdouble xa__,cdouble yb__, cdouble zc__ )
            :name(name__),p(p__),R(R__),xa(xa__),yb(yb__),zc(zc__)
        {
            Rp=std::pow(R,p);  info=name+" (R,p,a,b,c)=("+std::to_string(R)+","+
                                                        std::to_string(p)+","+
                                                        std::to_string(xa)+","+std::to_string(yb)+","+std::to_string(zc)+")";
        }

        virtual ~CSuperSphere() { }




        virtual bool isValid(const position &x, const position &y, const position &z)=0;
        const string & getInfo(){ return info;}

};

//................................................
class CCubic : public CSuperSphere
{

public:
        CCubic(const double p__, const double R__) : CSuperSphere("cubic structure",p__,R__) { }
        CCubic(const double p__, const double R__,cdouble xa__, cdouble yb__,cdouble zc__)
            : CSuperSphere("cubic structure",p__,R__,xa__,yb__,zc__) { }

        bool isValid(const position &x, const position &y, const position &z);

};

//................................................
class COctahedral: public CSuperSphere
{
public:
        COctahedral(const double p__, const double R__) : CSuperSphere("octahedral structure",p__,R__) { }
        COctahedral(const double p__, const double R__,cdouble xa__, cdouble yb__,cdouble zc__)
            : CSuperSphere("octahedral structure",p__,R__,xa__,yb__,zc__) { }

        bool isValid(const position &x, const position &y, const position &z);

};
//................................................
class CDodecahedral: public CSuperSphere
{
public:
        CDodecahedral(const double p__, const double R__) : CSuperSphere("octahedral structure",p__,R__) { }
        CDodecahedral(const double p__, const double R__,cdouble xa__, cdouble yb__,cdouble zc__)
            : CSuperSphere("octahedral structure",p__,R__,xa__,yb__,zc__) { }

        bool isValid(const position &x, const position &y, const position &z);
};

//................................................
class CPolyhedral: public CSuperSphere
{
private:
      const   double a,b;
      double iap,ibp;
public:
        CPolyhedral(const double p__, const double R__) :  CSuperSphere("polyhedral structure",p__,R__),a(1.0/0.0),b(1.0/0.0)
        {
            info+=", (a,b)=("+std::to_string(a)+","+std::to_string(b)+")";
        }
        CPolyhedral(const double p__, const double R__,const double a__,const double b__) ;
        CPolyhedral(const double p__, const double R__,const double a__,const double b__,cdouble xa__, cdouble yb__,cdouble zc__) ;

        bool isValid(const position &x, const position &y, const position &z);
};

//................................................
class CPolyhedral2D: public CSuperSphere
{
public:

    CPolyhedral2D(const double p__, const double R__):CSuperSphere("2D superellipse", p__,R__) { }

    bool isValid(const position &x, const position &y, const position &z);
} ;
//................................................
class CPolyhedral2D_HOD: public CSuperSphere
{
private:
      const   double a,b;
      double iap,ibp;
public:

    CPolyhedral2D_HOD(const double p__, const double R__,const double a__,const double b__);

    bool isValid(const position &x, const position &y, const position &z);
} ;

//-----------------------------------------------------------------------------
namespace NanoGrain{

enum Status{OK,ERR_FOPEN,ERR_FIN,ERR_FIN_EXC,ERR_FFORMAT,ERR_LP,ERR_RADII,ERR_GEOM,ERR_ATYPES,ERR_ANUM_ZERO,ERR_UNKSEQ,ERR_UNKSEQNUM,ERR_WRSHAPE,
                WARN_ANAMEDUP, WARN_ANAMEUNKMASS,ERR_UNKAN,ERR_LPLTZERO,ERR_RLTZERO,
                ERR_ABCPARSE,ERR_UNKSHAPE,ERR_HEIGHT,ERR_ST};




struct StAtom{
position x,y,z,r2;
static size_t idIterator;
size_t id;
size_t atype;  ///  position of an atom type in array
bool mark=false;
bool rdh=false;
vector<size_t> neighID;  /// ID of neigneighborhood atoms

        StAtom(){x=y=z=r2=0;setId();}
        StAtom(cpos &x__, cpos &y__, cpos &z__){x=x__;y=y__;z=z__;calcR2();setId();}
        StAtom(cpos &x__, cpos &y__, cpos &z__, const size_t &atype__){x=x__;y=y__;z=z__;calcR2();atype=atype__;setId();}
        StAtom(cpos &x__, cpos &y__, cpos &z__, const size_t &atype__,cpos &r2__){x=x__;y=y__;z=z__;r2=r2__;atype=atype__;setId();}
        //StAtom(cpos &x__, cpos &y__, cpos &z__,str * const aname){x=x__;y=y__;z=z__;p_name=aname;calcR2();}
        //StAtom(cpos &x__, cpos &y__, cpos &z__,str * const aname,cpos &r2__){x=x__;y=y__;z=z__;p_name=aname;r2=r2__;}
        void calcR2(){r2=x*x+y*y+z*z;}
        void rescale(const position &p){ x*=p;y*=p;z*=p;calcR2();}
        void rescale(const position &px, const position &py, const position &pz){ x*=px;y*=py;z*=pz;calcR2();}

        void operator*=(const position &p)
        {
            x*=p;y*=p;z*=p;
            r2*=p*p;
        }

        void operator+=(const StVector &v)
        {
            x+=v.x;
            y+=v.y;
            z+=v.z;
        }

        /*
        StVector operator+(const StVector &v)
        {
        return StVector(x+v.x,y+v.y,z+v.z);
        }*/

        StVector Pos()
        {
        return StVector(x,y,z);
        }

        void operator=(const StVector &v)
        {
            x=v.x;y=v.y;z=v.z;
        }

        void setId(){ id=idIterator; idIterator++; }

        static void resetId(){idIterator=0;}

        friend ostream & operator<<(ostream &,const NanoGrain::StAtom &a);
};


ostream & operator<<(ostream &,const NanoGrain::StAtom &a);
typedef vector<StAtom> vatoms;

//---------------------------------------------------------------------------

struct StAtomType{
std::string name,charge;

    StAtomType() { name ="?"; charge="?";  }
    StAtomType(const string name__):name (name__) {  }
    ~StAtomType(){  }

   // bool operator()(const string &a){ return a==name;}
    bool operator == (const StAtomType &a){
        return a.name==name;
    }

};



//---------------------------------------------------------------------------
struct StMinMax{
    union{
        struct {position xmin,xmax,ymin,ymax,zmin,zmax;};
        position xyzMinMax[6];
    };


    void searchForMinMax(const vatoms & atoms);
    position getMaxAbs();
    position getMinAbs();
    position getMinPos();
    position getWidth() {return xmax-xmin;}
    position getLength(){return ymax-ymin;}
    position getHeight(){return zmax-zmin;}


};
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
class StNanoGrain{

private:
    string numOfAtomsPush,saveFileStatusPush;

class CShells{
public:

    struct shell{
    position dev,radius;

        shell() {}
        shell(const string &dev__, const string &radius__)
        {
            dev =1+std::stod(dev__)/100.0;
            radius=std::stod(radius__);
        }

        shell(const double dev__, const double radius__)
        {
            dev=1+dev__/100.0;
            radius=radius__;
        }

    };


    vector<shell> prm;
    shell & operator [] (size_t i) {return prm[i];}
    position newR(const position &r);
} coreshell;


enum ESOURCE {AUTOCREATE,FILE} source;
enum EFCCTYPE{FCC,ZB,UO2};
enum EFAULTMODE{OFF,RANDOM,CUSTOM} faultmode;
enum EHCPSURF{sA,sB,sAB,sH} hcpsurf;
enum ECENTER {COFF,GEOM,CATOM,ID} center;

vatoms baseAtomsHcp;
position hcpa,hcpc;
position hcpRadius,hcpRadius2,hcpCylinderHeight;
vector<size_t> hcpFaultPos;
vector<string> replicate;
StVector a,b,c,d;
std::string hcpfillup,hcplaynum,hcpfault;
bool dispParams;
stdumap *ptr_uvar;
std::stringstream atomsRemoved;
//std::map<std::string, double > elcharge;


        void sc();
        void bcc();
        void fcc(EFCCTYPE fcctype=FCC);

        void hcp();
        void zb110();

        void feni();
        void buildFromUC();
        void buildFromTric();


        void buildHcpABC();
        void buildRandomABC(vector<string> &params);

        void buildHcpBaseHex();
        void buildHcpBaseRhm();
        void insertFaults();
        void buildHcpLayers();

        void buildHcpCylinder(const StAtom &atom,const StVector &vres);
        void buildHcpSphere  (const StAtom &atom,const StVector &vres);
        void buildHcpSuperSphere(const StAtom &atom,const StVector &vres);
        void (StNanoGrain::*fbuildHcpShape)(const StAtom &atom,const StVector &vres);
        void buildSuperSphere();

        void grainCentering();
        void catomCentering();
        void idCentering();
        void insertDisloc();

        static bool testSavedNumOfAtoms(const size_t );


        void saveDatFile(const string &fileName);
        void saveXYZFile(const string &fileName);
        void saveNXYZFile(const string &fileName);
        void saveMXYZFile(const string &fileName);
        void saveNDLFile(const string &fileName);
        void saveLammpsFile(const string &fileName);
        void saveHeader();


        void openFile();
        void openXYZFile();
        void openLMPFile();
        void openNDLFile();
        void sortAtomsByName();


        void build();
        void voids();
        void disperseN();
        void displaceAtoms();
        void rescale();
        void renameAtoms();
        void removeAtoms();
        void findMaxR();
        void findAtomNamesNumber();
        position getLP();
        position getRadius(cpos & lp);

        void dispParameters();

        struct StAtomDispFactor{
        string name,dispfact;
                StAtomDispFactor() { }
                StAtomDispFactor(const string &name__,const string &df__)
                {name=std::move(name__);dispfact=std::move(df__);}
        };

        /// Unit Cell atoms
        struct StUcAtom{
        position x,y,z;
        string name;size_t id;
        bool rdh=false;
        double rmProb=0;
                StUcAtom() { }
                StUcAtom(const string &name__):name(name__){ }
                bool operator()(const string &a){ return a==name;}
                void operator=(const StVector &v){x=v.x; y=v.y; z=v.z;}

        };

public:
std::string side,radius,clp,structure,scaleFactors;
std::string shape,shapePrm,shapePrm2D;
std::string disloc,dislocPlane,voidsprm;
std::string hcpu,hcpcs, hcpABC;
vector<StAtomType> atomTypes;
vector<size_t> atomNamesNumber;
vector<StAtomDispFactor> atomDisperse;
std::string fileNameIn,fileNameHeader;
vector<string> fileNameOut;
std::string lmpstyle,comment;
std::string rmatoms;
std::string threads,mthreads;
std::string rename;
std::string margins;
std::string catomType;
position maxR,maxX,maxY,maxZ; //maximal distance from (0,0,0)
vatoms atoms;
bool disperse,hcpsl,numOfAtomsTest;//,testSNA;
static list<size_t> savedNumOfAtoms;
CSuperSphere *ssShape=nullptr;
void (StNanoGrain::*callbackSetThreads)(std::string &threads);
    //************************************************
    struct StUnitCell{
        vector<StVector> vtrans;
        size_t rdhAtoms;      
        vector<StUcAtom> atoms;
        void clear() {rdhAtoms=0;vtrans.clear(); atoms.clear();}
        bool empty() {return vtrans.empty();}

    } uc;

    //************************************************
    struct StTric{
    string lpa,lpb,lpc;
    string alpha,beta,gamma;
    position xy, xz, yz;
    size_t rdhAtoms;
    vector<StUcAtom> atoms;
    enum ECOORDSYS{TRIC,CART} csys;   // coordination system

            StTric(){reset(); }

            void reset(){rdhAtoms=0;alpha=beta=gamma=std::string("90");  csys=ECOORDSYS::TRIC; lpa.clear(); lpb.clear(); lpc.clear(); }
            bool empty(){return lpa.empty() || lpb.empty() || lpc.empty(); }
    } tric;

    //************************************************
    struct StSaveOpt{
        size_t min,max;
        vector<string> lwh;
        bool   fileSaved;
        bool   lmpTric;

        StSaveOpt(){ }
        void reset(){ min=0;max=-1;  lwh.clear(); fileSaved=false;lmpTric=false;}

    } saveopt;

    //************************************************


    void resetPrms();
    bool parseCommands(vcmdlist &cmd, size_t &index, stdumap *uvar__);
    int findAtomName(const string &aname__);
    void addAtomName(const string &aname__);
    void saveToFile();

    friend ostream & operator<< (ostream &, StNanoGrain &grain) ;
};
//************************************************





ostream & operator<< (ostream &, StNanoGrain &grain) ;

} //end namespace :: Grain




//-----------------------------------------------------------------------------
#endif // NANOGRAIN_H

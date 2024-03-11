/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * affinemat.h
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
 
#ifndef AFFINEMAT_H
#define AFFINEMAT_H

#include <cmath>

typedef double position;
typedef const position cpos;

//-----------------------------------------------------------------------------

struct StVector{
    union{
    position x, _1;};

    union{
    position y, _2;};

    union{
    position z, _3;};


      StVector operator+(const StVector & op) const
      {
       StVector res;

                res.x=this->x+op.x;
                res.y=this->y+op.y;
                res.z=this->z+op.z;

       return res;
      }

      StVector operator-(const StVector & op) const
      {
       StVector res;

                res.x=this->x-op.x;
                res.y=this->y-op.y;
                res.z=this->z-op.z;

       return res;
      }

      StVector operator*(const double &val) const
      {
       StVector res;

                res.x=this->x*val;
                res.y=this->y*val;
                res.z=this->z*val;

       return res;
      }

      StVector  operator*=(const double &val)
      {
                x*=val;
                y*=val;
                z*=val;
       return *this;
      }


      //int operator = (const StGrainNode &node);
      //int operator +=(const StGrainNode &node);
      StVector() {  }
      StVector(cpos &x__, cpos &y__, cpos &z__):x(x__),y(y__),z(z__) { }


      double getModule()const { return sqrt(x*x+y*y+z*z);}
      bool isZero() const { return   (x==0) && (y==0) && (z==0); }

};


//---------------------------------------------------------------------------
struct StAxis{
union{
    struct{position a,b,c,xo,yo,zo;};
    position prm[6];
    };
bool on=false;
double tmax,tmin;


    StAxis(cpos a__,cpos b__, cpos c__, cpos xo__, cpos yo__, cpos zo__)
    : a(a__),b(b__),c(c__),xo(xo__),yo(yo__),zo(zo__) {  }

    StAxis(cpos &x, cpos &y, cpos &z){
        a=x;b=y;c=z;
    }

    StAxis(StVector v){
        a=v.x;
        b=v.y;
        c=v.z;
    }

    double getModule() const {return sqrt(a*a+b*b+c*c);}
    StVector getABC()  const {return StVector(a,b,c);}

};


class StRotationMatrix{
public:

union{
      struct{ position m11,m12,m13,m21,m22,m23,m31,m32,m33;};
      position m[9];
      position mm[3][3];
};

position ux,uy,uz;//axis of rotation
position theta;

bool on=false;

    void buildMatrix(const StAxis &axis_);
    void buildMatrix(const StAxis &axis_,cpos &sinA, cpos &cosA);
    StVector operator*(const StVector &v);
//StAtom operator*=(S);

    StRotationMatrix(){ ux=uy=uz=theta=0;}
    StRotationMatrix(const StAxis &axis_,cpos &sinA, cpos &cosA){
        buildMatrix(axis_,sinA,cosA);
    }
    void showMatrix();

private:
    void buildRotationAxis(const StAxis &axis_);
    void normUxyz();
};
//-----------------------------------------------------------------------------


cpos cosa(const StVector &a, const StVector &b);
StVector crossProduct(const StVector &a, const StVector &b);
StVector crossProductTriple(StVector &a, StVector &b, StVector &c);
cpos     tripleProduct(StVector &a, StVector &b, StVector &c);
position projLength(const StAxis &a, const StVector &b);
position pointAxisDistance(const StAxis &axis,const StVector &point);
position pointPlaneDistance(const StAxis &axis,const StVector &point);



#endif // AFFINEMAT_H

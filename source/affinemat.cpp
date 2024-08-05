/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * affinemat.cpp
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
 
#include "affinemat.h"
#include <iostream>
#include <iomanip>

#ifdef DEBUG
#define DB true
#else
#define DB false
#endif


inline cpos sqr(cpos &x){return x*x;}


void StRotationMatrix::buildMatrix(const StAxis &axis_)
{

}

///
///https://en.wikipedia.org/wiki/Rotation_matrix
///

void StRotationMatrix::buildMatrix(const StAxis &axis_, cpos &sinA, cpos &cosA)
{
cpos sinTheta=sinA;//sqrt(1-sqrCoord(cosTheta));
cpos cosTheta=cosA;
cpos oneMinCos=1-cosTheta;

        buildRotationAxis(axis_);

        ///rotation matrix elements
        m11=cosTheta+sqr(ux)*oneMinCos;  m12=ux*uy*oneMinCos-uz*sinTheta; m13=ux*uz*oneMinCos+uy*sinTheta;
        m21=uy*ux*oneMinCos+uz*sinTheta; m22=cosTheta+sqr(uy)*oneMinCos;  m23=uy*uz*oneMinCos-ux*sinTheta;
        m31=uz*ux*oneMinCos-uy*sinTheta; m32=uz*uy*oneMinCos+ux*sinTheta; m33=cosTheta+sqr(uz)*oneMinCos;

        on=true;
}

//-----------------------------------------------------------------------------
StVector StRotationMatrix::operator *(const StVector &a)
{
position bx,by,bz;

                bx=m11*a.x+m12*a.y+m13*a.z;
                by=m21*a.x+m22*a.y+m23*a.z;
                bz=m31*a.x+m32*a.y+m33*a.z;

return StVector(bx,by,bz);
}

void StRotationMatrix::showMatrix()
{
            std::cout<<"[ ";
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++)
                    std::cout<<" "<<std::setw(9)<<mm[i][j];

                if(i<2) std::cout<<std::endl<<"  ";
                else std::cout<<"   ]";
            }
}

//-----------------------------------------------------------------------------
void StRotationMatrix::buildRotationAxis(const StAxis &axis_)
{
        ux=axis_.a;
        uy=axis_.b;
        uz=axis_.c;

        normUxyz();
}

//-----------------------------------------------------------------------------
void StRotationMatrix::normUxyz()
{
cpos sumSq=sqr(ux)+sqr(uy)+sqr(uz);
            if(DB){
                if(sumSq<1e-9){
                    std::cerr<<__FILE__<<":"<<__LINE__<<"   ERROR: axis length <1e-9 "<<std::endl;
                    ux=1;
                    uy=uz=0;
                    return;
                }
            }

cpos norm=sqrt(1/sumSq);

            ux*=norm;
            uy*=norm;
            uz*=norm;
}
//-----------------------------------------------------------------------------
StVector crossProduct(const StVector &a, const StVector &b)
{
cpos cx=(a.y*b.z-a.z*b.y);
cpos cy=(a.z*b.x-a.x*b.z);
cpos cz=(a.x*b.y-a.y*b.x);

return StVector(cx,cy,cz);
}
//-----------------------------------------------------------------------------
position projLength(const StAxis &axis, const StVector &b)
{
cpos aDotb=(axis.a*(axis.xo-b.x)+axis.b*(axis.yo-b.y)+axis.c*(axis.zo-b.z));
return aDotb/axis.getModule();
}
//-----------------------------------------------------------------------------
position pointAxisDistance(const StAxis &axis, const StVector &point)
{
//StVector r{point.x-axis.xo, point.y-axis.yo, point.z-axis.zo};
/// h= |AxR|/|A|
///

return crossProduct(axis.getABC(),point).getModule()/axis.getModule();
}
//-----------------------------------------------------------------------------
position pointPlaneDistance(const StAxis &axis, const StVector &point)
{
StVector a(point.x-axis.xo,point.y-axis.yo,point.z-axis.zo);
StVector v(axis.a,axis.b,axis.c);
StVector a_x_v{crossProduct(a,v)};

return a_x_v.getModule()/v.getModule();
}
//-----------------------------------------------------------------------------


////  D=Ax(BxC)
/// {-a_2 b_2 c_1 - a_3 b_3 c_1 + a_2 b_1 c_2 + a_3 b_1 c_3,
/// a_1 b_2 c_1 - a_1 b_1 c_2 - a_3 b_3 c_2 + a_3 b_2 c_3,
/// a_1 b_3 c_1 + a_2 b_3 c_2 - a_1 b_1 c_3 - a_2 b_2 c_3}
StVector crossProductTriple(StVector &a, StVector &b, StVector &c)
{
cpos d_1= -a._2*b._2*c._1 - a._3*b._3*c._1 + a._2*b._1*c._2 + a._3*b._1*c._3;
cpos d_2=  a._1*b._2*c._1 - a._1*b._1*c._2 - a._3*b._3*c._2 + a._3*b._2*c._3;
cpos d_3=  a._1*b._3*c._1 + a._2*b._3*c._2 - a._1*b._1*c._3 - a._2*b._2*c._3;
        //std::cout<<"cpt " << d_1<<", "<<d_2<<", "<<d_3<<std::endl;

return StVector(d_1,d_2,d_3);
}

//-----------------------------------------------------------------------------

cpos cosa(const StVector &a, const StVector &b)
{
cpos sum=a.x*b.x+a.y*b.y+a.z*b.z;
return sum/(a.getModule()*b.getModule());
}

cpos tripleProduct(StVector &a, StVector &b, StVector &c)
{
return a.x*(b.y*c.z-b.z*c.y)+a.y*(b.z*c.x-b.x*c.z)+a.z*(b.x*c.y-b.y*c.x);
}

//-----------------------------------------------------------------------------
StVector StMatrix::operator *(const StVector &a)
{
StVector o;

        o.x=m11*a.x+m12*a.y+m13*a.z;
        o.y=m21*a.x+m22*a.y+m23*a.z;
        o.z=m31*a.x+m32*a.y+m33*a.z;

return o;
}

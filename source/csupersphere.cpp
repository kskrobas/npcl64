#include "nanograin.h"

/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * CSuperSphere.cpp
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


///  Nanomaterials 2016, 6, 27; doi:10.3390/nano6020027


//inline double pfun()

typedef const double cdouble;

inline double powf(cdouble &x, cdouble &p)
{
return std::pow(std::fabs(x),p);
}

//................................................
bool CCubic::isValid(const position &x, const position &y, const position &z)
{    
cdouble xp=powf(x/xa,p);
cdouble yp=powf(y/yb,p);
cdouble zp=powf(z/zc,p);

return (xp+yp+zp)<=Rp;
}
//................................................//................................................

bool COctahedral::isValid(const position &x, const position &y, const position &z)
{
cdouble xn=x/xa;
cdouble yn=y/yb;
cdouble zn=z/zc;


cdouble a=pow( xn+yn+zn,p);
cdouble b=pow(-xn+yn+zn,p);
cdouble c=pow( xn-yn+zn,p);
cdouble d=pow( xn+yn-zn,p);

return (a+b+c+d)<=Rp;
}
//................................................//................................................

bool CDodecahedral::isValid(const position &x, const position &y, const position &z)
{

cdouble xn=x/xa;
cdouble yn=y/yb;
cdouble zn=z/zc;

cdouble a=pow( xn+yn,p);
cdouble b=pow( xn-yn,p);
cdouble c=pow( yn+zn,p);
cdouble d=pow( yn-zn,p);
cdouble e=pow( xn+zn,p);
cdouble f=pow( xn-zn,p);

return (a+b+c+d+e+f)<=Rp;
}
//................................................//................................................
CPolyhedral::CPolyhedral(const double p__, const double R__,const double a__,const double b__) :
        CSuperSphere("polyhedral structure",p__,R__),a(a__),b(b__)
{
    iap=1.0/std::pow(a,p);
    ibp=1.0/std::pow(b,p);
    info+=", (a,b)=("+std::to_string(a)+","+std::to_string(b)+")";
}


//................................................//................................................
CPolyhedral::CPolyhedral(const double p__, const double R__,const double a__,const double b__,cdouble xa__, cdouble yb__,cdouble zc__) :
        CSuperSphere("polyhedral structure",p__,R__,xa__,yb__,zc__),a(a__),b(b__)
{
    iap=1.0/std::pow(a,p);
    ibp=1.0/std::pow(b,p);
    info+=", (a,b)=("+std::to_string(a)+","+std::to_string(b)+")";
}
//................................................//................................................
bool CPolyhedral::isValid(const position &x, const position &y, const position &z)
{

///cubic (hexa)
/// 
///

cdouble xn=x/xa;
cdouble yn=y/yb;
cdouble zn=z/zc;


cdouble xp{std::pow(fabs(xn),p)};
cdouble yp=std::pow(fabs(yn),p);
cdouble zp=std::pow(fabs(zn),p);
cdouble Hcubic=(xp+yp+zp);

/// octa

cdouble oa=pow( fabs(xn+yn+zn),p);
cdouble ob=pow( fabs(-xn+yn+zn),p);
cdouble oc=pow( fabs(xn-yn+zn),p);
cdouble od=pow( fabs(xn+yn-zn),p);
cdouble Hocta=(oa+ob+oc+od)*iap;

///doda

cdouble da=pow( fabs(xn+yn),p);
cdouble db=pow( fabs(xn-yn),p);
cdouble dc=pow( fabs(yn+zn),p);
cdouble dd=pow( fabs(yn-zn),p);
cdouble de=pow( fabs(xn+zn),p);
cdouble df=pow( fabs(xn-zn),p);
cdouble Hdoda=(da+db+dc+dd+de+df)*ibp;

///

return (Hcubic+Hocta+Hdoda)<=Rp;
}
//......................................................................................................
bool CPolyhedral2D::isValid(const position &x, const position &y, const position &z)
{
cdouble xp=powf(x/xa,p);
cdouble yp=powf(y/yb,p);

return (xp+yp)<=Rp;
}

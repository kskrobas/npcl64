#include "crandom.h"


/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * crandom.cpp
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


Crandom::Crandom(const double &A, const double &B): prmA(A),prmB(B)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    generator.seed(seed);
}

Crandom::~Crandom()
{ }


CrandomUni::CrandomUni(const double &A, const double &B): Crandom(A,B)
{
    std::uniform_real_distribution<double>::param_type params(A,B);

    rdistr.param(params);
}

CrandomNor::CrandomNor(const double &A, const double &B): Crandom(A,B)
{
    std::normal_distribution<double>::param_type params(A,B);

    rdistr.param(params);
}


CrandomLogNor::CrandomLogNor(const double &A, const double &B): Crandom(A,B)
{
    std::lognormal_distribution<double>::param_type params(A,B);

    rdistr.param(params);
}


//===================================================================================
inline double sqr(cdouble x){return x*x;}


double uniform(cdouble &x, cdouble &a, cdouble &b)
{
    if( a<=x  && x<=b) return 1.0;

    ////(b-a+1.0);

return 0;
}


double normal(cdouble &x, cdouble &mi, cdouble &sig)
{
cdouble dm=1.0/sig/std::sqrt(2.0*M_PI);
cdouble fexp=sqr((x-mi)/sig)*(-0.5);
return dm*std::exp(fexp);
}


double lognormal(cdouble &x, cdouble &m, cdouble &s)
{

cdouble dm=1.0/(s*x*std::sqrt(2*M_PI)) ;
cdouble fexp=sqr( (std::log(x)-m)/s)*(-0.5);

return dm*std::exp(fexp);
}

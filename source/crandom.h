/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * crandom.h
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

#ifndef CRANDOM_H
#define CRANDOM_H

#include <random>
#include <chrono>



using namespace  std;


typedef const double cdouble;

//===================================================================================
class Crandom
{

protected:

    const double prmA;
    const double prmB;
    std::default_random_engine generator;

public:
    Crandom():prmA(0),prmB(0) { }
    Crandom(const double &A, const double &B);

    virtual double randNumber()=0;
    virtual ~Crandom();
};
//.........................................................................................
class CrandomIgnore: public Crandom
{
private:

public:
    double randNumber(){ return 0;}
};

//.........................................................................................
class CrandomUni: public Crandom
{
private:

    std::uniform_real_distribution<double> rdistr;

public:

    CrandomUni(const double &A, const double &B);

    double randNumber(){ return rdistr(generator);}

};

//.........................................................................................
class CrandomNor: public Crandom
{
private:

    std::normal_distribution<double> rdistr;

public:

    CrandomNor(const double &A, const double &B);

    double randNumber(){ return rdistr(generator);}

};
//.........................................................................................
class CrandomLogNor: public Crandom
{
private:

    std::lognormal_distribution<double> rdistr;

public:

    CrandomLogNor(const double &A, const double &B);

    double randNumber(){ return rdistr(generator);}

};
//.........................................................................................



double uniform(cdouble &x, cdouble &a, cdouble &b);
double normal(cdouble &x, cdouble &mi, cdouble &sig);
double lognormal(cdouble &x, cdouble &m, cdouble &s);




#endif // CRANDOM_H

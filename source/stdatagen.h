/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * stdatagen.h
 * Copyright (C) 2019 Kazimierz Skrobas <kskrobas@unipress.waw.pl>
 *
 * npcl is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * mcp is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef STDATAGEN_H
#define STDATAGEN_H

#include <iostream>
using namespace std;

template <class T>
struct StDataGeneric
{

//#pragma pack(push, 1)
T *values=nullptr;
size_t size=0;
bool empty=true;


    StDataGeneric(){
        values=nullptr;
        size=0;
        empty=true;
    }

    StDataGeneric(const StDataGeneric &data)
    {
        allocMem(data.size);

        for(size_t i=0;i<data.size;i++)
            values[i]=data.values[i];
    }

    StDataGeneric( StDataGeneric &&data)
    {
        size=data.size;
        values=std::move(data.values);
        empty=false;

        data.size=0;
        data.values=nullptr;
        data.empty=true;
    }


    StDataGeneric(const size_t &size){        
        allocMem(size);
    }

    ~StDataGeneric(){
        freeMem();
    }


    StDataGeneric & operator=(StDataGeneric && data)
    {
        size=data.size;
        values=std::move(data.values);
        empty=false;

        data.size=0;
        data.values=nullptr;
        data.empty=true;
    return *this;
    }



    void freeMem()
    {
        if(size){
            delete [] values;
            size=0;
            values=0;
            empty=true;
        }
    }


    bool allocMem(const size_t & size)
    {
        freeMem();

        if(size){
            this->size=size;
            values=new T[size];
            empty=false;
            return true;
        }
        else
          return false;       
    }


    bool allocMem(const size_t &size,const T iniVal)
    {

         freeMem();

          if(size){
            this->size=size;
            values=new T[size];
            empty=false;

            for(size_t i=0;i<size;i++)
                values[i]=iniVal;

           // for (auto &v:values)
             //       v=iniVal;

          return true;
          }
          else
          return false;

    }


     T & operator[](const size_t el) {return values[el];}
//#pragma pack(push)
};


#endif // STDATAGEN_H

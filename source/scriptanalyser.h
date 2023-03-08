/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * scriptanalyser.h
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



#ifndef SCRIPTANALYSER_H
#define SCRIPTANALYSER_H

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>


using namespace std;
typedef std::pair<string, string> strpair;


class ClKeyValues{
private:

        vector<string> keyvalues;
        string *ptr_frontValue;
        string *ptr_lastValue;

        void fini();
        void fcon();
        void (ClKeyValues::*ptr_f)();

public:

        ClKeyValues(){ keyvalues.reserve(8);}
        ClKeyValues(const char *key){keyvalues.reserve(8);keyvalues.emplace_back(string(key));}
        ClKeyValues(const string key){keyvalues.reserve(1);keyvalues.emplace_back(string(key));}
        ClKeyValues(strpair keyValue){keyvalues.reserve(2);keyvalues.emplace_back(keyValue.first);keyvalues.emplace_back(keyValue.second);}

        ClKeyValues & operator << (string &kv);
        ClKeyValues & operator << (const char *);
        ClKeyValues & operator << (const vector<string> &vs ) ;
        ClKeyValues & operator >> (string &kv);
        ClKeyValues & operator >> (vector<string> &vs ) ;

        friend ostream& operator<< (ostream &o, const ClKeyValues &kv);

        ClKeyValues operator = (const ClKeyValues &b);
       // ClKeyValues operator = (ClKeyValues &&);
        ClKeyValues(const ClKeyValues &);
        ClKeyValues( ClKeyValues &&);

        void dispValues();
        int numOfKeyValues(){return keyvalues.size();}
        void clear(){ keyvalues.clear();}
        bool empty(){ return keyvalues.empty();}


        bool isKey(std::string str__)
        {   return  keyvalues[0].find(str__)!=std::string::npos;}

        const std::string & getKey(){return keyvalues[0];}
        std::string & getValue(const size_t value=1){return keyvalues[value];} // dopuszczalna modyfikacja wartosci

        bool operator<(const ClKeyValues &b)
            {return keyvalues[0]<b.keyvalues[0];}

        bool operator==(const std::string &s)
            {return keyvalues[0]==s;}
        bool operator==(const char *key__)
            {std::string skey(key__);
             return keyvalues[0]==skey;}
        bool operator == (const ClKeyValues &b)
            {return keyvalues[0]==b.keyvalues[0];}


        bool operator!=(const char *key__)
            {std::string skey(key__);
             return keyvalues[0]!=skey;}

        bool operator!=(const std::string &s)
            {return keyvalues[0]!=s;}

        string & operator [](size_t i)
        {    return keyvalues[i];}

};


typedef vector<ClKeyValues> vcmdlist;
typedef vector<ClKeyValues>::iterator iter_vcmdlist;
typedef vector<ClKeyValues> stdumap;

//-----------------------------------------------------------------------------

namespace Script{

enum Result{OK,WARN,ERR0,ERR1,BLOCK_ERR,MIS_UNK_ERR,DUP_ERR,ERR_EXPR,ERR_OUTOFRANGE,
            ERR_UNK_VAR,ENDRET,ENDELSE,ERR_MATH_EXP,ERR_PRT,ERR_VAL_0,ERR_NO_NUMBER,
            ERR_INV_PAR,ERR_INF_LOOP,ERR_TRVUC,ERR_NAUC,ERR_EOF};
enum ResultRepVar{EMPTY,SUCC,FAIL};

Result scriptParsing(fstream &file, size_t &cline, vcmdlist *ptr_cl,stdumap *ptr_uvar,const size_t options);
//ResultRepVar findVariable(stdumap *ptr_uvar);
ResultRepVar replaceVars(stdumap *ptr_uvar,std::string &value);

};


//-----------------------------------------------------------------------------
// trim from start
static inline std::string &ltrim(std::string &s)
{
auto isNotSpace=[](char c){return !std::isspace(c);};
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), isNotSpace));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s)
{
auto isNotSpace=[](char c){return !std::isspace(c);};
        s.erase(std::find_if(s.rbegin(), s.rend(), isNotSpace).base(), s.end());
        return s;
}




// trim from both ends
static inline std::string &trim(std::string &s)
{
        return ltrim(rtrim(s));
}

//-----------------------------------------------------------------------------
template<typename T>
vector<T>  split(const T & str, const T & delimiters)
{
vector<T> v;
typename T::size_type start = 0;
auto pos = str.find_first_of(delimiters, start);

        while(pos != T::npos) {
            if(pos != start) // ignore empty tokens
                v.emplace_back(str, start, pos - start);
            start = pos + 1;
            pos = str.find_first_of(delimiters, start);
        }
        if(start < str.length()) // ignore trailing delimiter
            v.emplace_back(str, start, str.length() - start); // add what's left of the string
return v;
}
//-----------------------------------------------------------------------------


#endif // SCRIPTANALYSER_H

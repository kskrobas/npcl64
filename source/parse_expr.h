/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * parse-expr.h
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



#ifndef PARSE_EXPR_H
#define PARSE_EXPR_H


#include <string>
#include <iostream>
using namespace  std;


//---------------------------------------------------------------------------

class CParseABC{
private:

        struct StExpr{
        size_t start=0,stop=0;
        std::string token;

            StExpr(){ }
            int Length(){return stop-start;}
        };


        std::string expression;
        size_t exprLength;

        //---------------------------------------------------------------------------
         StExpr buildABCexpand(const unsigned startPos);
         StExpr buildABC(const unsigned startPos=0);
        //---------------------------------------------------------------------------

public:
        enum Results{OK,ERR_POS,ERR_BRA,ERR_ILLCHAR};
        bool run(const string &expr,string &result);
};

#endif // PARSE_EXPR_H

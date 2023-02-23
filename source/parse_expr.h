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

using namespace  std;


//---------------------------------------------------------------------------
//class TExprException: public exception
//{
//private:
//        const char *descr[4]={"unknown","unbalanced brackets","invalid expression","invalid pointer position"};
//        unsigned ce=0;
//        unsigned pos=0;
//public:
//        enum ErrorCode{UNKNOWN,UNBRACKETS,EXPRESSION,POSITION};
//        TExprException(){ }
//        TExprException(const ErrorCode ce_){ce=ce_;}
//        TExprException(const ErrorCode ce_,const unsigned pos_){ce=ce_;pos=pos_;}

//        virtual const char* what() const throw()
//        {
//        return descr[ce];
//        }

//        const	unsigned & getPosition() {return pos;}
//} ;


class CParseABC{
private:


            struct StExpr{
            size_t start=0,stop=0;
            std::string token;

                StExpr(){ }
                int Length(){return stop-start;}
            };


enum Results{OK,ERR_POS,ERR_BRA};

        std::string expression;
        size_t exprLength;

            //---------------------------------------------------------------------------
            StExpr buildABCexpand(const unsigned startPos)
            {
            const string &str=expression;
            size_t pos=startPos;
            StExpr expr;
            string digits;

                    expr.start=startPos;

                    while(pos<exprLength && isdigit(str[pos]) )
                        digits+=str[pos++];

            string letters;

                    if(isalpha(str[pos]))
                        letters=str[pos];
                    else
                        if(str[pos]=='(') {
                        StExpr exprUp=buildABC(pos);
                                pos+=exprUp.Length();
                                letters=exprUp.token;

                        }
                        else
                        throw Results::ERR_POS;

                    expr.stop=pos;

            const size_t intNum=std::stoi(digits);

                        expr.token="";
                        for(size_t i=0;i<intNum;i++)
                            expr.token+=letters;

            return expr;
            }

            //---------------------------------------------------------------------------
            StExpr buildABC(const unsigned startPos=0)
            {
            const string &str=expression;
            size_t pos=startPos;
            StExpr expr;

                    expr.start=pos;
                    pos++;

                    while( pos<exprLength && str[pos]!=')'){

                        if(isalpha(str[pos]))
                            expr.token+=str[pos];
                        else{
                            if(isdigit(str[pos])){
                            StExpr exprUp=buildABCexpand(pos);
                                pos+=exprUp.Length();
                                expr.token+=exprUp.token;
                            }
                            else{
                            StExpr exprUp=buildABC(pos);
                                pos+=exprUp.Length();
                                expr.token+=exprUp.token;
                            }
                        }
                        pos++;
                    }

                    if(pos>exprLength)
                    throw Results::ERR_POS;

                    expr.stop=pos;

            return	expr;
            }
            //---------------------------------------------------------------------------


public:

            bool run(const string &expr,string &result){

                    expression="("+expr+")";
                    exprLength=expression.length();

                    try{
                    StExpr expr=buildABC();

                        if( expr.stop != exprLength-1) //unbalanced braces detection
                            throw  Results::ERR_BRA;

                         result=expr.token;

                        return true;
                    }
                    catch(Results r){
                        cerr<<"error: ABC expression fault ";

                        switch(r){
                        case Results::ERR_BRA : cerr<<"unbalanced brackets"<<endl;break;
                        case Results::ERR_POS : cerr<<""<<endl;break;
                        default: cerr<<"OK"<<endl;
                        }

                        return false;
                    }
            }

};

#endif // PARSE_EXPR_H

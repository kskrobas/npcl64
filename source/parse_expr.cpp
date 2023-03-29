#include "parse_expr.h"


//---------------------------------------------------------------------------
bool  CParseABC::run(const string &expr,string &result){

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
        cerr<<"ERROR: ABC expression fault "<<expr<<endl;

        switch(r){
        case Results::ERR_BRA :    cerr<<"        unbalanced brackets"<<endl;break;
        case Results::ERR_POS :    cerr<<"        position failure "<<endl;break;
        case Results::ERR_ILLCHAR: cerr<<"        illegal character"<<endl;break;
        default: cerr<<" unidentified "<<endl;
        }

        return false;
    }
}

//---------------------------------------------------------------------------
CParseABC::StExpr CParseABC::buildABCexpand(const unsigned startPos)
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
CParseABC::StExpr CParseABC::buildABC(const unsigned startPos)
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
                    if(str[pos]=='('){
                    StExpr exprUp=buildABC(pos);
                        pos+=exprUp.Length();
                        expr.token+=exprUp.token;
                    }
                    else{

                    throw CParseABC::Results::ERR_ILLCHAR;
                    }

                }
            }
            pos++;
        }

        if(pos>exprLength)
        throw Results::ERR_POS;

        expr.stop=pos;

return	expr;
}


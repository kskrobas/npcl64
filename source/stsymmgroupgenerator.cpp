#include "stsymmgroupgenerator.h"

position StAtomFracPostion::tol=1e-4;

//-----------------------------------------------------------------------------
inline position sqr(const position x){return x*x; }

//-----------------------------------------------------------------------------
struct StRowShift
{
        union{
            struct{position row[3],shift;};
            position rs[4];
        };

        StRowShift(){ rs[0]=rs[1]=rs[2]=rs[3]=0;}
};
//-----------------------------------------------------------------------------


StSymmGroupGenerator::StSymmGroupGenerator()
{
    for(size_t i=0;i<12;i++) mx12[i]=0;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------

float decodeNumber(const string &formula,const size_t & i)
{
string v(formula.substr(i,3));

        if(v=="1/4")   return 0.25;
        if(v=="1/3")   return 0.3333333;
        if(v=="1/2")   return 0.5;
        if(v=="2/3")   return 0.6666666;
        if(v=="3/4")   return 0.75;

return std::stof(v.substr(0))/std::stof(v.substr(2));
}
//-----------------------------------------------------------------------------
StRowShift decodeRow(const string & formula)
{
StRowShift rowSh;
auto &m=rowSh.row;
position &shift=rowSh.shift;
int pm=1;
char ch;


            for(size_t i=0;i<formula.length();i++){

                    ch=formula[i];
                    if(ch==' ' || ch=='\'') continue;
                    if(ch=='+' || ch=='-') { pm=','-ch; continue;}
                    if(std::isdigit(ch)) { shift=pm*decodeNumber(formula,i); i+=2; continue;}

                    m[ch-'x']=pm;
            }

return rowSh;
}
//-----------------------------------------------------------------------------

StSymmGroupGenerator buildGenerator(const string &formula)
{
vector<string> tokens(split<string>(formula,","));
StSymmGroupGenerator gnr;
auto & velem=gnr.mx12;
size_t i=0,j;

            for(auto & token: tokens){
                for(j=0;j<4;j++,i++)
                    velem[i]=decodeRow(token).rs[j];
            }
return gnr;
}

bool StAtomFracPostion::operator ()(const StAtomFracPostion &v)
{
auto r2=sqr(v.x-x)+sqr(v.y-y)+sqr(v.z-z);
return r2<tol;
}

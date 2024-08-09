#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

typedef vector<int> vint;
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
struct StRowShift
{
    union{
        struct{float row[3],shift;};
        float rs[4];
    };

    StRowShift(){ rs[0]=rs[1]=rs[2]=rs[3]=0;}
};
//-----------------------------------------------------------------------------
StRowShift decodeRow(const string & formula)
{
StRowShift rowSh;
auto &m=rowSh.row;
float &shift=rowSh.shift;
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
inline float sqr(const float x){return x*x; }


struct StVector
{
static float tol;
float x,y,z;
size_t id;

        StVector(){ x=y=z=0;}
        StVector(const float x__, const float y__, const float z__)
        {x=x__;y=y__;z=z__;}
        StVector(const float x__, const float y__, const float z__, const size_t id__)
        {x=x__;y=y__;z=z__;id=id__;}


        friend ostream & operator <<(ostream &, StVector &);

        void reduct()
        {
            if(x<0) x+=1;
            if(y<0) y+=1;
            if(z<0) z+=1;
        }



        bool operator == (const StVector & vin) const
        {
        auto r2=sqr(vin.x-x)+sqr(vin.y-y)+sqr(vin.z-z);
        return r2<tol;
        }

        bool operator () (StVector &v) { return *this == v; }
};

float StVector::tol=1e-4;

const int sw=6;
const int sp=3;
ostream & operator <<(ostream & o, StVector &stv)
{
    o<<":"<<stv.id<<" "<<
       "[ " <<setw(sw)<<setprecision(sp)<<stv.x<<
       "," <<setw(sw)<<setprecision(sp)<<stv.y<<
       "," <<setw(sw)<<setprecision(sp)<<stv.z<<" ]";
return o;
}

//-----------------------------------------------------------------------------

struct StGenerator
{
union{
    float elems3x4[3][4];
    float elems12[12];
    struct{ float m11,m12,m13,m14;
            float m21,m22,m23,m24;
            float m31,m32,m33,m34;};
    };

float lpa,lpb,lpc;

    StGenerator(){ for(size_t i=0;i<12;i++) elems12[i]=0;}

    StVector ucAtomPos(StVector & ap)
    {
    StVector v;
            v.x=ap.x*m11+ap.y*m12+ap.z*m13+m14;
            v.y=ap.x*m21+ap.y*m22+ap.z*m23+m24;
            v.z=ap.z*m31+ap.y*m32+ap.z*m33+m34;
    return v;
    }

    friend ostream & operator<<(ostream &, StGenerator &);

};


ostream &operator <<(ostream & o, StGenerator &sg)
{
    /*for(size_t i=0;i<3;i++){
        o<<endl;
        for(size_t j=0;j<4;j++)
            cout<<"  "<<sg.elems3x4[i][j];
    }*/

    o<<"┌ "<<setw(27)<<" "<<" ┐\n";
    o<<"│ ";for(size_t i=0;i<3;i++) o<<setw(sw)<<setprecision(sp)<<sg.elems3x4[0][i]<< ",";  o<<setw(sw)<<setprecision(sp)<<sg.elems3x4[0][3]<<" │\n";
    o<<"│ ";for(size_t i=0;i<3;i++) o<<setw(sw)<<setprecision(sp)<<sg.elems3x4[1][i]<< ",";  o<<setw(sw)<<setprecision(sp)<<sg.elems3x4[1][3]<<" │\n";
    o<<"│ ";for(size_t i=0;i<3;i++) o<<setw(sw)<<setprecision(sp)<<sg.elems3x4[2][i]<< ",";  o<<setw(sw)<<setprecision(sp)<<sg.elems3x4[2][3]<<" │\n";
    o<<"└ "<<setw(27)<<" "<<" ┘";

return o;
}




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


StGenerator buildGenerator(const string &formula)
{
//const string trimFormula(formula.substr(1,formula.length()-2));
vector<string> tokens(split<string>(formula,","));
StGenerator gnr;
auto & velem=gnr.elems12;
size_t i=0,j;

                for(auto & token: tokens){
                    for(j=0;j<4;j++,i++)
                        velem[i]=decodeRow(token).rs[j];
                }
return gnr;
}

//-----------------------------------------------------------------------------


int main()
{

fstream fileIn("rutile.txt",ios::in);
                if(!fileIn){
                    cerr<<" error , file not found "<<endl;
                return -1;
                }

vector<string> sgOperators;
vector<StVector> atomPosSet;
size_t numOfatoms;
string fline;



                fileIn>>numOfatoms;std::getline(fileIn,fline);

                atomPosSet.resize(numOfatoms);

                for(size_t i=0;i<numOfatoms;i++){
                    std::getline(fileIn,fline);

                vector<string> tokens(split<string>(fline," "));
                    if(tokens.size()!=3) {
                        cerr<<" error, wrong format, line: "<<i+1<<endl;
                    return -1;
                    }

                auto &apos=atomPosSet[i];

                    apos.x=std::stof(tokens[0]);
                    apos.y=std::stof(tokens[1]);
                    apos.z=std::stof(tokens[2]);
                    apos.id=i;

                }


                sgOperators.reserve(200);

                cout<<" sg operators ";
                while(!fileIn.eof() ){

                    std::getline(fileIn,fline);
                    if(fline.empty()) continue;

                    sgOperators.push_back(fline);
                    cout<<"\n  "<<fline;
                    fline.clear();
                }

                cout<<endl;
                fileIn.close();
                sgOperators.shrink_to_fit();


StVector ucAtomPos;
vector<StVector> acceptedAtoms;
//auto cmpFunc=[&ucAtomPos](const StVector &v){ return v==ucAtomPos; };


                    //atomPosSet.emplace_back(StVector(0,0,0,0));
                    //atomPosSet.emplace_back(StVector(1.0f/3.0f,2.0f/3.0f,0));
                    //atomPosSet.emplace_back(StVector(0.25,0.25,0.25,1));

                    acceptedAtoms.reserve((sgOperators.size())*atomPosSet.size());


                    for(auto &ap: atomPosSet){
                        for(auto & sgFormula: sgOperators){
                        StGenerator stg(buildGenerator(sgFormula));


                            cout<<"\n----------------------------"
                                <<"\n sg formula: "<<sgFormula
                                <<"\n"<<stg;

                            cout<<"\n atomPositons for: "<<ap<<" -> ";

                            ucAtomPos=stg.ucAtomPos(ap);
                            cout<<ucAtomPos<<endl;

                            ucAtomPos.reduct();
                            auto itr=std::find_if(acceptedAtoms.begin(),acceptedAtoms.end(),ucAtomPos);

                            if(itr==acceptedAtoms.end()){
                                acceptedAtoms.push_back(ucAtomPos);
                                acceptedAtoms.back().id=ap.id;
                            }
                        }
                    }

                acceptedAtoms.shrink_to_fit();

                cout<<"-------------------------------------"<<endl;
                cout<<"accepted positions: "<<acceptedAtoms.size()<<endl;
                for(auto &ap: acceptedAtoms){
                    cout<<ap<<endl;
                }


//            buildGenerator(genFormula);


return 0;
}

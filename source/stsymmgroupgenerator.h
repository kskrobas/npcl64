#ifndef STSYMMGROUPGENERATOR_H
#define STSYMMGROUPGENERATOR_H

#include "nanograin.h"

struct StAtomFracPostion
{
static position tol;
position x,y,z;
string name;

        StAtomFracPostion() { }
        StAtomFracPostion(const position x__, const position y__, const position z__)
            :x(x__),y(y__),z(z__) { }

        void reduct()
        {
            if(x<0) x+=1; if(x>1) x-=1;
            if(y<0) y+=1; if(y>1) y-=1;
            if(z<0) z+=1; if(z>1) z-=1;
        }

        bool operator () (const  StAtomFracPostion &v) ;
      //  bool operator () (const StVector v) ;
};






class StSymmGroupGenerator
{
public:


    union{
        float mx3x4[3][4];
        float mx12[12];
        struct{ float m11,m12,m13,m14;
                float m21,m22,m23,m24;
                float m31,m32,m33,m34;};
        };

        StSymmGroupGenerator();

        StAtomFracPostion ucAtomPos(StAtomFracPostion  ap)
        {
        StAtomFracPostion v;
                v.x=ap.x*m11+ap.y*m12+ap.z*m13+m14;
                v.y=ap.x*m21+ap.y*m22+ap.z*m23+m24;
                v.z=ap.z*m31+ap.y*m32+ap.z*m33+m34;
        return v;
        }

     //   friend ostream & operator<<(ostream &, StSymmGroupGenerator &);
};




StSymmGroupGenerator buildGenerator(const string &formula);

#endif // STSYMMGROUPGENERATOR_H

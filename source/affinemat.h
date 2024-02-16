#ifndef AFFINEMAT_H
#define AFFINEMAT_H

#include "nanograin.h"

typedef double position;
typedef const position cpos;

//-----------------------------------------------------------------------------
//---------------------------------------------------------------------------
struct StAxis{
union{
    struct{
    position a,b,c,xo,yo,zo;
    };
    position prm[6];
    };
bool on=false;
double tmax,tmin;

//StVector3 begin,end;

    StAxis(cpos &x, cpos &y, cpos &z){
        a=x;b=y;c=z;
    }

};


class StRotationMatrix{
public:

union{
      struct{ position m11,m12,m13,m21,m22,m23,m31,m32,m33;};
      position m[9];
      position mm[3][3];
};

position ux,uy,uz;//axis of rotation
position theta;

bool on=false;

    void buildMatrix(const StAxis &axis_);
    void buildMatrix(const StAxis &axis_,cpos &sinA);
    StVector operator*(const StVector &v);
//StAtom operator*=(S);

    StRotationMatrix(){ ux=uy=uz=theta=0;}
    StRotationMatrix(const StAxis &axis_,cpos &sinA){
        buildMatrix(axis_,sinA);
    }

private:
    void buildRotationAxis(const StAxis &axis_);
    void normUxyz();
};
//-----------------------------------------------------------------------------


#endif // AFFINEMAT_H

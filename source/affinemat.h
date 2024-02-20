#ifndef AFFINEMAT_H
#define AFFINEMAT_H

#include <cmath>

typedef double position;
typedef const position cpos;

//-----------------------------------------------------------------------------

struct StVector{
position x,y,z;


      StVector operator+(const StVector & op) const
      {
       StVector res;

                res.x=this->x+op.x;
                res.y=this->y+op.y;
                res.z=this->z+op.z;

       return res;
      }

      StVector operator-(const StVector & op) const
      {
       StVector res;

                res.x=this->x-op.x;
                res.y=this->y-op.y;
                res.z=this->z-op.z;

       return res;
      }

      StVector operator*(const double &val) const
      {
       StVector res;

                res.x=this->x*val;
                res.y=this->y*val;
                res.z=this->z*val;

       return res;
      }


      //int operator = (const StGrainNode &node);
      //int operator +=(const StGrainNode &node);
      StVector() {  }
      StVector(cpos &x__, cpos &y__, cpos &z__):x(x__),y(y__),z(z__) { }


      double getModule()const { return sqrt(x*x+y*y+z*z);}

};


//---------------------------------------------------------------------------
struct StAxis{
union{
    struct{position a,b,c,xo,yo,zo;};
    position prm[6];
    };
bool on=false;
double tmax,tmin;


    StAxis(cpos a__,cpos b__, cpos c__, cpos xo__, cpos yo__, cpos zo__)
    : a(a__),b(b__),c(c__),xo(xo__),yo(yo__),zo(zo__) {  }

    StAxis(cpos &x, cpos &y, cpos &z){
        a=x;b=y;c=z;
    }

    double getModule() const {return sqrt(a*a+b*b+c*c);}

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


StVector crossProduct(StVector &a,StVector &b);
position projHeight(const StAxis &a, const StVector &b);
position pointPlaneDistance(const StAxis &axis,const StVector &point);



#endif // AFFINEMAT_H

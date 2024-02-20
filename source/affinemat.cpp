#include "affinemat.h"

inline cpos sqr(cpos &x){return x*x;}


void StRotationMatrix::buildMatrix(const StAxis &axis_)
{

}

///
///https://en.wikipedia.org/wiki/Rotation_matrix
///

void StRotationMatrix::buildMatrix(const StAxis &axis_, cpos &sinA)
{
cpos sinTheta=sinA;//sqrt(1-sqrCoord(cosTheta));
cpos cosTheta=std::sqrt(1-sqr(sinTheta));
cpos oneMinCos=1-cosTheta;

        buildRotationAxis(axis_);

        ///rotation matrix elements
        m11=cosTheta+sqr(ux)*oneMinCos; m12=ux*uy*oneMinCos-uz*sinTheta; m13=ux*uz*oneMinCos+uy*sinTheta;
        m21=uy*ux*oneMinCos+uz*sinTheta; m22=cosTheta+sqr(uy)*oneMinCos; m23=uy*uz*oneMinCos-ux*sinTheta;
        m31=uz*ux*oneMinCos-uy*sinTheta; m32=uz*uy*oneMinCos+ux*sinTheta; m33=cosTheta+sqr(uz)*oneMinCos;

        on=true;
}

//-----------------------------------------------------------------------------
StVector StRotationMatrix::operator*(const StVector &a)
{
position bx,by,bz;

                bx=m11*a.x+m12*a.y+m13*a.z;
                by=m21*a.x+m22*a.y+m23*a.z;
                bz=m31*a.x+m32*a.y+m33*a.z;

return StVector(bx,by,bz);
}

//-----------------------------------------------------------------------------
void StRotationMatrix::buildRotationAxis(const StAxis &axis_)
{
        ux=axis_.a;
        uy=axis_.b;
        uz=axis_.c;

        normUxyz();
}

//-----------------------------------------------------------------------------
void StRotationMatrix::normUxyz()
{
cpos sumSq=sqr(ux)+sqr(uy)+sqr(uz);
cpos norm=sqrt(1/sumSq);

            ux*=norm;
            uy*=norm;
            uz*=norm;
}
//-----------------------------------------------------------------------------
StVector crossProduct(StVector &a, StVector &b)
{
cpos cx=(a.y*b.z-a.z*b.y);
cpos cy=(a.z*b.x-a.x*b.z);
cpos cz=(a.x*b.y-a.y*b.x);

return StVector(cx,cy,cz);
}
//-----------------------------------------------------------------------------
position projHeight(const StAxis &axis, const StVector &b)
{
cpos aDotb=std::abs(axis.a*(axis.xo-b.x)+axis.b*(axis.yo-b.y)+axis.c*(axis.zo-b.z));
return aDotb/axis.getModule();
}
//-----------------------------------------------------------------------------
position pointPlaneDistance(const StAxis &axis, const StVector &point)
{
StVector a(point.x-axis.xo,point.y-axis.yo,point.z-axis.zo);
StVector v(axis.a,axis.b,axis.c);
StVector a_x_v{crossProduct(a,v)};

return a_x_v.getModule()/v.getModule();
}


/*--------------------iunit.c-------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include "iunit.h" 
#define PI 3.1415




struct point makepoint(float initx, float inity, float initz)
{  
    struct point p;

    p.x=initx;
    p.y=inity;
    p.z=initz;
    return p;
}

void printpoint(struct point p)
{
   printf(" %3f",p.x);printf(" %3f",p.y); printf(" %3f",p.z);printf("\n");   
}

struct matrix3d makematrix3d(float x11, float x12, float x13,
float x21, float x22, float x23,float x31, float x32, float x33)
{  
    struct matrix3d m;

    m.x11=x11; m.x12=x12; m.x13=x13;
    m.x21=x21; m.x22=x22; m.x23=x23;
    m.x31=x31; m.x32=x32; m.x33=x33;
    return m;
}

struct point matrixmult(struct matrix3d m,struct point p)
{
   return makepoint(m.x11*p.x+m.x12*p.y+m.x13*p.z,
                    m.x21*p.x+m.x22*p.y+m.x23*p.z,
                    m.x31*p.x+m.x32*p.y+m.x33*p.z);
}

struct point rotateeuler(struct point p,float psi,float theta,float phi)
{
   struct matrix3d m; 
   float sinpsi,cospsi,sintheta,costheta,sinphi,cosphi;

  sinpsi=sin(psi);cospsi=cos(psi);
  sintheta=sin(theta);costheta=cos(theta);
  sinphi=sin(phi);cosphi=cos(phi);
  m=makematrix3d(cosphi*cospsi-sinphi*sinpsi*costheta,
             -(sinphi*cospsi+cosphi*sinpsi*costheta),
               sinpsi*sintheta,
               cosphi*sinpsi+sinphi*cospsi*costheta,
              -sinphi*sinpsi+cosphi*cospsi*costheta,
              -cospsi*sintheta,
               sinphi*sintheta,
               cosphi*sintheta,
               costheta);
    return matrixmult(m,p);
}

struct matrix3d rotmatrix(float psi,float theta,float phi)
{
  /* This is the inverse euler matrix A-1 */
  /* The rotation is in the bodies coordinates system */
   float sinpsi,cospsi,sintheta,costheta,sinphi,cosphi;

  sinpsi=sin(psi);cospsi=cos(psi);
  sintheta=sin(theta);costheta=cos(theta);
  sinphi=sin(phi);cosphi=cos(phi);
  return makematrix3d(cosphi*cospsi-sinphi*sinpsi*costheta,
             -(sinphi*cospsi+cosphi*sinpsi*costheta),
               sinpsi*sintheta,
               cosphi*sinpsi+sinphi*cospsi*costheta,
              -sinphi*sinpsi+cosphi*cospsi*costheta,
              -cospsi*sintheta,
               sinphi*sintheta,
               cosphi*sintheta,
               costheta);
}

struct matrix3d INVrotmatrix(float phi,float theta,float psi)
{
   float sinpsi,cospsi,sintheta,costheta,sinphi,cosphi;

  sinpsi=sin(psi);cospsi=cos(psi);
  sintheta=sin(theta);costheta=cos(theta);
  sinphi=sin(phi);cosphi=cos(phi);
  return makematrix3d(
		      cosphi*cospsi-sinphi*sinpsi*costheta,
             cospsi*sinphi+costheta*cosphi*sinpsi,
               sinpsi*sintheta,
               -sinpsi*cosphi-costheta*sinphi*cospsi,
              -sinphi*sinpsi+cosphi*cospsi*costheta,
	       cospsi*sintheta,
               sinphi*sintheta,
               -cosphi*sintheta,
               costheta);
}

struct point mirror_x(struct point vec)
{
  struct point mirror;
  mirror.x = -vec.x;
  mirror.y = vec.y;
  mirror.z = vec.z;

  return mirror;
}
struct point mirror_z(struct point vec)
{
  struct point mirror;
  mirror.x = vec.x;
  mirror.y = vec.y;
  mirror.z = -vec.z;

  return mirror;
}


struct matrix3d rotz(float teta)
{
   float st,ct;

   st=sin(teta);
   ct=cos(teta);
  return makematrix3d(
		      1,0,0,
                      0,ct,st,
		      0,-st,ct);
}

struct matrix3d roty(float teta)
{
   float st,ct;

   st=sin(teta);
   ct=cos(teta);
  return makematrix3d(
		      ct,0,st,
                      0,1,0,
		      -st,0,ct);
}

struct matrix3d rotx(float teta)
{
   float st,ct;

   st=sin(teta);
   ct=cos(teta);
  return makematrix3d(
		      ct,st,0,
                      -st,ct,0,
		      0,0,1);
}

float rad(float angle)
{
    return (PI*angle)/180;
}

struct point add(struct point p1,struct point p2)
{   
  return makepoint(p1.x+p2.x,p1.y+p2.y,p1.z+p2.z);
}

struct point sub(struct point p1,struct point p2)
{   
  return makepoint(p1.x-p2.x,p1.y-p2.y,p1.z-p2.z);
}

struct point scalarmult(struct point p1,float c)
{   
    return makepoint(p1.x*c,p1.y*c,p1.z*c);
}

float scalarprod(struct point p1,struct point p2)
{      
    return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;    
}

struct point vectorproduct(struct point p1,struct point p2)
{   
    return makepoint(p1.y*p2.z-p1.z*p2.y,
            p1.z*p2.x-p1.x*p2.z,
            p1.x*p2.y-p1.y*p2.x);
}

float sqr(float c)
{
   return c*c;
} 

float betrag(struct point p1,struct point p2)
{
    return sqrt(sqr(p1.x-p2.x)+sqr(p1.y-p2.y)+sqr(p1.z-p2.z));
} 

int mymessage()
{
    printf(" das ist mein test \n");
    return (0);
}


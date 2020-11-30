/*--------------------iunit.h-------------------------------------------------*/
struct point {
  float x;
  float y;
  float z;
};

struct matrix3d {
  float x11,x12,x13,x21,x22,x23,x31,x32,x33;
};

extern struct point makepoint(float initx, float inity, float initz);

extern void printpoint(struct point p);

extern struct matrix3d makematrix3d(float x11, float x12, float x13,
float x21, float x22, float x23,float x31, float x32, float x33);

extern struct point matrixmult(struct matrix3d m,struct point p);

extern struct point rotateeuler(struct point p,float psi,float theta,float phi);
extern struct point mirror_x(struct point vec);
extern struct point mirror_z(struct point vec);
extern struct matrix3d rotz(float teta);
extern struct matrix3d roty(float teta);
extern struct matrix3d rotx(float teta);

extern struct matrix3d rotmatrix(float psi,float theta,float phi);
extern struct matrix3d INVrotmatrix(float phi,float theta,float psi);

extern float rad(float angle);

extern struct point add(struct point p1,struct point p2);

extern struct point sub(struct point p1,struct point p2);

extern struct point scalarmult(struct point p1,float c);

extern float scalarprod(struct point p1,struct point p2);

extern struct point vectorproduct(struct point p1,struct point p2);

extern float sqr(float c);

extern float betrag(struct point p1,struct point p2);

extern int mymessage();



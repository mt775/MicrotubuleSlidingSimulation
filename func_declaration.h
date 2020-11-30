#ifndef _func_declaration_f
#define _func_declaration_f


/* Function declaration */                             /* filename */


void usage_err(char *argv[]);                          /* functions.c  */ 
void rdata(int argc, char *argv[]);                    /* functions.c  */
char *filept(char *filename);                          /* functions.c  */
void init();                                           /* init.c       */
void overlapMap();                                     /* overlapMap.c */
void mkMTs(int newMTs);                                /* mkMTs.c      */
float overlap(int ia, int ib);                         /* overlap.c    */
void motorvec(int ia, int ib, struct point *motorhead, 
	      struct point *motorleg);                 /* overlap.c    */
float maxf(float a, float b);                          /* functions.c  */
float maxi(int a, int b);                              /* functions.c  */
float minf(float a, float b);                          /* functions.c  */
float mini(int a, int b);                              /* functions.c  */
struct point force(int iMT);                           /* force.c      */
struct point Fmotor(int iMT);                          /* force.c      */
struct point Fdrag(int iMT);                           /* force.c      */
struct point Frand(int iMT);                           /* force.c      */ 
void equations();                                      /* equations.c  */
void ludcmp(float **a, int n, int *indx, float *d);    /* ludcmp.c     */
void lubksb(float **a, int n, int *indx, float b[]);   /* lubksb.c     */
void mprove(float **a, float **alud, int n, int indx[], 
	    float b[], float x[]);                     /* mprove.c     */          
void avcm2(struct point *av);
void avvel2(struct point *av);
int iMTmax();
int nMTright(float x, float dx);
int nMTleft(float x, float dx);

/* Functions for writing XMGR files */
void agr_std(FILE *fp, const char *Title, const char *Xaxis, const char *Yaxis, float x, float y);
void agr_begset(FILE *fp, int k, int lncolor, int lnstyle, float lnwidth, int symbol, float symbsz, int symbcl);
void agr_endset(FILE *fp);
void agr_comment(FILE *fp, const char *comment);
void agr_wcoor(FILE *fp, float x, float y);
void agr_wcoor3D(FILE *fp, float x, float y, float z);
void print_line(FILE *output);
void short_line(FILE *output);
void agr_legend(FILE *fp, int k, const char *name);
void agr_cyl(FILE *fp,float xcm,float zcm,float r,int lw,int lc,int fpt,int fc);void agr_arrowline(FILE *fp, float pendx, float pendz, float mendx, float mendz, int icolor);       
         
#endif

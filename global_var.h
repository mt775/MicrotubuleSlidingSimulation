/* Constants */
//#define MXMTS 600
//#define MXMOTOR 100000
//#define MXITER 70010
/*For Valgrind (debugger) substitute by following*/
//#define MXMTS 302
//#define MXMOTOR 1001
//#define MXITER 10001
//#define MXITER1 1
//#define MXBOX 150
//#define MXBOX2 150
#define microM 1
#define pN 1
#define mM 1000*microM
#define nM 0.001*microM
#define kB 1.3806503E-5 /* pN microM/k */
#define PI 3.14159265359
#define EPSI 0.001 /*minimal overlap length*/

#define before printf("Before, iter=%d\n",iter)
#define after printf("Aefore, iter=%d\n",iter)



float perc_prob;
float degprob;
int mxmts,mxmotors;

/******************CHANGE  to 1 for debugging********************/
int DEBUG=0;
/**********************************************************/

int twopop;     /* conditional index to genrate one/two MT populations */
float popdist;  /* distance between populations' CM */
float t_length; /*Tubulin length*/
enum boolean {FALSE, TRUE} manual,FixedDensity,lengthmix,TotFluxConst, MOTOR, D3_output, ACTIN
  ,AvNeighboredMTs;

/* Model for adding new MTs (Difference between LB LBPOP:
 LB makes new MTs show up with their tip first at the LB
 LBPOP however lets MTs popup at the left bound with full length ) */
enum model {LB,LBPOP,RB,RLB,CENTER,ALL} addmodel;
enum polaritymodel {FIXED,LOCAL,GC,GCLOCAL} polaritymodel;


struct point {
  float x,y,z;
  };

/* Input/Output */
char *dirname; /* directory name for input/output files */

/* Seed for Random Number Generator */
/* long idum; */

/* Microtubul definitions */
int MTfirst,newMTs0,newfreq;
struct MTdata {
  /* MT.cm -- center of mass */
  /* MT.rcm -- real position of cm; used in case of periodic b.c. */
  /* MT.pend -- plus end of MT */
  /* MT.mend -- minus end of MT */
  /* MT.direct -- direction of MT; unit vector */
  /* MT.vel -- velocity */
  /* MT.rxmin/MT.rxmax -- real positions of xmax and xmin; used in periodic
                          b.c         */
  /* MT.iter0 -- is the iteration number at which that MT entered the system */

  struct point *cm,*rcm,*pend, *mend, *direct, *vel;
  int *grow_direct;//1 for growing -1 for shrinking 0 for idle
  float *length;
  float maxL;
  float minL;
  int number;
  int *is_in_net;
  int *index,*iter0, *cluster;
  int imin,imax,izmin,izmax,iymin,iymax,rimin,rimax;
    int imin_real,imax_real;
    float xmin_real, xmax_real;
  float xmin,xmax,zmin,zmax,ymin,ymax;
  float rxmin,rxmax;
  float radi; /* hard-core radius of MT */
} MT;

struct segdata {
  struct point *direct,*cm,*mend,*pend;
  float *length;
  int number;
} seg;

struct sumdata {
  float *xr,*xr2,*xl,*xl2,*nr,*nl;
} sum;

struct Trackdata {
  int *nMT;
  int **index;
  int *index_count;/*this integer array is set to zero at the beginning
                            and then switched to 1 if index is used. SO it takes notice of used
                            indexes*/
  float **vel;
  float **xcm;
  float **ycm;
  float **zcm;
  float **length;
    float **xcm_vis; //all for the visualization
    float **ycm_vis;
    float **zcm_vis;

    float **mathematica_arrows; // Table to print out mathematica arrows see update.c and write_output.c
    float *arrow_number; // For better output: Number of Motors per step
    float *lbound; //boundaries for visuals
    float *rbound;
  float **direct;
  int **it0;
    int **cross_type;
} track;
long corr_ti,corr_tf; /* initial/final time for printing
			  time correlation function */

struct LNET { float association; float dissociation;  } Lnet;

struct BOUNDARY {
  /* boundary 0 is the left boundary */
  /* boundary 1 is the right boundary */
  enum btype {NON,ABSORB,REFLECT,PERIODIC,RFIXREF,LFIXREF,POPUP,FORCE,WALL,STICKY, NET, CONFLX} type[2];
} bound;
float reflectb0; /* position of reflcting boundary */

int wfreq;



float PolarityRatio0; /* The fraction of Left oriented and
                        Right oriented minus ends in first pull of MTs*/
float PolarityRatioNew; /* The fraction of Left oriented and
                        Right oriented minus ends in newly added MTs*/
float ProbBipolar; /* Fraction of bipolar motors -- at this point bipolar
		      motors move to plus-end like myosin, kinesin */
float ProbBundling; /* Fraction of overalping regions with Bundling
		       proteins */
float ProbLegDown;     /* Out of all Unipolar overlaping regions, this is the
                          Probability for Motor legs to point downwards */
float ProbActive; /* Probability that the Motors in overlap region are active,
                     or fraction of active motors */
float ProbKinesin;

float pFlp;

int neqsteps; /* Number of equilibration steps */
float MinOvlpDist;   /* Cutoff for overalpping distance between MTs */
float DomainX, DomainZ; /* Domain size */
float fstall_k, fstall_d, fstall_bp, vel0_k, vel0_d, vel0_bp;   /* stall force and maximum velocity */
float xim0,xisol,sticky_drag;/* friction by motors/solution/membrane, respectively */
float lamb_dyn; /* line density of minus end direct motors - numbeqr of motors per unit length */
float *dist_dyn;
float lamb_kin; /* line density of plus end direct motors - number of motors per unit length */
float *dist_kin;
float *ProbBi_dist;
float *ProbKi_dist;
float *ProbDy_dist;
float *ProbAct_dist;
float vel_flux;
int iter,niter;  /* iteration index */
int iav,nav; /* ensemble averaging iterations */
/* left and right boundaries, for first MTs, for new MTs */
float lbound0,rbound0,lboundNew,rboundNew;
float dt; /* time step */
float exclude; /* 2*exclude is the typical cross-bridge distance */
int nbox,nbox2; /* number of grid boxes in density profiles */
/* minimum and maximum values for velocity histograms            */
float hist_vel_min,hist_vel_max,hist_absvel_min,hist_absvel_max;
float randav,randstd; /* average and std of random force */

int NmaxMTs; /*Maximum number of MTs per length of MT in bundle*/
float velpol; /* Polarization velocity -- so far used only for presentation
                 of velocity histogram! No polarization yet implemented */
float veldepol;
float ext_fr; // The force applied by the right boundary
float ext_fl; // The force applied by the left boundary
float touch_depth; //The length of penetration for FORCE
float spring_r;
float spring_l;
int force_it0;
/*#define ZERO 0*/
/*#define LEGUP 1*/
/*#define LEGDOWN 2*/
struct OVLPdata {
    enum motortype {ZERO, LEGUP,LEGDOWN,BIPOLAR,TWOMOTORS,BUNDLING} **type;
    int **motor_direction;// +1 for motor walking towards - end and -1 for plus end
  float size;
    int *actin;
    int *actin_it0;
  float **iter0;
} ovlp;
struct ACTINdata{
    float attach, deattach;
} actin;



struct average {
  float **nr, **nl;
  float **rnr, **rnl;
  float **rmin, **rmax;
  float *ntot, *xmin, *xmax, *rxmin, *rxmax;
  float *zmin, *zmax, *ymin, *ymax;
  float *xr, *xr2, *xl, *xl2;
  float *xt, *xt2;
  float *rxr,*rxr2, *rxl, *rxl2;
  float *rxt, *rxt2;
  float *velt, *vell, *velr;
  float *velt2 ,*vell2 , *velr2 ;
  float *neighbors ;
  float *ovlp ;
  /* quantities for velocity ditributions of right and left handed MTs */
  float *absvmin ,*absvmax ;
  float *vmin ,*vmax ;
  float **absvel_nr, **absvel_nl, *absvel_ntot ;
  float **vel_nr, **vel_nl, *vel_ntot ;
  /* average velocity along x-axis of right and left handed MTs */
  float **xvelr , **xvell;
  float **polarcor;
  float **polarcor_ryz;
  float **velcor;       //velocity correlation
  float **thick;
  float *ovlpsize_perMT;
  float *ovlpsize_perOV;
  float **velcort,**velcort_pp,**velcort_pm,**velcort_mm, **tubulin_dens;
  float **ovlp_kin, **ovlp_dyn, **ovlp_act, **ovlp_bi;
  float *xmax_real,*xmin_real;
  float **tubulin_r, **tubulin_l;
  float **nv_cluster;/*average number of MT moving with certain velocity*/
  /* av.zeros is the average number of MTs with nzero neighbors */
  float *nMT , *nMTr , *nMTl , *nzeros ;
  float *ryz ;
  float *op ;
  /*float x2tau[MXITER],xr2tau[MXITER],xl2tau[MXITER];*/
  float *xltau,*xl2tau,*xrtau,*xr2tau,*xtau,*x2tau,*nrtau,*nltau;
} av;

struct BundleProp {
  float rad,length,density;
  float xi0,xi;
} bundle;

float maxvel;
long idum,idum2;
float kBT;
int lin_density;
int forcepol;
int forcelength;
int MCfreq; /* frequency of MC compaction */


/* --------------------------------------------------------------------- */


/* Function declaration */                             /* filename */
int getjMT(int i,int it2);                             /* getjMT.c */
int iabsVelmin();                                      /* calcav.c     */
int iabsVelmax();                                      /* calcav.c     */
int iVelmin();                                         /* calcav.c     */
int iVelmax();                                         /* calcav.c     */
int iMTmin();                                          /* calcav.c     */
int iMTmax();                                          /* calcav.c     */
void usage_err(char *argv[]);                          /* functions.c  */
void rdata(int argc, char *argv[]);                    /* functions.c  */
char *filept(char *filename);                          /* functions.c  */
char *filept_mathematica(char *filename);              /* functions.c  */
char *filept_max(char *filename);              /* functions.c  */
void flipMTs();						/*MTflip.c*/
void init();                                           /* init.c       */
void overlapMap();                                     /* overlapMap.c */
void percolation();                                     /*percolation.c*/
void mkMTs(int newMTs);                                /* mkMTs.c      */
void mkseg(int iMT);                                   /* boundary.c   */
void boundary();                                   /* boundary.c   */
void nrerror(char error_text[]);                    /* nrutil.c*/
double overlap(int ia, int ib);                         /* overlap.c    */
double xoverlap(int ia, int ib);                 /* overlap.c    */
double soverlap(int ia, int ib);                        /* overlap.c    */
double sxoverlap(double xa, double xb,
		 double lena, double lenb);            /* overlap.c    */
void motorvec(int ia, int ib, struct point *motorhead,
	      struct point *motorleg);                 /* overlap.c    */
float maxf(float a, float b);                          /* functions.c  */
float maxi(int a, int b);                              /* functions.c  */
float minf(float a, float b);                          /* functions.c  */
float mini(int a, int b);                              /* functions.c  */
void fmotor(int iMT, struct point *fm);                /* force.c      */
void fmdrag(int iMT, struct point *fd);                /* force.c      */
void fsol(int iMT, struct point *fd);                  /* force.c      */
void ffrand(int iMT, struct point *fr);                 /* force.c      */
void equations();                                      /* equations.c  */
void ludcmp(float **a, int n, int *indx, float *d);    /* ludcmp.c     */
void lubksb(float **a, int n, int *indx, float b[]);   /* lubksb.c     */
void mprove(float **a, float **alud, int n, int indx[],
	    float b[], float x[]);                     /* mprove.c     */
float min(float A, float B);
float max(float A, float B);
void add_bounds(int iMT,float *rb, float *lb);
int get_mathematica_arrow2D();
int get_mathematica_arrow3D();
float gasdev(long *idum);
void wdata();
void wav_density (int filenum);
void wav_rdensity (int filenum);
void wav_MTorder(int filenum);
void wav_nzeros();
void wav_OrderParam();
void wav_absveldensity (int filenum);
void wav_veldensity (int filenum);
void wav_polarcor (int filenum);
void wav_polarcor_ryz (int filenum);
void wav_velcor (int filenum);
void wav_x();
void wav_x2();
void wav_x2tau();
void wav_Rg();
void wav_V2();
void wav_V();
void wav_neighbors();
void wav_ovlp();
void wav_end2end();
void wav_RgYZ();
void wav_stallforce ();
void wav_absvel();
void wav_percprob ();
void wav_percprob_bi ();
void wav_avovlpsize_perMT();
void wav_avovlpsize_perOV();
void wav_avstrands_tot();
void wav_avstrands(int filenum);
void avcm2(struct point *av);
void avvel2(struct point *av);
int iMTcmmax();
int iMTcmmin();
int iMTycmmax();
int iMTycmmin();
int iMTzcmmax();
int iMTzcmmin();
int iMTxmax();
int iMTxmin();
int riMTcmmax();
int riMTcmmin();
int riMTxmax();
int riMTxmin();
int nMTright(float x, float dx);
int nMTleft(float x, float dx);
void wdensity(int filenum, int nbox);
void wmap_xz(int filenum);
void wmap_yz(int filenum);
void wvelmap(int filenum);
void wbundle_shape (int filenum);
void wav_totabsveldensity ();
void wav_totveldensity ();
void wav_nMT();
void wav_motordens(int filenum);
void write_force_fstall();
void wav_stalllength_extf ();
void write_force_polarity();
void write_force_motorfraction();
void normalize();
void unnormalize();
void update_av();
void newpos();
char *itoa(int i);
void compact();
double vdw(struct point *coor, float len[]);
void timecorr();
void wav_timevelcor (int filenum);
int segindex(int iseg);
void wav_xtraject();
void wav_xvelprof (int filenum);
void wav_xvelprofToT();
void wav_traject_2Dmathematica();
void wav_traject_3Dmathematica();
void empty(float *array, int l);
void emptyint(int *array, int l, int r);
void m_empty (float **A, int l,int r);
void m_emptyint (int **A, int l,int r);
float highestvalue(float *array, int size);
float Interval_overlap(float Pos1, float Pos2, float length1, float length2);
struct point *pvector(long nl, long nh);
void free_pvector(struct point *v, long nl, long nh);
void Initial_dist();
int special_plots();                          /* special_plots.c */
void initialize_av_variables();
void initialize_it_variables();
void free_MTdat();
void free_all();
void grow();
void actin_overlapMap();
void wav_tubulin_dens();
void wav_tubulin_dens_it( int filenum);
void wav_tubulin_op();
void wav_clustervel();
void prob_update();
float bundle_rad();

void write_force_spring_r() ;
/* Functions for writing XMGR files */
void agr_string(FILE *fp, int gn, float xcm, float zcm, int num);
void agr_std(FILE *fp, int gn, const char *Title, const char *Xaxis, const char *Yaxis, float xmin, float ymin, float xmax, float ymax);
void agr_begset(FILE *fp, int gn, int k, int lncolor, int lnstyle, float lnwidth, int symbol, float symbsz, int symbcl);
void agr_endset(FILE *fp);
void agr_comment(FILE *fp, const char *comment);
void agr_wcoor(FILE *fp, float x, float y);
void agr_wcoor3D(FILE *fp, float x, float y, float z);
void print_line(FILE *output);
void short_line(FILE *output);
void agr_legend(FILE *fp, int k, const char *name);
void agr_cyl(FILE *fp,float xcm,float zcm,float r,int lw,int lc,int fpt,int fc);void agr_arrowline(FILE *fp, int gn, float lw, int lc, float pendx, float pendz, float mendx, float mendz, int arrow);
void agr_box(FILE *fp,float xcm,float zcm,float MTlen, float MTwidth,
	     float lw,int lc,int fpt,int fc);

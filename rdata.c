/* #################################################################### */
/* Open Input file, Extract data, and Echo it to a log file             */
#include <stdio.h>
#include "global_var.h"
#include "ran1.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>


void rdata(int argc, char *argv[])
{

  FILE *fp,*outdat;
  char string[60],str[10];
  int i,iMT,jMT,nrandom;
  long int idum3,idum4;
  float temp;
  float avMTlength;

  /* Open input data file and output data file */
  if ((fp=fopen(filept("data"),"r"))==NULL)
    {
      printf("COULDNT OPEN INPUT DATA FILE data \n");
      exit(0);
    }

  if ((outdat=fopen(filept("data.dat"),"w"))==NULL)
    {
      printf("COULDNT OPEN OUTPUT DATA FILE data.dat \n");
      exit(0);
    }

  /* defalt data gets overwritten later */
  MT.maxL=10.0*microM; /* micrometer */
  MT.minL=10.0*microM;
  rbound0=2.0*MT.maxL;
  lbound0=0.0;
    reflectb0=0.0;
    twopop=0;
    popdist=0.0;
    touch_depth=0;
    t_length=0.008; /*Tubulin length 8nm*/
  bundle.rad=0.01; /* Bundle radius MicroM */
    bundle.xi0=0.2;
    PolarityRatio0=0.5;
  PolarityRatioNew=0.5;
  ProbLegDown=0.5;
    /*initial calls to nrand*/
    nrandom=1;
  xim0=0.0; /* pN/(microM/sec)/motor */
  xisol=0.01; /* pN/(microM/sec) */
  lamb_dyn=1.0; /* motors/microM */
  lamb_kin=1.0; /* motors/microM */
  dt=0.2;   /* sec */
  fstall_k=1; /* pN/motor */
  vel0_k=0.5; /* microM/sec */
  niter=100; /* 1hour -- each iter 1 sec */
  exclude=0.07; /* excluded width of a MT */
  nbox=300;
  nav=1;
  temp=30;
    wfreq=500;
    MCfreq=0;
  ext_fr=0.0;
  ext_fl=0.0;
  ProbActive=1.0;
  ProbBipolar=0.0;
    MinOvlpDist=0.154;
  TotFluxConst=0;
    corr_tf=2;
    corr_ti=1;
  FixedDensity=FALSE;
  neqsteps=500; /* number of equilibration steps */

  idum=-1; /* initialize random number generator */
  ran1(&idum);
  /*MT.radi=0.015;*/ /* 15nm hard-core radius of MT */
  MT.radi=0.0;
  /* printf("rdata:> minimal interaction distance set to zero!!!\n"); */
/*   getchar(); */
  lengthmix=FALSE; /* TRUE, implies a mixter of 2 MT lengths */
  AvNeighboredMTs=TRUE; /* Ensemble averaging will include only MTs
			   that have at least one neighbor MT */

  /* printf("-----------\n\n"); */
  /*   printf("idum3=%d:>\n",idum3); */
  /*   for (i=1;i<10;i++) */
  /*     { */
  /*       float r; */
  /*       r=ran1(&idum3); */
  /*       printf("r=%f  idum3=%d\n",r,idum3); */
  /*     } */

  /*   idum4=-1;ran1(&idum4);idum4=11;  */
  /*   printf("-----------\n\n"); */
  /*   printf("idum4=%d:>\n",idum4); */
  /*   for (i=1;i<10;i++) */
  /*     { */
  /*       float r; */
  /*       r=ran1(&idum4); */
  /*       printf("r=%f  idum4=%d\n",r,idum4); */
  /*     } */
  /*   getchar(); */



  /* Read data from data file */

  /* Boundary Conditions */
  /* Left boundary */
  fscanf(fp,"%s %s\n",&string,&str);
  if (strcmp(str,"NON")==0)
    {
      bound.type[0]=NON;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"ABSORB")==0)
    {
      bound.type[0]=ABSORB;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"POPUP")==0)
    {
      bound.type[0]=POPUP;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"REFLECT")==0)
    {
      /* not implemented */
      bound.type[0]=REFLECT;
      fprintf(outdat,"%s %s\n",string,str);
      fprintf(outdat,"rdata:> bound.type[0]=%d;\n",bound.type[0]);
      fprintf(outdat,"rdata:> Reflecting left bound not implemented!\n");
      exit(0);
    }
  else if(strcmp(str,"PERIODIC")==0)
    {
      bound.type[0]=PERIODIC;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"LFIXREF")==0)
    {
      /* not implemented */
      bound.type[0]=LFIXREF;
      fprintf(outdat,"%s %s\n",string,str);
      fprintf(outdat,"rdata:> bound.type[0]=%d;\n",bound.type[0]);
      fprintf(outdat,"rdata:> Fixed Reflecting left bound not implemented!\n");
    }
  else if(strcmp(str,"WALL")==0){
      bound.type[0]=WALL;
      fprintf(outdat,"%s %s\n",string,str);
  }
  else if(strcmp(str,"STICKY")==0){
      bound.type[0]=STICKY;
      fprintf(outdat,"%s %s\n",string,str);
  }
  else if(strcmp(str,"FORCE")==0){
      bound.type[0]=FORCE;
      fprintf(outdat,"%s %s\n",string,str);
  }
  else if(strcmp(str,"NET")==0){
      bound.type[0]=NET;
      fprintf(outdat,"%s %s\n",string,str);
  }
  else if(strcmp(str,"CONFLX")==0){
      bound.type[0]=CONFLX;
      fprintf(outdat,"%s %s\n",string,str);
  }

  else{
      printf("Chosen left bound is not valid");
      exit(0);
  }


  /* Right boundary */
  fscanf(fp,"%s %s\n",&string,&str);
  if (strcmp(str,"NON")==0)
    {
      bound.type[1]=NON;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"ABSORB")==0)
    {
      bound.type[1]=ABSORB;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"POPUP")==0)
    {
      bound.type[1]=POPUP;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"REFLECT")==0)
    {
      bound.type[1]=REFLECT;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"PERIODIC")==0)
    {
      bound.type[1]=PERIODIC;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"RFIXREF")==0)
    {
      bound.type[1]=RFIXREF;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"FORCE")==0){
      bound.type[1]=FORCE;
      fprintf(outdat,"%s %s\n",string,str);
  }

  else{
      printf("Chosen right bound is not valid\n");
      exit(0);
  }

  fscanf(fp,"%s %s\n",&string,&str);
  if (strcmp(str,"RB")==0)
    {
      addmodel=RB;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"LB")==0)
    {
      addmodel=LB;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"LBPOP")==0)
  {
      addmodel=LBPOP;
      fprintf(outdat,"%s %s\n",string,str);
  }
  else if(strcmp(str,"ALL")==0)
    {
      addmodel=ALL;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"RLB")==0)
    {
      addmodel=RLB;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"CENTER")==0)
    {
      addmodel=CENTER;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else
    {
      printf("rdata:> Erorr: addmodel not selected!\a\n");
      exit(0);
    }
//POLARITY MODEL
    fscanf(fp,"%s %s\n",&string,&str);
    if (strcmp(str,"FIXED")==0)
      {
        polaritymodel=FIXED;
      }
    else if(strcmp(str,"LOCAL")==0)
      {
        polaritymodel=LOCAL;
      }
    else if(strcmp(str,"GC")==0)
        {
          polaritymodel=GC;
        }
    else if(strcmp(str,"GCLOCAL")==0)
        {
          polaritymodel=GCLOCAL;
        }
    else{
      printf("rdata:> Erorr: polaritymodel not selected!\a\n");
      exit(0);
    }
    fprintf(outdat,"%s %s\n",string,str);

  fscanf(fp,"%s %s\n",&string,&str);
  if (strcmp(str,"TRUE")==0)
    {
      TotFluxConst=TRUE;
      fprintf(outdat,"%s %s\n",string,str);
    }
  else if(strcmp(str,"FALSE")==0)
    {
      TotFluxConst=FALSE;
      fprintf(outdat,"%s %s\n",string,str);
    }


  fscanf(fp,"%s %s\n",&string,&str);

  if (strcmp(str,"TRUE")==0)
  {
      D3_output=TRUE;
      //printf("TRUE\n");
      fprintf(outdat,"%s %s\n",string,str);
  }
  else if(strcmp(str,"FALSE")==0)
  {
      D3_output=FALSE;
     // printf("FALSE\n");
      fprintf(outdat,"%s %s\n",string,str);
  }
  else {
      printf("cannot read data file 3D visuals WHATTHEFUCK?\nSETTING DEFAULT: 3D=TRUE\n");
      D3_output=TRUE;
  }

    fscanf(fp,"%s %s\n",&string,&str);

    if (strcmp(str,"TRUE")==0)
    {
        MOTOR=TRUE;
        //printf("TRUE\n");
        fprintf(outdat,"%s %s\n",string,str);
    }
    else if(strcmp(str,"FALSE")==0)
    {
        MOTOR=FALSE;
        // printf("FALSE\n");
        fprintf(outdat,"%s %s\n",string,str);
    }
    else {
        printf("cannot read data file MOTOR WHATTHEFUCK?\nSETTING DEFAULT: 3D=TRUE\n");
        exit(0);
    }

    fscanf(fp,"%s %s\n",&string,&str);
    if (strcmp(str,"TRUE")==0)
    {
        ACTIN=TRUE;
        //printf("TRUE\n");
        fprintf(outdat,"%s %s\n",string,str);
    }
    else if(strcmp(str,"FALSE")==0)
    {
        ACTIN=FALSE;
        // printf("FALSE\n");
        fprintf(outdat,"%s %s\n",string,str);
    }
    else {
        printf("cannot read data file ACTIN WHATTHEFUCK?\n");
        exit(0);
    }
    fscanf(fp,"%s %f\n",&string,&dt);
    fprintf(outdat,"%s %f\n",string,dt);

    fscanf(fp,"%s %f\n",&string,&ext_fr);
    fprintf(outdat,"%s %f\n",string,ext_fr);
    fscanf(fp,"%s %f\n",&string,&ext_fl);
    fprintf(outdat,"%s %f\n",string,ext_fl);
  fscanf(fp,"%s %d\n",&string,&niter);
  fprintf(outdat,"%s %d\n",string,niter);
  fscanf(fp,"%s %d\n",&string,&neqsteps);
  fprintf(outdat,"%s %d\n",string,neqsteps);

  fscanf(fp,"%s %d\n",&string,&nav);
  fprintf(outdat,"%s %d\n",string,nav);


  fscanf(fp,"%s %d\n",&string,&MTfirst);
  fprintf(outdat,"%s %d\n",string,MTfirst);
  fscanf(fp,"%s %d\n",&string,&newfreq);
  fprintf(outdat,"%s %d\n",string,newfreq);
  fscanf(fp,"%s %d\n",&string,&newMTs0);
  fprintf(outdat,"%s %d\n",string,newMTs0);
  fscanf(fp,"%s %f\n",&string,&lbound0);
  fprintf(outdat,"%s %f\n",string,lbound0);
  fscanf(fp,"%s %f\n",&string,&rbound0);
  fprintf(outdat,"%s %f\n",string,rbound0);

  fscanf(fp,"%s %f\n",&string,&MT.maxL);
  fprintf(outdat,"%s %f\n",string,MT.maxL);
  fscanf(fp,"%s %f\n",&string,&MT.minL);
  fprintf(outdat,"%s %f\n",string,MT.minL);

  fscanf(fp,"%s %f\n",&string,&PolarityRatio0);
  fprintf(outdat,"%s %f\n",string,PolarityRatio0);
  fscanf(fp,"%s %f\n",&string,&PolarityRatioNew);
  fprintf(outdat,"%s %f\n",string,PolarityRatioNew);
  fscanf(fp,"%s %f\n",&string,&ProbBipolar);
  fprintf(outdat,"%s %f\n",string,ProbBipolar);

  fscanf(fp,"%s %f\n",&string,&ProbBundling);
  fprintf(outdat,"%s %f\n",string,ProbBundling);
  fscanf(fp,"%s %f\n",&string,&ProbLegDown);
  fprintf(outdat,"%s %f\n",string,ProbLegDown);
    fscanf(fp,"%s %f\n",&string,&ProbKinesin);
    fprintf(outdat,"%s %f\n",string,ProbKinesin);
  fscanf(fp,"%s %f\n",&string,&ProbActive);
  fprintf(outdat,"%s %f\n",string,ProbActive);
  fscanf(fp,"%s %f\n",&string,&fstall_k);
  fprintf(outdat,"%s %f\n",string,fstall_k);
    fscanf(fp,"%s %f\n",&string,&fstall_d);
    fprintf(outdat,"%s %f\n",string,fstall_d);
    fscanf(fp,"%s %f\n",&string,&fstall_bp);
    fprintf(outdat,"%s %f\n",string,fstall_bp);
    fscanf(fp,"%s %f\n",&string,&lamb_dyn);
    fprintf(outdat,"%s %f\n",string,lamb_dyn);
    fscanf(fp,"%s %f\n",&string,&lamb_kin);
    fprintf(outdat,"%s %f\n",string,lamb_kin);

  fscanf(fp,"%s %f\n",&string,&vel0_k);
    fprintf(outdat,"%s %f\n",string,vel0_k);
    fscanf(fp,"%s %f\n",&string,&vel0_d);
    fprintf(outdat,"%s %f\n",string,vel0_d);
    fscanf(fp,"%s %f\n",&string,&vel0_bp);
    fprintf(outdat,"%s %f\n",string,vel0_bp);

    fscanf(fp,"%s %f\n",&string,&vel_flux);
    fprintf(outdat,"%s %f\n",string,vel_flux);
  fscanf(fp,"%s %f\n",&string,&xisol);
  fprintf(outdat,"%s %f\n",string,xisol);

    fscanf(fp,"%s %f\n",&string,&sticky_drag);
    fprintf(outdat,"%s %f\n",string,sticky_drag);

  fscanf(fp,"%s %d\n",&string,&nbox);
  fprintf(outdat,"%s %d\n",string,nbox);



    fscanf(fp,"%s %f\n",&string,&spring_r);
    fprintf(outdat,"%s %f\n",string,spring_r);
    fscanf(fp,"%s %f\n",&string,&spring_l);
    fprintf(outdat,"%s %f\n",string,spring_l);
    fscanf(fp,"%s %i\n",&string,&force_it0);
    fprintf(outdat,"%s %i\n",string,force_it0);
    fscanf(fp,"%s %f\n",&string,&Lnet.association);
    fprintf(outdat,"%s %f\n",string,Lnet.association);

    fscanf(fp,"%s %f\n",&string,&Lnet.dissociation);
    fprintf(outdat,"%s %f\n",string,Lnet.dissociation);
    fscanf(fp,"%s %f\n",&string,&actin.attach);
    fprintf(outdat,"%s %f\n",string,actin.attach);
    fscanf(fp,"%s %f\n",&string,&actin.deattach);
    fprintf(outdat,"%s %f\n",string,actin.deattach);
    fscanf(fp,"%s %f\n",&string,&pFlp);
    fprintf(outdat,"%s %f\n",string,pFlp);

    fscanf(fp,"%s %f\n",&string,&degprob);
    fprintf(outdat,"%s %f\n",string,degprob);

  if (spring_r>0){
    ext_fr=0.;
  }
  if (spring_l>0){
    ext_fl=0.;
}
  /* check consistancy of variables */
  if ((bound.type[0]==PERIODIC && bound.type[1]!=PERIODIC) ||
      (bound.type[1]==PERIODIC && bound.type[0]!=PERIODIC))
    {
      fprintf(stdout,"Error in radata:\n");
      fprintf(stdout,"only one of two boundaries was asigned PERIODIC:\n");
      exit(0);
    }
  if (bound.type[1]==PERIODIC && bound.type[1]==PERIODIC &&
      MT.maxL>=(rbound0-lbound0))
    {
      fprintf(stdout,"Error in radata:\n");
      fprintf(stdout,"Using periodic boundary cond.\n");
      fprintf(stdout,"MT length must be smaller\n");
      fprintf(stdout,"than domain size rbound0-lbound0\n");
      fprintf(stdout,"MaxL=%f rbound0=%f, lbound0=%f\n",
	      MT.maxL,rbound0,lbound0);
      exit(0);
    }




  if ((int)fmod(bundle.rad,2*exclude)!=0)
    {
      fprintf(stdout,"Error in rdata:\n");
      fprintf(stdout,"Bundle radius should be integral number\n");
      fprintf(stdout,"of excluded volume diameter\n");
      fprintf(stdout,"bundle.rad=%f 2exclude=%f\n",bundle.rad,2*exclude);
      fprintf(stdout,"bundle.rad % 2exclude=%f\n",fmod(bundle.rad,2*exclude),2*exclude);
      exit(0);
    }
    if (MT.maxL >= (rbound0-lbound0-1.1) || MT.maxL >= (rbound0-lbound0-1.1) ) {
        printf("ERROR: maximal MT length larger than system\n");
        exit(0);
    }
    if (bound.type[0]==CONFLX && newMTs0==0 && TotFluxConst==FALSE) {
        printf("ERROR: b.c. CONFLX needs MT influx\n");
        exit(0);
    }
    if (corr_ti==0 || corr_tf==corr_ti) {
        printf("Invalid correlation time initial/final step (corr_ti corr_tf)\n");
        exit(0);
    }
    if (ProbBipolar+ProbKinesin>1){
        printf("all motor probabilities sum up to larger than 1!\n");
        exit(0);
    }


  /* Calculate quantities from data */

  kBT=kB*temp;
    if (newMTs0==0) {
        maxvel= 4.;
        //maxvel= ((float)MTfirst-1.)/2.*vel0_k;
    }
    else{
        printf("Warning: Max number of MTs not defined! => maximum velocity relies on educated guess...\n");

        maxvel =4.;
        //maxvel= ((float)niter/(float)newfreq*(float)newMTs0-1.)/2.*vel0_k;
    }
 //calculate probabilities from data
 //Lnet.dissociation=1-exp(-1*Lnet.dissociation*dt);
 Lnet.dissociation=Lnet.dissociation*dt;

//pFlp=1-exp(-1*pFlp*dt);
 pFlp=pFlp*dt;
 //Lnet.association=1-exp(-1*Lnet.association*dt);
 Lnet.association=Lnet.association*dt;
 newfreq=(int)round((float)newfreq/dt);
 force_it0=(int)round((float)force_it0/dt);
if ((bound.type[0]==NET && (Lnet.dissociation>1 || Lnet.association>1)) || pFlp>1){
 	printf("please choose smaller rates (>1)!");
 	exit(0);
 }

printf("Lnetdiss= %f \n Lnetasso= %f \n pFl=  %f \n newfreq= %i \n",Lnet.dissociation, Lnet.association, pFlp, newfreq);

 if (MTfirst<=0)
    MTfirst=(int)floor(0.5+bundle.xi0*(rbound0-lbound0)*
		       (3*(0.5*bundle.rad/exclude+1)*
			(0.5*bundle.rad/exclude)+1)/
		       ((MT.maxL+MT.minL)/2.0));

  fprintf(outdat,"------------------------------\n");
  fprintf(outdat,"Diffution constant for a 1 microM MT:\nD=%e (microM)^2/sec = %e cm^2/sec\n",kBT/xisol,kBT/xisol*1E-8);


  for(i=1;i<=nrandom;i++)
    ran1(&idum);

  if (AvNeighboredMTs)
    {
      printf("rdata:>Ensemble Averaging will include only MTs that have\n       One neighbor at least!\n");
      fprintf(outdat,"rdata:>Ensemble Averaging will include only MTs that have\n       One neighbor at least!\n");
    }

  if (ACTIN==TRUE)
      printf("Motor overlaps are also possible to the actin myosin membrane\n");

  if ((bound.type[1]==REFLECT || bound.type[1]==RFIXREF) &&
      ProbBipolar>0.0)
    {
      printf("rdata:> Bipolar motors are not yet implemented with\nReflective boundary condition!\n");
      printf("rdata:> ProbBipolar=%f\n",ProbBipolar);
      exit(0);
    }

  if (neqsteps>=niter)
    {
      printf("Number of equilibration steps exceed total number of steps");
      printf("neqsteps=%d   nstep=%d",neqsteps,niter);
      exit(0);
    }

  if (corr_tf<corr_ti || corr_tf>niter)
    {
      printf("Final correlation time is smaller than Initial value");
      printf("or, final corr-time is larger than niter");
      printf("corr_ti=%d    corr_tf=%d    niter=%d",corr_ti,corr_tf,niter);
      exit(0);
    }


  /* minimum and maximum values for velocity histograms */


  hist_vel_min = -20*vel0_d;
  hist_vel_max = 20*vel0_d;
  hist_absvel_min = 0.0;
  hist_absvel_min = 20*vel0_d;


  /* nbox2 - alternative number of grid boxes - typically smaller division */
  nbox2=40;

/* Update input according to command line arguments. */

  fclose(fp);
  fclose(outdat);

    if (ProbKinesin<1.) {
        printf("ASSAF OUTPUT WILL NORMALIZE AV VALUES WITH RESPECT TO DYNEIN VELOCITY\n");
    }


}

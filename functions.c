#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "global_var.h"
#include <math.h>

#define NR_END 1
#define FREE_ARG char*


void empty(float *array, int l)
{
    int i;
    for (i=0 ; i<=l ; i++)
    {
        array[i]=0;
    }
    
}

void emptyint(int *array, int l, int r)
{
    int i;
    for (i=l ; i<=r ; i++)
    {
        array[i]=0;
    }
    
}

void m_empty (float **A, int l,int r){
    int i,j;
    for (i=l; i<=r; i++) {
        for ( j=l; j<=r; j++) {
            A[i][j]=0.0;
        }
    }
}

void m_emptyint (int **A, int l,int r){
    int i,j;
    for (i=l; i<=r; i++) {
        for ( j=l; j<=r; j++) {
            A[i][j]=0;
        }
    }
}


float highestvalue(float *array, int size) {
    float p=array[0];
    int i;
    for (i=1; i<size; i++) {
        if (p<array[i]) p=array[i];
    }
    return(p);
}

float lowestvalue(float *array, int size) {
    float p=array[0];
    int i;
    for ( i=1; i<size; i++) {
        if (p>array[i]) p=array[i];
    }
    return(p);
}


float bundle_rad(){
    int iMT;
    float rad, rad_max=0;
    for (iMT=1; iMT<=MT.number; iMT++) {
        rad=sqrtf(MT.cm[iMT].y*MT.cm[iMT].y+MT.cm[iMT].z*MT.cm[iMT].z);
        if (rad_max<rad) {
            rad_max=rad;
        }
    }
    return(rad_max);
}


float max(float A, float B)
{
    return (A>B)?A:B;
}

float min(float A, float B)
{
    return (A<B)?A:B;
}

float Interval_overlap(float Pos1, float Pos2, float length1, float length2)
{
    return ((float)max(0, min((double)Pos1+(double)length1, (double)Pos2+(double)length2) - max((double)Pos1,(double)Pos2) ));
}

void usage_err(char *argv[])
{
  printf("Usage Error:>\n");
  printf("Use: %s -outdir <output dir name>\n",argv[0]);
  printf("\n+Optional arguments:\n");
    printf("[-lindist for linearily rising/falling lamda]\n");
    printf("[-forcepol for for plotting force against polarity for given parameters]\n");
        printf("[-forcelength for for plotting bundle length against time for different forces (only for const ext force)]\n");
        printf("[-MTflip followed by a flipping rate in #/sec to specify MT flip probability if desired different from 0]\n");
  exit(0);
}

char *filept(char *filename)
{
  /* This function concatonates the file name to the dirname
     given as an argument to main() */
  static char *pt;
  int length=strlen(filename)+strlen(dirname)+1;
  pt=(char *)calloc(length,sizeof(char));
  strcat(pt,dirname);
  strcat(pt,filename);
  return(pt);
}

char *filept_max(char *filename)
{
    /* This function concatonates the file name to the Max_output
     given as an argument to main() */
    
    static char *pt;
    int length=strlen(filename)+strlen(dirname)+strlen("Max_output/")+1;
    pt=(char *)calloc(length,sizeof(char));
    strcat(pt,dirname);
    strcat(pt,"Max_output/");
    strcat(pt,filename);
    return(pt);}

char *filept_mathematica(char *filename)
{
    /* This function concatonates the file name to the dirname
     given as an argument to main() */

    static char *pt;
    int length=strlen(filename)+strlen(dirname)+strlen("Mathematica/")+1;
    pt=(char *)calloc(length,sizeof(char));
    strcat(pt,dirname);
    strcat(pt,"Mathematica/");
    strcat(pt,filename);
    return(pt);
}

void Initial_dist( ){
    
    int i;
    //linear or const
    if (lin_density) {
        for ( i=1; i<=nbox; i++) {
            dist_dyn[i]=lamb_dyn;
            dist_kin[i]=lamb_kin*i/nbox;
            ProbKi_dist[i]=ProbKinesin;
            ProbDy_dist[i]=1-ProbBipolar-ProbKinesin;
            ProbBi_dist[i]=ProbBipolar;
            ProbAct_dist[i]=ProbActive;
        }
    }
    else{
        for ( i=1; i<=nbox; i++) {
            dist_dyn[i]=lamb_dyn;
            dist_kin[i]=lamb_kin;
            ProbKi_dist[i]=ProbKinesin;
            ProbDy_dist[i]=1-ProbBipolar-ProbKinesin;
            ProbBi_dist[i]=ProbBipolar;
            ProbAct_dist[i]=ProbActive;
        }
    }
}


float maxf(float a, float b)
{
  float max;
  max = (a>b) ? a : b;
  return max;
}

float maxi(int a, int b)
{
  int max;
  max = (a>b) ? a : b;
  return max;
}

float minf(float a, float b)
{
  float min;
  min = (a<b) ? a : b;
  return min;
}

float mini(int a, int b)
{
  int min;
  min = (a<b) ? a : b;
  return min;
}

char *itoa(int i)
{
  static char *string;
  if (i>-10 && i<0)
    {
      string=(char *)calloc(4,sizeof(char));
      string[0]='-';
      string[1]='0';
      string[2]='0'-i;
      string[3]='\0';
      return (string);
    } 
  if (i>=0 && i<10)
    {
      string=(char *)calloc(3,sizeof(char));
      string[0]='0';
      string[1]='0'+i;
      string[2]='\0';
      return (string);
    }
  else if (i>=10 && i<100)
    {
      string=(char *)calloc(3,sizeof(char));
      string[0]='0'+i/10;
      string[1]='0'+i%10;
      string[2]='\0';
      return (string);
    }
  else if (i>=100 && i<1000)
    {
      string=(char *)calloc(4,sizeof(char));
      string[0]='0'+i/100;
      string[1]='0'+i%100/10;
      string[2]='0'+i%100%10;
      string[3]='\0';
      return(string);
    }
  else if (i>=1000 && i<10000)
    {
      string=(char *)calloc(5,sizeof(char));
      string[0]='0'+i/1000;
      string[1]='0'+i%1000/100;
      string[2]='0'+i%1000%100/10;
      string[3]='0'+i%1000%100%10;
      string[4]='\0';
      return(string);
    }
  else if (i>=10000 && i<100000)
    {
      string=(char *)calloc(6,sizeof(char));
      string[0]='0'+i/10000;
      string[1]='0'+i%10000/1000;
      string[2]='0'+i%10000%1000/100;
      string[3]='0'+i%10000%1000%100/10;
      string[4]='0'+i%10000%1000%100%10;
      string[5]='\0';
      return(string);
    }
  else
    {
      string=(char *)calloc(7,sizeof(char));
      
      string[0]='0'+i/100000;
      string[1]='0'+i%100000/10000;
      string[2]='0'+i%100000%10000/1000;
      string[3]='0'+i%100000%10000%1000/100;
      string[4]='0'+i%100000%10000%1000%100/10;
      string[5]='0'+i%100000%10000%1000%100%10;
      string[6]='\0';
      return(string);
    }
}
 
/* ----------------------------------------------------------------- */

struct point *pvector(long nl, long nh)
{
 struct point *v;
  v=(struct point *)malloc((size_t) (nh-nl+1+NR_END)*sizeof(struct point));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl+NR_END;
} 

void free_pvector(struct point *v, long nl, long nh)
/* free a vector of struct point allocated with pvector()*/
{
  free((FREE_ARG) (v+nl-NR_END));  
}

/*-------------------------------------------------------*/

int segindex(int iseg)
{
  int iMT;

  if (iseg%2==0) 
    iMT=iseg/2;
  else
    iMT=(iseg+1)/2;

  return iMT;
}

float sign(float x){
    if (x>0)
        return(fabs(x)/x);
    else
        return(0.);
}

float arctan(float x, float y){
    if (y>0)
        return(atan (x/y));
    else
        return(PI/2);
}

int get_mathematica_arrow2D(){
    
    float sz, actin_z,actin_y;
    int jMT, iMT, ihigh, jhigh, tracker=1, typetrack=1, iseg, jseg;
    for (iseg=1; iseg<=2*MT.number; iseg++) {
        iMT=segindex(iseg);
        
        /*Actin motor connection*/
        if (ACTIN==TRUE) {
            if (ovlp.actin[iMT]==1 ) {
                actin_y= (cos(arctan( MT.cm[iMT].z,MT.cm[iMT].y ))*exclude)*sign(MT.cm[iMT].y)+MT.cm[iMT].y;
                actin_z= (sin(arctan( MT.cm[iMT].z,MT.cm[iMT].y ))*exclude)*sign(MT.cm[iMT].z)+MT.cm[iMT].z;
                track.mathematica_arrows[iter][tracker]=MT.cm[iMT].x;
                track.mathematica_arrows[iter][tracker+1]=100*(actin_z+actin_y)/sqrt(2.);
                track.mathematica_arrows[iter][tracker+2]=MT.cm[iMT].x;
                track.mathematica_arrows[iter][tracker+3]=100*(MT.cm[iMT].z+MT.cm[iMT].y)/sqrt(2.);
                tracker+=4;
                track.cross_type[iter][typetrack]=1;
                typetrack++;
            }
            else if (ovlp.actin[iMT]==-1 ) {
                actin_y= (cos(arctan( MT.cm[iMT].z,MT.cm[iMT].y ))*exclude)*sign(MT.cm[iMT].y)+MT.cm[iMT].y;
                actin_z= (sin(arctan( MT.cm[iMT].z,MT.cm[iMT].y ))*exclude)*sign(MT.cm[iMT].z)+MT.cm[iMT].z;
                track.mathematica_arrows[iter][tracker]=MT.cm[iMT].x;
                track.mathematica_arrows[iter][tracker+1]=100*(actin_z+actin_y)/sqrt(2.);
                track.mathematica_arrows[iter][tracker+2]=MT.cm[iMT].x;
                track.mathematica_arrows[iter][tracker+3]=100*(MT.cm[iMT].z+MT.cm[iMT].y)/sqrt(2.);
                tracker+=4;
                track.cross_type[iter][typetrack]=-1;
                typetrack++;
            }
        }
        
        for (jseg=iseg; jseg<=2*MT.number; jseg++) {
            jMT=segindex(jseg);
            // Only arrows to different filaments
            if (jMT != iMT) {
                //only if not zero and overlap unequal to zero
                sz=soverlap(iseg,jseg);
                if (sz>EPSI && ovlp.type[iseg][jseg]!=ZERO){
                    
                    if (ovlp.type[iseg][jseg]==BIPOLAR || ovlp.type[iseg][jseg]==BUNDLING) /*Bipolar Motors (double ended arrow)*/{
                        track.mathematica_arrows[iter][tracker]=seg.cm[iseg].x;
                        track.mathematica_arrows[iter][tracker+1]=100*(seg.cm[iseg].y+seg.cm[iseg].z)/sqrt(2.);
                        track.mathematica_arrows[iter][tracker+2]=seg.cm[jseg].x;
                        track.mathematica_arrows[iter][tracker+3]=100*(seg.cm[jseg].y+seg.cm[jseg].z)/sqrt(2.);
                        track.mathematica_arrows[iter][tracker+4]=seg.cm[jseg].x;
                        track.mathematica_arrows[iter][tracker+5]=100*(seg.cm[jseg].y+seg.cm[jseg].z)/sqrt(2.);
                        track.mathematica_arrows[iter][tracker+6]=seg.cm[iseg].x;
                        track.mathematica_arrows[iter][tracker+7]=100*(seg.cm[iseg].y+seg.cm[iseg].z)/sqrt(2.);
                        track.cross_type[iter][typetrack]=-1;
                        typetrack++;
                        track.cross_type[iter][typetrack]=-1;
                        typetrack++;
                        tracker+=8;
                    }
                    else /*UNIPOLAR*/{
                        
                        /* Determine respective orientation of Motors and MTs */
                        if (MT.cm[iMT].z>MT.cm[jMT].z
                            || ((MT.cm[iMT].z==MT.cm[jMT].z) &&
                                (MT.cm[iMT].y>MT.cm[jMT].y))
                            /* or, arrbitrary choice (iMT>jMT) if they overlap */
                            || ((MT.cm[iMT].z==MT.cm[jMT].z) &&
                                (MT.cm[iMT].y==MT.cm[jMT].y) && iseg>jseg))
                            ihigh=iMT;
                        else
                            ihigh=jMT;
                        
                        //IF LEGS on iMT and HEAD on jMT
                        if ((iMT==ihigh && ovlp.type[iseg][jseg]==LEGUP) || (jMT==ihigh && ovlp.type[iseg][jseg]==LEGDOWN)) {
                            
                            track.mathematica_arrows[iter][tracker]=seg.cm[jseg].x;
                            track.mathematica_arrows[iter][tracker+1]=100*(seg.cm[jseg].y+seg.cm[jseg].z)/sqrt(2.);
                            track.mathematica_arrows[iter][tracker+2]=seg.cm[iseg].x;
                            track.mathematica_arrows[iter][tracker+3]=100*(seg.cm[iseg].y+seg.cm[iseg].z)/sqrt(2.);
                            tracker+=4;
                            if (ovlp.motor_direction[iseg][jseg]==1){
                                track.cross_type[iter][typetrack]=1;
                                typetrack++;
                            }
                            else if (ovlp.motor_direction[iseg][jseg]==-1){
                                track.cross_type[iter][typetrack]=-1;
                                typetrack++;
                            }
                            else{
                                printf("MOTOR not specified in mathematica arrows!\n");
                                exit(0);
                            }
                        }
                        //IF LEGS on jMT and HEAD on iMT
                        else if ((jMT==ihigh && ovlp.type[iseg][jseg]==LEGUP) || (iMT==ihigh && ovlp.type[iseg][jseg]==LEGDOWN)) {
                            track.mathematica_arrows[iter][tracker]=seg.cm[iseg].x;
                            track.mathematica_arrows[iter][tracker+1]=100*(seg.cm[iseg].y+seg.cm[iseg].z)/sqrt(2.);
                            track.mathematica_arrows[iter][tracker+2]=seg.cm[jseg].x;
                            track.mathematica_arrows[iter][tracker+3]=100*(seg.cm[jseg].y+seg.cm[jseg].z)/sqrt(2.);
                            tracker+=4;
                            if (ovlp.motor_direction[iseg][jseg]==1){
                                track.cross_type[iter][typetrack]=1;
                                typetrack++;
                            }
                            else if (ovlp.motor_direction[iseg][jseg]==-1){
                                track.cross_type[iter][typetrack]=-1;
                                typetrack++;
                            }
                            else{
                                printf("MOTOR not specified in mathematica arrows!\n");
                                exit(0);
                            }
                        }
                        
                        else {
                            printf("ERRROR in get_mathematica_arrow\n");
                            if (ovlp.type[iseg][jseg]==BUNDLING) {
                                printf("CRITICAL!\n");
                            }
                        }
                    }
                    
                }
            }
        }
    }
    if (typetrack>=mxmotors){
        printf("MORE MOTORS THAN ALLOWED IN THE SYSTEM typetrack!\n");
        printf("mxmotors=%i motors=%i\n",mxmotors,typetrack);
        
        exit(0);
    }
    if (tracker>=mxmotors*12) {
        printf("MORE MOTORS THAN ALLOWED IN THE SYSTEM!\n");
        printf("mxmotors=%i motors=%i\n",mxmotors,typetrack);
        exit(0);
    }
    
    return ( (int)(((float)tracker-1)/4) );//return the number of arrows
}

int get_mathematica_arrow3D(){
    
    float sz ,actin_y,actin_z;
    int jMT, iMT, ihigh, jhigh, tracker=1, typetrack=1, iseg, jseg;
    for (iseg=1; iseg<=2*MT.number; iseg++) {
        iMT=segindex(iseg);
        /*Actin motor connection*/
        if (ACTIN==TRUE) {
            if (ovlp.actin[iMT]==1 ) {
                actin_y= (cos(arctan( MT.cm[iMT].z,MT.cm[iMT].y ))*exclude)*sign(MT.cm[iMT].y)+MT.cm[iMT].y;
                actin_z= (sin(arctan( MT.cm[iMT].z,MT.cm[iMT].y ))*exclude)*sign(MT.cm[iMT].z)+MT.cm[iMT].z;
                track.mathematica_arrows[iter][tracker]=MT.cm[iMT].x;
                track.mathematica_arrows[iter][tracker+1]=100*actin_y;
                track.mathematica_arrows[iter][tracker+2]=100*actin_z;
                track.mathematica_arrows[iter][tracker+3]=MT.cm[iMT].x;
                track.mathematica_arrows[iter][tracker+4]=100*MT.cm[iMT].y;
                track.mathematica_arrows[iter][tracker+5]=100*MT.cm[iMT].z;
                tracker+=6;
                track.cross_type[iter][typetrack]=1;
                typetrack++;
                
            }
            if (ovlp.actin[iMT]==-1 ) {
                actin_y= (cos(arctan( MT.cm[iMT].z,MT.cm[iMT].y ))*exclude)*sign(MT.cm[iMT].y)+MT.cm[iMT].y;
                actin_z= (sin(arctan( MT.cm[iMT].z,MT.cm[iMT].y ))*exclude)*sign(MT.cm[iMT].z)+MT.cm[iMT].z;
                track.mathematica_arrows[iter][tracker]=MT.cm[iMT].x;
                track.mathematica_arrows[iter][tracker+1]=100*actin_y;
                track.mathematica_arrows[iter][tracker+2]=100*actin_z;
                track.mathematica_arrows[iter][tracker+3]=MT.cm[iMT].x;
                track.mathematica_arrows[iter][tracker+4]=100*MT.cm[iMT].y;
                track.mathematica_arrows[iter][tracker+5]=100*MT.cm[iMT].z;
                tracker+=6;
                track.cross_type[iter][typetrack]=-1;
                typetrack++;
            }
            
        }

        for (jseg=iseg; jseg<=2*MT.number; jseg++) {
            jMT=segindex(jseg);
            // Only arrows to different filaments
            if (jMT != iMT) {
                //only if not zero and overlap unequal to zero
                sz=soverlap(iseg,jseg);
                if (sz>EPSI && ovlp.type[iseg][jseg]!=ZERO){
                    
                    if (ovlp.type[iseg][jseg]==BIPOLAR || ovlp.type[iseg][jseg]==BUNDLING) /*Bipolar Motors (double ended arrow)*/{
                        track.mathematica_arrows[iter][tracker]=seg.cm[iseg].x;
                        track.mathematica_arrows[iter][tracker+1]=100*seg.cm[iseg].y;
                        track.mathematica_arrows[iter][tracker+2]=100*seg.cm[iseg].z;
                        track.mathematica_arrows[iter][tracker+3]=seg.cm[jseg].x;
                        track.mathematica_arrows[iter][tracker+4]=100*seg.cm[jseg].y;
                        track.mathematica_arrows[iter][tracker+5]=100*seg.cm[jseg].z;
                        track.mathematica_arrows[iter][tracker+6]=seg.cm[jseg].x;
                        track.mathematica_arrows[iter][tracker+7]=100*seg.cm[jseg].y;
                        track.mathematica_arrows[iter][tracker+8]=100*seg.cm[jseg].z;
                        track.mathematica_arrows[iter][tracker+9]=seg.cm[iseg].x;
                        track.mathematica_arrows[iter][tracker+10]=100*seg.cm[iseg].y;
                        track.mathematica_arrows[iter][tracker+11]=100*seg.cm[iseg].z;
                        tracker+=12;
                        track.cross_type[iter][typetrack]=-1;
                        typetrack++;
                        track.cross_type[iter][typetrack]=-1;
                        typetrack++;
                    }

                    else /*UNIPOLAR*/{
                        
                        /* Determine respective orientation of Motors and MTs */
                        if (MT.cm[iMT].z>MT.cm[jMT].z
                            || ((MT.cm[iMT].z==MT.cm[jMT].z) &&
                                (MT.cm[iMT].y>MT.cm[jMT].y))
                            /* or, arrbitrary choice (iMT>jMT) if they overlap */
                            || ((MT.cm[iMT].z==MT.cm[jMT].z) &&
                                (MT.cm[iMT].y==MT.cm[jMT].y) && iseg>jseg))
                            ihigh=iMT;
                        else
                            ihigh=jMT;
                        
                        //IF LEGS on iMT and HEAD on jMT
                        if ((iMT==ihigh && ovlp.type[iseg][jseg]==LEGUP) || (jMT==ihigh && ovlp.type[iseg][jseg]==LEGDOWN)) {
                            
                            track.mathematica_arrows[iter][tracker]=seg.cm[jseg].x;
                            track.mathematica_arrows[iter][tracker+1]=100*seg.cm[jseg].y;
                            track.mathematica_arrows[iter][tracker+2]=100*seg.cm[jseg].z;
                            track.mathematica_arrows[iter][tracker+3]=seg.cm[iseg].x;
                            track.mathematica_arrows[iter][tracker+4]=100*seg.cm[iseg].y;
                            track.mathematica_arrows[iter][tracker+5]=100*seg.cm[iseg].z;
                            tracker+=6;
                            if (ovlp.motor_direction[iseg][jseg]==1){
                                track.cross_type[iter][typetrack]=1;
                                typetrack++;
                            }
                            else if (ovlp.motor_direction[iseg][jseg]==-1){
                                track.cross_type[iter][typetrack]=-1;
                                typetrack++;
                            }
                            else{
                                printf("MOTOR not specified in mathematica arrows!\n");
                                exit(0);
                            }
                        }
                        //IF LEGS on jMT and HEAD on iMT
                        else if ((jMT==ihigh && ovlp.type[iseg][jseg]==LEGUP) || (iMT==ihigh && ovlp.type[iseg][jseg]==LEGDOWN)) {
                            
                            track.mathematica_arrows[iter][tracker]=seg.cm[iseg].x;
                            track.mathematica_arrows[iter][tracker+1]=100*seg.cm[iseg].y;
                            track.mathematica_arrows[iter][tracker+2]=100*seg.cm[iseg].z;
                            track.mathematica_arrows[iter][tracker+3]=seg.cm[jseg].x;
                            track.mathematica_arrows[iter][tracker+4]=100*seg.cm[jseg].y;
                            track.mathematica_arrows[iter][tracker+5]=100*seg.cm[jseg].z;
                            tracker+=6;
                            if (ovlp.motor_direction[iseg][jseg]==1){
                                track.cross_type[iter][typetrack]=1;
                                typetrack++;
                            }
                            else if (ovlp.motor_direction[iseg][jseg]==-1){
                                track.cross_type[iter][typetrack]=-1;
                                typetrack++;
                            }
                            else{
                                printf("MOTOR not specified in mathematica arrows!\n");
                                exit(0);
                            }
                        }
                        
                        else {
                            printf("ERRROR in get_mathematica_arrow\n");
                            if (ovlp.type[iseg][jseg]==BUNDLING) {
                                printf("CRITICAL!\n");
                            }
                        }
                    }
                    
                }
            }
        }
    }
    
    return ( (int)(((float)tracker-1)/6) );//return the number of arrows
}



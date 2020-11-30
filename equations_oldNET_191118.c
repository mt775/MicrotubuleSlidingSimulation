/* #################################################################### */
/* equations  -- generates a linear set of equations for MTs velocities */
/*               and solves them numerically using Numerical Recipes 
                 routines.                                              */

#include <stdio.h>
#include "global_var.h"
#include "nrutil.h" 
#include "ran1.h"
#include <time.h>
#include <math.h>
//#include <vecLib/clapack.h>
#include <stdlib.h>

#define REF (bound.type[1]==REFLECT || bound.type[1]==RFIXREF) 
#define SI (iseg-2*iMT+2)   /* seg = {1,2} */
#define SJ (jseg-2*jMT+2)

/*Define LAPACK function*/
//int sgesv_(int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);

/* ################################################################## */
/* ################################################################## */
/* ##########################FUNCTIONS DEFINITIONS################### */
/* ################################################################## */


// This function scans the right boundary region and returns the microtubule indices
// in the array r_MTs and the total number of touching Microtubules IMPORTANT r_MT starts at zero

int MT_rightb(int *r_MTs, float rxmax){
    int count_r=0, iMT;// Number of MTs touching rbound
    for ( iMT=1; iMT<=MT.number; iMT++) {
        //If MT number iMT penetrates the force wall to touch_depth it counts as a touch
        if (spring_r<0.000001) {
            if (  MT.cm[iMT].x + MT.length[iMT]/2 + touch_depth >= rxmax){
                count_r++;
                r_MTs[count_r]=iMT;
            }
        }
        else if (spring_r>0.000001){
            //if ( fabs(spring_r* (rbound0 - rxmax)) >= 0.1){
                if (  MT.cm[iMT].x + MT.length[iMT]/2 + touch_depth >= rxmax){
                    count_r++;
                    r_MTs[count_r]=iMT;
                }
            //}
        }
        else{
            printf("equations.c >> FORCE boundary not properly specified!\n");
            exit(0);
        }
    }
    return(count_r);// return number of MT in bound region
}

int MT_leftb_FORCE(int *l_MTs, float lxmin){
    int count_l=0, iMT;// Number of MTs touching rbound
    for ( iMT=1; iMT<=MT.number; iMT++) {
        //If MT number iMT penetrates the force wall to touch_depth it counts as a touch
        if (spring_l<0.000001) {
            if (  MT.cm[iMT].x - MT.length[iMT]/2 - touch_depth <= lxmin){
                count_l++;
                l_MTs[count_l]=iMT;
            }
        }
        else if (spring_l>0.000001){
            //if ( fabs(spring_l* (lbound0 - lxmin)) >= 0.1){
                if (  MT.cm[iMT].x - MT.length[iMT]/2 - touch_depth <= lxmin){
                    count_l++;
                    l_MTs[count_l]=iMT;
                }
            //}
        }
        else{
            printf("equations.c >> FORCE boundary not properly specified!\n");
            exit(0);
        }
    }
    
    return(count_l);// return number of MT in bound region
}


/*This function calculates the number of MT connected to the left NET
 The model runs as follows: MT in LB region may become connected if they are
 new to the region only ONCE
 It also saves the index of the bound MT to lMTs
 */
////////////////////////////////////////////////////////////////////////////
void NET_connect_lb( int *l_MTs ){
    
    float rand;
    int count_l=0, iMT;
    
    /*Walk through all filaments*/
    /*#############################################################*/
    for (iMT=1 ; iMT<=MT.number; iMT++) {

        /*Look if MT touches bound*/
        if ((MT.cm[iMT].x- 0.5*MT.length[iMT])-touch_depth < lbound0) {
            rand=ran1(&idum);
            
            /*If MT is new to region (0) it connects to the NET*/
            
            if ( MT.is_in_net[iMT] == 0 ) {
                    MT.is_in_net[iMT]=2;
                    count_l++;
                    l_MTs[count_l]=iMT;
            }
            /*IF MT is connected (2) it may be released and can not attach again*/
            else if (MT.is_in_net[iMT] == 2 ){
                if (rand <= Lnet.dissociation){
                    MT.is_in_net[iMT]=1;

                }
                /*If filament stays connected count as touching bound*/
                else{
                    count_l++;
                    l_MTs[count_l]=iMT;
                }
            }
                
                
        }
        
        /*IF MT is away from bound region its status gets reset to 0
         so that it may connect to NET again upon entering bound region*/
        else MT.is_in_net[iMT]=0;
    }
    /*#############################################################*/
    if (!(iter%wfreq)){
        printf("There are %i MTs connected to the Cortex (NET)\n",count_l);
    }
}

void ALLNET_connect_lb( int *l_MTs ){
    
    float rand;
    int count_l=0, iMT;
    
    /*Walk through all filaments*/
    /*#############################################################*/
    for (iMT=1 ; iMT<=MT.number; iMT++) {
        
        rand=ran1(&idum);
        
        /*If MT is new to region (0) it may connect to the NET
         OR stay unconnected (1) until it leaves the bound region*/
        
        if ( MT.is_in_net[iMT] == 0 ) {
            /*If filament becomes connected count it as touchin the NET
             */
            if (rand <= Lnet.association ) {
                MT.is_in_net[iMT]=2;
                count_l++;
                l_MTs[count_l]=iMT;
            }
        }
        /*IF MT is connected (2) it may be released*/
        else if (MT.is_in_net[iMT] == 2 ){
            if (rand <= Lnet.dissociation)
                MT.is_in_net[iMT]=0;
            /*If filament stays connected count as touching bound*/
            else{
                count_l++;
                l_MTs[count_l]=iMT;
            }
        }
    }
    /*#############################################################*/
    if (!(iter%wfreq)){
        printf("There are %i MTs connected to the Cortex (NET)\n",count_l);
    }
}

////////////////////////////////////////////////////////////////////////////
/*The following function calculates the number of MTs subjected to CONFLX bound
 This means those MTs will become fixed on some velocity vel_flux if they are new to the system*/
////////////////////////////////////////////////////////////////////////////
int MT_leftb_CONFLX (int *l_MTs){
    int count_l=0,iMT;
    /*Walk through all filaments*/
    /*#############################################################*/
    for (iMT=1 ; iMT<=MT.number; iMT++) {
        
        /*IF MTs originate from last MT pull make them subject to bound velocity*/
        if ( (iter-MT.iter0[iMT]) < newfreq ) {
            count_l++;
            l_MTs[count_l]=iMT;
        }
    
    }
    /*#############################################################*/

    return (count_l);
}
////////////////////////////////////////////////////////////////////////////
                    
// The same for left bound without conditions!
////////////////////////////////////////////////////////////////////////////
void MT_leftb(int *l_MTs){

    int count_l=0, iMT;// Number of MTs touching lbound
    
    for (iMT=1; iMT<=MT.number; iMT++) {
        //If MT number iMT penetrates wall  it counts as a touch
        if (  MT.cm[iMT].x - MT.length[iMT]/2 <= lbound0){
            count_l++;
            l_MTs[count_l]=iMT;
        }
    }
}
////////////////////////////////////////////////////////////////////////////

void checkMT(int *l_MTs, int *r_MTs, int *count_l, int *count_r){
    int iMT=1;
    int jMT=1;
    int foo=0;

    while (l_MTs[iMT]) {
        while (r_MTs[jMT]) {
            if (r_MTs[jMT] == l_MTs[iMT]) {
                printf("Filament touching both boundaries detected %i\n",l_MTs[iMT]);
                *count_r=*count_r-1;
                *count_l=*count_l-1;
                for (foo=jMT; foo<=*count_r; foo++) {
                    r_MTs[foo]=r_MTs[foo+1];
                }
                r_MTs[*count_r+1]=0;
                for (foo=iMT; foo<= *count_l; foo++) {
                    l_MTs[foo]=l_MTs[foo+1];
                }
                l_MTs[*count_l+1]=0;
                
            }
            jMT++;
        }
        iMT++;
    }
    
}
/*#############################################################*/

/////////THIS FUNCTIONS CALCULATES THE MOTOR INTERACTION MATRIX//////
/*Written for the general, but not implemented case of varying lambda*/
void get_motor_forces (float **A, float *B){

    float boxlength= (MT.xmax_real+0.5*MT.length[MT.imax_real]-MT.xmin_real+0.5*MT.length[MT.imin_real])/(float)nbox;
    struct BOXsegmentdata {
        float length_i,length_j, leftend_i, leftend_j;
    } boxseg;
    float imfii,imfij,imfbi; /* image factors */
    int n,ihigh,iMT,jMT,iseg,jseg,ibox;
    n=MT.number;
    float f0,vel0;
    float sz; //overlapsize,start of overlap
    float left_domain_end=MT.xmin_real-0.5*MT.length[MT.imin_real]; //the left end of the domain
    float f0mult=1.;

    /* default values, image factors  for no reflecting bound*/
    imfii=1; imfij=1;imfbi=1;
    
    /*Scan through all boxes behold that the length of the domain is NOT constant
     in the case of Force bound or NON....*/
    /*#################################################################################*/
    for (ibox=1; ibox<=nbox; ibox++) {
        /*#################################################################################*/
        for (iseg=1;iseg<=2*n;iseg++) {
            
            // only carry on in this loop if the segment is existing and overlapping with ibox
            if (seg.length[iseg]<0.0001) continue;
            boxseg.length_i = Interval_overlap( (ibox-1)*boxlength+left_domain_end ,seg.cm[iseg].x-seg.length[iseg]/2. , boxlength , seg.length[iseg] );
            if (boxseg.length_i<0.0001) continue;
            boxseg.leftend_i = max( (ibox-1)*boxlength+left_domain_end ,  seg.cm[iseg].x-seg.length[iseg]/2. );
            
            
            iMT=segindex(iseg);// seg index returns returns iMT number of the segment
            
            /*If we also simulate actin we need to add one term to A[iMT][iMT]
             and one to B[iMT]*/
            if (ACTIN==TRUE) {
                if (ovlp.actin[iMT]==1) {

                    A[iMT][iMT]+=boxseg.length_i*dist_dyn[ibox]*fstall_d/vel0_d;
                    B[iMT]-=boxseg.length_i*fstall_d*dist_dyn[ibox]*MT.direct[iMT].x*ovlp.actin[iMT];
                }
                if (ovlp.actin[iMT]==-1) {

                    A[iMT][iMT]+=boxseg.length_i*dist_kin[ibox]*fstall_k/vel0_k;
                    B[iMT]-=boxseg.length_i*fstall_k*dist_kin[ibox]*MT.direct[iMT].x*ovlp.actin[iMT];
                }
            }
            
            /*#################################################################################*/
            for (jseg=1;jseg<=2*n;jseg++){
                
                // only carry on in this loop if the segment is existing and overlapping with ibox
                if (seg.length[jseg]<0.0001) continue;
                boxseg.length_j = Interval_overlap( (ibox-1)*boxlength+left_domain_end ,seg.cm[jseg].x-seg.length[jseg]/2., boxlength , seg.length[jseg] );
                if (boxseg.length_j<0.0001) continue;
                boxseg.leftend_j = max( (ibox-1)*boxlength+left_domain_end ,  seg.cm[jseg].x-seg.length[jseg]/2. );


                jMT=segindex(jseg);
                
                // Calculate the overlap of the segments in ibox
                if (iseg==jseg) sz=0;
                else sz=Interval_overlap(boxseg.leftend_i,  boxseg.leftend_j, boxseg.length_i, boxseg.length_j);
                //
                
                if (sz>EPSI && ovlp.type[iseg][jseg]!=ZERO) {
                    //Ucomment next line for non ovlp length dependence
                    //sz=1.;
                    
                    //Determine OVLP motor direction type and set parameters respectively
                    /*#################################################################################*/
                    if (ovlp.motor_direction[iseg][jseg] == 1) {
                        f0=dist_dyn[ibox]*fstall_d;
                        vel0=vel0_d;
                    }
                    else if (ovlp.motor_direction[iseg][jseg] == -1) {
                        f0=dist_kin[ibox]*fstall_k;
                        vel0=vel0_k;
                    }
                    else{
                        printf("direct: %i iMT:%i jMT:%i type:%i\n", ovlp.motor_direction[iseg][jseg], iMT, jMT, ovlp.type[iseg][jseg]);
                        printf("\nERROR in motor direction exiting....\n");
                        exit(0);
                    }
                    /*#################################################################################*/
                    
                    //SCAN FOR MOTOR-TYPES
                    /*#################################################################################*/
                    if (ovlp.type[iseg][jseg]==BIPOLAR) /*Bipolar Motors */{
                        
                        if (0){/*Lingering bipolar motor*/
                            if (seg.direct[iseg].x*seg.direct[jseg].x>0) {
                                
                                A[iMT][iMT]+=sz*(xim0+f0/vel0);
                                A[iMT][jMT]-=sz*f0/vel0;
                                /*plusend directed bipol motors will bundle +ends*/
                                
                                    if ( seg.direct[iseg].x*ovlp.motor_direction[iseg][jseg]>0 && fabs(MT.cm[iMT].x+MT.length[iMT]/2 - MT.cm[jMT].x - MT.length[jMT]/2)>0.1 ){
                                        if ((MT.cm[iMT].x+MT.length[iMT]/2) > (MT.cm[jMT].x+MT.length[jMT]/2))
                                            B[iMT]-=sz*f0;
                                        else
                                            B[iMT]+=sz*f0;
                                    }
                                    else if (seg.direct[iseg].x*ovlp.motor_direction[iseg][jseg]<0 && fabs(MT.cm[iMT].x-MT.length[iMT]/2 - MT.cm[jMT].x + MT.length[jMT]/2)>0.1){
                                        if ((MT.cm[iMT].x-MT.length[iMT]/2) > (MT.cm[jMT].x-MT.length[jMT]/2))
                                            B[iMT]-=sz*f0;
                                        else
                                            B[iMT]+=sz*f0;
                                    }
                                
                            }
                            
                            else {
                        /* Anti-parallel pair */
                                A[iMT][iMT]+=sz*(xim0+0.5*f0/vel0_k);
                                A[iMT][jMT]-=0.5*sz*f0/vel0_k;
                                B[iMT]-=sz*f0*seg.direct[iseg].x*ovlp.motor_direction[iseg][jseg];
                            }
                        }
                        else{/*non lingering bipol*/
                            if (seg.direct[iseg].x*seg.direct[jseg].x>0) {
                                A[iMT][iMT]+=sz*(xim0+f0mult*0.5*f0/vel0_k);
                                A[iMT][jMT]-=0.5*f0mult*sz*f0/vel0_k;
                            }
                            else {
                                A[iMT][iMT]+=sz*(xim0+0.5*f0mult*f0/vel0_k);
                                A[iMT][jMT]-=0.5*f0mult*sz*f0/vel0_k;
                                B[iMT]-=sz*f0mult*f0*seg.direct[iseg].x*ovlp.motor_direction[iseg][jseg];
                            }
                            
                        }
                    }
                    /* BUNDLING passive BIPOL motor */
                    else if (ovlp.type[iseg][jseg]==BUNDLING) {
                        /*Bundling only*/
                            A[iMT][iMT]+=sz*100*(xim0+0.5*f0/vel0_k);
                            A[iMT][jMT]-=0.5*sz*100*f0/vel0_k;
                    }
                    else if (ovlp.type[iseg][jseg]==TWOMOTORS){
                        
                        A[iMT][iMT]+=sz*(xim0+f0/vel0);
                        A[iMT][jMT]-=sz*f0/vel0;
                        B[iMT]-=sz*f0*seg.direct[iseg].x;
                
                        A[iMT][iMT]+=sz*(xim0+f0/vel0);
                        A[iMT][jMT]-=sz*f0/vel0;
                        B[iMT]+=sz*f0*seg.direct[jseg].x;
                    }
                    else /* Unipolar Motors */{
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
                    
                        if ((iMT==ihigh && ovlp.type[iseg][jseg]==LEGUP)
                            || (jMT==ihigh && ovlp.type[iseg][jseg]==LEGDOWN)) {
                    /* Motor legs are on iseg/iMT, and their heads on jseg/jMT
                     */
                        /*----------------------------------*/
                        /* In case of a reflection boundary */
                        /* Determine sign of image forces   */
                            if (REF){
                                
                                if(((SI==1 && SJ==2)) && iMT==jMT){
                                    imfii=4.0;
                                    imfij=0.0;
                                    imfbi=2.0;
                                }
                                else if (((SI==1 && SJ==1) || (SI==2 && SJ==2)
                                    || (SI==2 && SJ==1)) &&
                                     iMT==jMT) {
                                    imfii=0.0;
                                    imfij=0.0;
                                    imfbi=0.0;
                                }
                                else if (SI==2 && SJ==1 && iMT!=jMT)
                            /* im,b: "leg on image" j==head */{
                                    imfii=1.0;
                                    imfij=-1.0;
                                    imfbi=-1.0;
                                }
                                else if (SI==1 && SJ==2 && iMT!=jMT)
                            /* im,a: "head on image" j==head */{
                                    imfii=1.0;
                                    imfij=-1.0;
                                    imfbi=1.0;
                                }
                                else if (SI==2 && SJ==2 && iMT!=jMT){
                            /* im,c: "both on image" j==head */
                                    imfii=1.0;
                                    imfij=1.0;
                                    imfbi=-1.0;
                                }
                                else if (SI==1 && SJ==1 && iMT!=jMT){
                                    
                                    imfii=1.0;
                                    imfij=1.0;
                                    imfbi=1.0;
                                }
                            }
                        /* ------- End of Image conditions ---- */
			   	 
                            A[iMT][iMT]+=sz*(xim0+f0/vel0)*imfii;
                            A[iMT][jMT]-=sz*f0/vel0*imfij;
                            B[iMT]-=sz*f0*seg.direct[iseg].x*imfbi*ovlp.motor_direction[iseg][jseg];

                            //printf("A[%i][%i]=%f B[%i]=%f iseg=%i jseg=%i\n",iMT,iMT,A[iMT][iMT],iMT,B[iMT],iseg,jseg);
                            
                        }
                        else /* Motor heads are on iseg/iMT, and their legs on
                          jseg/jMT */{
                        /*----------------------------------*/
                        /* In case of a reflection boundary */
                        /* Determine sign of image forces   */ 
                            if (REF) {
                                if (((SI==1 && SJ==2)) && iMT==jMT){
                                    imfii=4.0;
                                    imfij=0.0;
                                    imfbi=2.0;
                                }
                                else if (((SI==1 && SJ==1) || (SI==2 && SJ==2)
                                      || (SI==2 && SJ==1))
                                     && iMT==jMT){
                                    imfii=0.0;
                                    imfij=0.0;
                                    imfbi=0.0;
                                }
                                else if (SI==2 && SJ==1 && iMT!=jMT)
                            /* im,a: "head on image" j==leg */{
                                    imfii=1.0;
                                    imfij=-1.0;
                                    imfbi=-1.0;
                                }
                                else if (SI==1 && SJ==2 && iMT!=jMT)
                            /* im,b: "leg on image" j==leg */{
                                    imfii=1.0;
                                    imfij=-1.0;
                                    imfbi=1.0;
                                }
                                else if (SI==2 && SJ==2 && iMT!=jMT)
                            /* im,c: "both on image" j==leg */{
                                    imfii=1.0;
                                    imfij=1.0;
                                    imfbi=-1.0;
                                }
                                else if (SI==1 && SJ==1 && iMT!=jMT){
                                /* both segments are real j==leg*/
                                    imfii=1.0;
                                    imfij=1.0;
                                    imfbi=1.0;
                                }
                            }
                        /* ------- End of Image conditions ---- */
                        //legs on

                            A[iMT][iMT]+=sz*(xim0+f0/vel0)*imfii;
                            A[iMT][jMT]-=sz*f0/vel0*imfij;
                            B[iMT]+=sz*f0*seg.direct[jseg].x*imfbi*ovlp.motor_direction[iseg][jseg];
                              
                                  //printf("A[%i][%i]=%f B[%i]=%f iseg=%i jseg=%i\n",iMT,iMT,A[iMT][iMT],iMT,B[iMT],iseg,jseg);
                              
                        
                        }
                        /*#################################################################################*/
                        //end motor type scan
                    }
                }
            }
            /*#################################################################################*/
        }
        /*#################################################################################*/
    }
    /*#################################################################################*/
}
/////////THIS FUNCTION INITIALIZES THE EXT FORCE////////////////////
float set_external_force_r(float rxmax)
{
    /*IF k>0 spring is applied*/
    if (spring_r>0.000001) {
        if (iter > force_it0) {
        
            //printf("force in iter%i =%f\n",iter,-spring_r);
            return(-spring_r*(rxmax-rbound0));
        }
        else return(0.);
    }
    else {
        if (iter > force_it0) {/* To avoid the case where all MTs leave the system */
            
            return(ext_fr);
        }
        else return(0.);
    }
}

float set_external_force_l(float lxmin)
{
    if (spring_l>0.000001 ) {
        if (iter > force_it0  && fabs(spring_l*(lxmin-lbound0))>1) {
            
            return(-spring_l*(lxmin-lbound0));
        }
        else return(0.);
    }
    else {
        if (iter > force_it0 ) {/* To avoid the case where all MTs leave the system */
            //printf("extf=%f\n",ext_fl);
            return(ext_fl);
        }
        else return(0.);
    }
}

///////THIS FUNCTION ADDS THE BOUNDARY FORCES TO THE MATRIX/////////
void add_boundary_force(float **A, float *B, int *r_MTs, int *l_MTs, int count_r, int count_l, float f_r, float f_l) {
    
    int iMT=1;
    
    // first add the boundary force from the right to each force velocity relation
    // Like xinsol*v(1)=(...)+f1 !!
    
    while (r_MTs[iMT])
    {
        
        A  [ r_MTs[iMT]  ][MT.number + iMT] = -1;
        iMT++;
    }
   
    
    //GENERATE EQUATIONS
    iMT=1;
    // ALL FORCES INDUCED BY BOUND MUST SUM UP TO EXTERNAL FORCE f1+f2+..=Fext

    while (r_MTs[iMT])
    {
        A [MT.number+1][ MT.number + iMT ] = 1;
        iMT++;
    }
        //If at least one is touching invoke ext force
    if (r_MTs[1]) B[MT.number+1] = f_r;
    
    // if more than one MT touches boundary their velociries have to syncronize
    // V1-Vn=0 for all n iMT=2 because these equations only make sense if more than one touches bound
    if (r_MTs[1]) {
        iMT=2;
        
        while(r_MTs[iMT]) {
            
            A [MT.number+iMT][ r_MTs[1] ]   = - 1 ;
            
            A [MT.number+iMT][ r_MTs[iMT] ] =  1 ;
            iMT++;
        }
    }

    
    /*Max implemented 2 different possibilities of making the filaments get
     caught in the NET. One is below and raises the matrix equations by the number of MT 
     in left bound. The other just raises the drag coeff of the MTs in NET so that they cannot move 
     (a lot better for performance)*/
    if (0) {
        iMT=1;
        while (l_MTs[iMT]) {
            A [MT.number+count_r+iMT][l_MTs[iMT] ] = 1;
            B[MT.number+count_r+iMT] = 0;
            iMT++;
        }
        
    }
    
    /*For CONFLX bound force the boundary MTs to vel_flux*/
    if (bound.type[0]==CONFLX) {
        iMT=1;
        while (l_MTs[iMT]) {
            A [MT.number+count_r+iMT][l_MTs[iMT] ] = 1;
            B[MT.number+count_r+iMT] = vel_flux;
            iMT++;
        }
        
    }
    //Same as above for left force bound
    else if (bound.type[0]==FORCE) {
        iMT=1;
        while (l_MTs[iMT]){
            A  [ l_MTs[iMT] ][MT.number+ count_r + iMT] = -1;
            A [MT.number+count_r+1][ MT.number +count_r + iMT ] = 1;
            iMT++;
        }
        if (l_MTs[1]) B[MT.number+count_r+1] = f_l;
        
        iMT=2;
        while(l_MTs[iMT]) {
            A [MT.number+count_r+iMT][ l_MTs[1] ]   = - 1 ;
            
            A [MT.number+count_r+iMT][ l_MTs[iMT] ] =  1 ;
            iMT++;
        }
    }
}

//THIS FUNCTION CAN SOLVE A LIN EQUATION
void solve_lin_eq ( float **A, float *B, int n) {
    
    // VARIABLES NEEDED
    int iMT,jMT;
    float d;
    float **svA;
    float *svB;
    int *index;
    svA=matrix(1,n,1,n);
    svB=vector(1,n);
    index=ivector(1,n);
    
    for (iMT=1;iMT<=n;iMT++)
    {
        index[iMT]=0.0;
        svB[iMT]=0.0;
        for (jMT=1;jMT<=n;jMT++) svA[iMT][jMT]=0.0;
    }
    
    /* Solve linear set of equations using LU decomposition routine */
    /* Numerical Recipes functions: ludcmp to get the LU decomposition
     and then lubksb to get the answer.                           */
    /* We first copy A and B to A1 and B1 since these are changed
     by the above routines. */
    
    for (iMT=1;iMT<=n;iMT++)
    {
        svB[iMT]=B[iMT];
        for (jMT=1;jMT<=n;jMT++)
            svA[iMT][jMT]=A[iMT][jMT];
    }
    
    /* solve the equations */
    ludcmp(A,n,index,&d);        /* Decompose matrix A1[n][n] */
    /* index, d, is output */
    lubksb(A,n,index,B);        /* A is input, as the result of ludcmp */
    /* B is input, and returns with the solution*/
    
    mprove(svA,A,n,index,svB,B); /* Improve the solution which returns with B */
    
    free_vector(svB,1,n);
    free_matrix(svA,1,n,1,n);
    free_ivector(index,1,n);

}

int solve_lin_eq_LAPACK ( float **A, float *B, int n ){
    
    float Afortran[n*n], Btmp[n];
    int i,j, nb=1, check, Bpivot[n];
    int fail=0;
    
    for (i=0; i<n; i++){		/* to call a Fortran routine from C we */
                                /* have to transform the matrix to an array*/
        Btmp[i]=B[i+1];
        for(j=0; j<n; j++){
            Afortran[j+n*i]=A[j+1][i+1];
        }
    }

    /*here for general matrix, but if A is symmetric (no-force case) one could
     actually use a special subroutine for those like ssysv_(...) (LAPACK)*/
    //sgesv_(&n, &nb, Afortran, &n, Bpivot, Btmp, &n, &check);
    
    if (check!=0) {
        //printf("sgesv failed!\n");
        //printf("check value=%i\n",check);
        fail=1;
    }
    
    if (fail==0){
        for (i=0; i<n; i++)
            B[i+1]=Btmp[i];
    }
    return(fail);
    
}

// THIS FUNCTION CHECKS IF EXT FORCES ARE CALCULATED CORRECTLY
int check_external_forces( float *B, int count_l, int count_r, float f_r, float f_l ){
    int iMT;
    
    
    float opp_force = 0;         // holds the size of force that is calculated in incorrect form.
    int   max_opp_force_i  = 0;  // holds the index of force that was calculated in incorrect form.
    // "max_opp_force_i" is  returned from the function;
    // if all forces are correct, 0 is returned, and no forther
    // correction will be done;
    
    float max_opp_force    = 0;  // hold the max incorrect force, to be fixed first.
    
    
    /*Check at right boundary*/
    for (iMT = MT.number+1; iMT <=MT.number+count_r; iMT++){
        //if ext force is calculated to be in wrong direction
        if(  ( B[iMT]>0.1 && f_r<-0.1) || ( B[iMT]<-0.1 && f_r>0.1) ){

            opp_force =   B[iMT];
            if (opp_force < 0 )  opp_force = opp_force *(-1);
        }
        if (opp_force > max_opp_force) {
            max_opp_force = opp_force ;
            max_opp_force_i = iMT-MT.number;
        }
    }
    /*Check at left boundary*/
    for (iMT = MT.number+count_r+1; iMT <=MT.number+count_r+count_l; iMT++){
        //if ext force is calculated to be in wrong direction
        if(  ( B[iMT]>0.1 && f_l<-0.1 ) || ( B[iMT]<-0.1 && f_l>0.1 ) ){
            
            opp_force =   B[iMT];
            if (opp_force < 0 )  opp_force = opp_force *(-1);
        }
        if (opp_force > max_opp_force) {
            max_opp_force = opp_force ;
            max_opp_force_i = iMT-MT.number;
        }
    }
    
    if ( max_opp_force_i  && (1 )){ // FOR DEBUGGING, change to (1)
        printf ("correcting opposite external forces in iter:%i\n",iter);
        printf("for iMT%i cm=%f f_l=%f f_r=%f maxoppforce=%f \n",max_opp_force_i, MT.cm[max_opp_force_i].x, f_l,f_r,max_opp_force);
    }
    return max_opp_force_i;
}

//DEBUG FUCNTIONS
void print_external_forces(int *r_MTs, float *B)
{
    int iMT = 1;
    while (r_MTs[iMT])
	{
        printf("r_ex_f on MT %i= %f\n",r_MTs[iMT], B[ MT.number  + iMT]);
        iMT++;
	}
}

void print_linear_equation(int n, float **A, float *B)
{
    int i,j;
    printf ("\n              Linear Equation  A*V=B\n");
    printf("MT.number=%d iter[%d]\n",MT.number,iter);
    printf("n=%d\n",n);
    printf("           ---------------------------------\n");
    for (i=1;i<=n;i++)
    {
        
        printf("[");
        if ( (A[i][1]> -10)  && (A[i][1]<10) )
            printf("%f\t",A[i][1]);
        else
            printf("%f\t",A[i][1]);
        printf("\t\b");
        
        for (j=2;j<=n;j++)
        {
            if ( (A[i][j])>-10  && (A[i][j])<10 )
                printf("%f\t",A[i][j]);
            else
                printf("%f\t",A[i][j]);
            printf("\t\b\b");
        }
        
        //printf ("\t\b\b\b]");
        if (i <=MT.number) printf ("[V%d] = ",i);

        if ( i > MT.number)
            printf ("[Rf%d]= ",i - MT.number);
        
        
        //printf ("B[%d] ",i);
        printf ("\t\b\b\b\b\b[%f]\n",B[i]);
    } 
}

void print_MT_velocities()
{
    int iMT;
    printf("MT's Velocities\n");
    for (iMT=1;iMT<=MT.number;iMT++)
    {
        printf("MT[%i]\t%.3f.\t ",iMT,MT.vel[iMT].x);
        if (!(iMT%5))
            printf("\n");
    }
    printf("\n");
	
}


/* ################################################################## */
/* ################################################################## */
/* ################################################################## */
////FUNCTIONS END

/* ################################################################## */
/* The equations have the form A x = B, where A is a matrix of 
   cooeficiants, which will be derived here, x=vel/vel0 is the 
   vector of velocities of all MTs, and B is another vector of 
   constants that will be determined here.                            */



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
///////////////////////////MAIN EQUATION//////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void equations() 
{
    int n, nb=1,ihigh,iMT,jMT,iseg,jseg, fail;
    float **A, *B, *norm;
    float sz; /* overlapping size; Lij - in my notes */
    float f0, f_l,f_r; /* Total stall force per unit length */
    struct point fr;
    float d;
    float avforce;
    float lxmin, rxmax;
    clock_t t;

    int *l_MTs, *r_MTs; //arrays giving the index number of the ith MT which touches boundary
    int count_r=0, count_l=0; //simple counters to count number of MT in bound region
    int *faulty_MT; //array which saves the indexing of the microtubules pushing against the WALL
    int wrong_ext_F, counter, move_into_wall;
    int inf_loop_abort; // temporary WALL solution
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
    //just in case everybody is touching/faulty full array(the easiest not the fastest)
    // "+1" is needed so that the while >0 loop aborts if it is outside the MT number domain
    //Otherwise it will not as for iMT> MT.number r_MTs is not defined There could be more
    //elegant solutions but I will stick with this
    l_MTs=ivector(1,MT.number+1);
    r_MTs=ivector(1,MT.number+1);
    faulty_MT=ivector(1,MT.number+1);
    emptyint(l_MTs,1,MT.number+1);
    emptyint(r_MTs,1,MT.number+1);


    emptyint(faulty_MT,1,MT.number+1);
    
    //Get most extreme coordinate values values
    lxmin=MT.xmin_real-0.5*MT.length[MT.imin_real];
    rxmax=MT.xmax_real+0.5*MT.length[MT.imax_real];
//////////////////////////////////////////////////////////////////////////////
//////////////////////////Get numbers of MT in boundary if FORCE or WALL /////
//////////////////////////////////////////////////////////////////////////////
    /*RIGHT BOUND*/
    if (bound.type[1]==FORCE)                            count_r = MT_rightb( r_MTs , rxmax);
        // get number and index of MT in bound region
    
    /*LEFT BOUND*/
    if (bound.type[0]== WALL || bound.type[0]==STICKY )  MT_leftb( l_MTs );
    
    else if (bound.type[0]==NET)                         NET_connect_lb( l_MTs );
    
    else if (bound.type[0]==CONFLX)                      count_l=MT_leftb_CONFLX(l_MTs);
    
    else if (bound.type[0]==FORCE)                       count_l=MT_leftb_FORCE(l_MTs, lxmin);
    
    //check if one MT touches BOTH boundaries and delete then
    
    //checkMT(l_MTs,r_MTs, &count_l, &count_r);
    //printf("cl=%i\n",count_l);
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
    inf_loop_abort=0;
    
        //THE MAIN LOOP
    do{
       
        
        inf_loop_abort++;
        //Loop switch of by default. will be set to one if boundary WALL is on
        move_into_wall=0;
        /*################################################################*/
        
        /*################################################################*/
        /*Note count_l != 0 only for CONFLX*/
        n=MT.number+count_r+count_l; //every boundary condition raises the equation number
        A=matrix(1,n,1,n);    /* allocate memory for matrices */
        norm=vector(1,n);
        B=vector(1,n);
        

        for (iMT=1;iMT<=n;iMT++)
        {
            B[iMT]=0.0;
            for (jMT=1;jMT<=n;jMT++) A[iMT][jMT]=0.0;
        }

        /*################################################################*/
        /* Generate equations */
        /* The loop runs over all segments, some MTs have two segments
         because of application of periodic or reflecting boundary conditions.*/
        /* Interactions of different segments of one MT with those of another
         MT are additive. */
        // The following loop adds the motor forces to the A Matrix each with prefac f0/v0
       
        /*################################################################*/
        get_motor_forces(A,B);
        /*################################################################*/

        
        //THE NEXT FUNCTION FILLS THE MATRIX WITH THE FORCE BOUNDARY CONDITIONS
        /*################################################################*/
        f_r=set_external_force_r(rxmax);
        f_l=set_external_force_l(lxmin);
        //printf("\n n=%i f_r=%f f_l=%f\n",n,f_r,f_l);
        add_boundary_force( A, B, r_MTs, l_MTs, count_r, count_l, f_r , f_l );
        
        /*################################################################*/

        /*Add drag coeff to the matrix diagonal*/
        /*and add random temperature force*/
        for (iMT=1;iMT<=MT.number;iMT++) {
            A[iMT][iMT]+=xisol*MT.length[iMT];
            B[iMT]+=(2*ran1(&idum)-1)*kBT;
        }



        /*################################################################*/
        /*Raise drag coeff for NET and STICKY (there exists another approach for NET,
         look up in add boundary force, but this one here is quicker) */
        if ( bound.type[0]==STICKY) {
            iMT=1;
            while (l_MTs[iMT]) {
                A[ l_MTs[iMT] ][ l_MTs[iMT] ]-= xisol*MT.length[ l_MTs[iMT] ];
                A[ l_MTs[iMT] ][ l_MTs[iMT] ]+= sticky_drag*(lbound0 - MT.cm[l_MTs[iMT]].x + MT.length[ l_MTs[iMT] ]/2);
                iMT++;
            }
        }
        else if ( bound.type[0]==NET ) {
            iMT=1;
            while (l_MTs[iMT]) {
                A[ l_MTs[iMT] ][ l_MTs[iMT] ]-= xisol*MT.length[ l_MTs[iMT] ];
                A[ l_MTs[iMT] ][ l_MTs[iMT] ]+= sticky_drag*MT.length[ l_MTs[iMT] ];
                iMT++;
            }
        }
        /*################################################################*/
        // The next step looks for filaments which move towards WALL and raises their drag coeff
        //iMT=1;
        //while (faulty_MT[iMT]) {
           
            //If a filament is touching the wall and its velocity is directed in negative direction raise drag coeff so that filament does not move
          //  A[ faulty_MT[iMT] ][ faulty_MT[iMT] ] += 99999*xisol*(lbound0 - MT.cm[faulty_MT[iMT]].x + MT.length[ faulty_MT[iMT] ]/2);
            //Set to zero so as its not faulty anymore
            //faulty_MT[iMT]=0;
            //iMT++;
        
        //}
        /*################################################################*/
        
        /* So far A has units of force/velocity, and B has units
         of force. We now normalize the parameters Aij by
         Aii, and B by Aii */
        
        
        for (iMT=1;iMT<=n-count_l-count_r;iMT++)
        {
            norm[iMT]=A[iMT][iMT];
            B[iMT]/=(norm[iMT]);
            for (jMT=1;jMT<=n;jMT++)
            {
                A[iMT][jMT]/=norm[iMT];
            }
        }

        
        
                //THE NEXT FUNCTION SOLVES the LINEAR EQUATIONS
        //The result will be stored in the B vector so that the first MT.number rows are the
        //resulting velocities and the last r_count number of rows the forces induced by bound
        /*################################################################*/

        
        
        fail=1;
        //print_linear_equation(n,A,B);
        
        //printf("\n %i",iter);
        //printf("\n external forces\n");
        //print_external_forces(r_MTs, B);
        //fail=solve_lin_eq_LAPACK(A,B,n); /*CONSIDERABLY QUICKER (TESTED) but it seems to break down  for large MT numbers*/
        /*If LAPACK failed try nrutil*/
        if (fail==1){
            //printf("LAPACK failed... Trying NRutil...");
            //print_linear_equation(n,A,B);
            solve_lin_eq( A, B, n );
            
            //exit(0);
        }
            /*###################################*/
        //print_linear_equation(n,A,B);
        
        

        //CHECK IF FORCES WERE CALCUTATED CORRECTLY
        //In some cases the calculated forces do not have the same sign as the boundary force
        //which is unphysical (returns 0 if everything ok) these MT get deleted
        //printf("\n velocities after solving\n");
        //print_external_forces(r_MTs, B);

        wrong_ext_F = check_external_forces(B,count_l,count_r, f_r, f_l);
       
        if (wrong_ext_F>0) {
            if (wrong_ext_F <= count_r) {
                r_MTs[ wrong_ext_F ] = r_MTs[count_r];
                r_MTs[count_r] = 0;
                count_r --;
                
            }
            else if (wrong_ext_F > count_r) {
                l_MTs[ wrong_ext_F-count_r ] = l_MTs[count_l];
                l_MTs[count_l] = 0;
                count_l --;
            }
            else{
                printf("equations.c >> ERROR\n");
                exit(0);
            }
            
        }
        
        
        
        /* Now copy the solution to MT.vel */
        
        for (iMT=1;iMT<=MT.number;iMT++){
            MT.vel[iMT].x=B[iMT];
            if (fabs(MT.vel[iMT].x)>1) {
                printf("Too large velocity! Setting to vmax=10 maybe decrease timestep? %f\n",MT.vel[iMT].x);
                if (   MT.vel[iMT].x>0)
                    MT.vel[iMT].x=1;
                else
                    MT.vel[iMT].x=-1;

            }
            // EDIT 120718 MAX: Sometimes the equation solving gives ridiculous results. In that case set all nans to 0
            if isnan(MT.vel[iMT].x) {
                printf("error in solving equation:: setting all nans to 0!\n");
                MT.vel[iMT].x=0.;
            }
        }
        
        
        
    
        /*################################################################*/
        //Free memory for next round or final
        free_vector(B,1,n);
        free_vector(norm,1,n);
        free_matrix(A,1,n,1,n);
        

        //If caught in infinite loop because raising drag coeff leads to another MT inverting velocity
        // So far no other solution on my mind....
        if (inf_loop_abort >= 5 && 1 ) {
            if ( bound.type[0]==WALL){
                move_into_wall=0;
                printf("it%i Inf loop abort necessary...\n",iter);
            }
            else if (bound.type[0]==FORCE){
                printf("it%i takes long to correct...\n",iter);
                
            }
            else{
                printf("it%i takes long to correct...WTF no force\n",iter);
            }
        }
    
    }while (wrong_ext_F || move_into_wall); //redo until no force is calculated wrongly and all WALL touching MTs are slowed doen

    //print_MT_velocities();

    free_ivector(l_MTs ,1,MT.number+1);
    free_ivector(r_MTs ,1,MT.number+1);
    free_ivector(faulty_MT ,1,MT.number+1);
    
}

    

  

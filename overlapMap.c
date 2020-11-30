/* #################################################################### */
/* OverlapMap generates a map of overlapping MTs                         */
/* The map is stored in the struct ovlp                                 */
/* Each OVLP region is assigend an ovlp.type[iseg][jseg]: 
   {ZERO,LEGUP,LEGDOWN,BIPOLAR}                                         */
/* LEGUP -- if motors have their legs pointing upwards in the z-direction */    /* LEGDOWN -- vise versa */
/* BIPOLAR -- motors are bipolar (similarly directed motor heads)   */
/* ZERO -- if there is no overlap between iseg and jseg */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ran1.h"
#include "global_var.h"
#define SI1 (iseg1-2*iMT+2)   /* seg = {1,2} */
#define SJ1 (jseg1-2*jMT+2)
#define SI2 (iseg2-2*iMT+2)   /* seg = {1,2} */
#define SJ2 (jseg2-2*jMT+2)
#define REF (bound.type[1]==RFIXREF || bound.type[1]==REFLECT)

/*This functions calculates prbabilities due to a probability distribution*/
void get_Prob_active( int iseg, int jseg, float *prob, float dx){
    float iseg_l,jseg_l;
    iseg_l=seg.cm[iseg].x-seg.length[iseg]/2;
    jseg_l=seg.cm[jseg].x-seg.length[jseg]/2;
    float min_ovlp=max(iseg_l, jseg_l);
    float ov=Interval_overlap(iseg_l, jseg_l, seg.length[iseg], seg.length[jseg]);
    float length;
    
    int ibox=(int)floor((min_ovlp-MT.xmin_real+MT.length[MT.imin_real]/2)/dx)+1;
    /*a floating point error may lead to ibox=0,nbox has to be corrected*/
    if (ibox<=0 || ibox>nbox) {
        if (ibox==0) {
            ibox=1;
        }
        else if (ibox==nbox){
            ibox=nbox-1;
        }
        else{
            printf("FATAL ERROR in overlapmap!!\n");
            printf("minovlp=%f\n", min_ovlp);
            printf("dx=%f\n", dx);
            printf("box=%i\n", ibox);
            printf("MT.xmin=%f\n", MT.xmin_real-MT.length[MT.imax_real]/2);
            printf("MT.xmax=%f\n\n", MT.xmax_real+MT.length[MT.imin_real]/2);
            printf("MT.length=%f\n", MT.length[MT.imin_real]);
            printf("isegl=%f\n", iseg_l);
            printf("isegl=%f\n", jseg_l);
            printf("iseg=%i\n", iseg);
            printf("MTcm=%f\n",MT.cm[segindex(iseg)].x);
            printf("MT left end=%f\n",MT.cm[segindex(iseg)].x-MT.length[segindex(iseg)]/2);
            printf("checksum=%f\n",(min_ovlp-MT.xmin_real+MT.length[MT.imin_real]/2)/dx);
            exit(0);
        }
    }
  
    
    int i;
    for (i=0; i<=3; i++)
        prob[i]=0.;
    
    /*get ovlp with boxes*/
    length=Interval_overlap(min_ovlp, ((float)ibox-1)*dx+MT.xmin_real-MT.length[MT.imin_real]/2, ov, dx );

    do{
        prob[0]+=length*ProbAct_dist[ibox];
        prob[1]+=length*ProbBi_dist[ibox];
        prob[2]+=length*ProbKi_dist[ibox];
        prob[3]+=length*ProbDy_dist[ibox];
        ibox++;

    
        length=Interval_overlap(min_ovlp, ((float)ibox-1.)*dx+MT.xmin_real-MT.length[MT.imin_real]/2, ov, dx);
        
    }while (length>0.00000001 /*&& ibox<nbox*/);

    /*normalize*/
    prob[0]/=ov;
    prob[1]/=ov;
    prob[2]/=ov;
    prob[3]/=ov;
   // printf("prob active=%f\n", prob[0]);
   /*if(prob[1]>0)
       printf("iter=%i prob bipolar=%f\n",iter, prob[1]);*/
   // printf("prob kinesin=%f\n", prob[2]);
    
}

void free_ovlps ( float active, float active_ovlps, float bipol, float unipol_dyn, float unipol_kin, float bundl){
    int iseg,jseg,a,b, x,i;
    float sz,foo, ProbDynein=1-ProbKinesin-ProbBipolar-ProbBundling;
    /*Allowed error in the occupancy approx 6*NMT ovlps possible*/
    float delta=0.5/MT.number;
    /*bipol*/

    if (!(iter%wfreq))
    {
    printf("kinesin=%f dynein=%f bipolar=%f\n",unipol_kin/active, unipol_dyn/active, bipol/active);
    }

    if (bipol/active>ProbBipolar+delta) {
       
        x=(int)floor((bipol/active-ProbBipolar)*active_ovlps/delta);
        i=1;
        
        while (i<=x){
           
            sz=1000000;
            for (iseg=1;iseg<=2*MT.number-1;iseg++){
                for(jseg=iseg+1;jseg<=2*MT.number;jseg++) {
                    if (ovlp.type[iseg][jseg]==BIPOLAR){
                      
                        foo=soverlap(iseg,jseg);
                        if (foo<sz) {
                            sz=foo;
                            a=iseg;
                            b=jseg;
                        }
                    }
                }
            }
            ovlp.type[a][b]=ZERO;
            ovlp.type[b][a]=ZERO;
            i++;
        }
    }
    /*Kinesins +end walker*/
    if (unipol_kin/active>ProbKinesin+delta) {
        
        x=(int)floor((unipol_kin/active-ProbKinesin)*active_ovlps/delta);
        i=1;
        while (i<=x){
            
            sz=1000000;
            for (iseg=1;iseg<=2*MT.number-1;iseg++){
                for(jseg=iseg+1;jseg<=2*MT.number;jseg++) {
                    if ((ovlp.type[iseg][jseg]==LEGDOWN || ovlp.type[iseg][jseg]==LEGUP) && ovlp.motor_direction[iseg][jseg]==-1){
                        
                        foo=soverlap(iseg,jseg);
                        if (foo<sz) {
                            sz=foo;
                            a=iseg;
                            b=jseg;
                        }
                    }
                }
            }
            ovlp.type[a][b]=ZERO;
            ovlp.type[b][a]=ZERO;
            i++;
        }
    }
    /*Dynein -end walker*/
    
    if (unipol_dyn/active>ProbDynein+delta) {
        x=(int)floor((unipol_dyn/active-ProbDynein)*active_ovlps/delta);
        i=1;
        while (i<=x){
            
            sz=1000000;
            for (iseg=1;iseg<=2*MT.number-1;iseg++){
                for(jseg=iseg+1;jseg<=2*MT.number;jseg++) {
                    if ((ovlp.type[iseg][jseg]==LEGDOWN || ovlp.type[iseg][jseg]==LEGUP) && ovlp.motor_direction[iseg][jseg]==1){
                        
                        foo=soverlap(iseg,jseg);
                        if (foo<sz) {
                            sz=foo;
                            a=iseg;
                            b=jseg;
                        }
                    }
                }
            }
            ovlp.type[a][b]=ZERO;
            ovlp.type[b][a]=ZERO;
            i++;
        }
    }
    /*bundling prot*/
    if (bundl/active>ProbBundling+delta) {
        x=(int)floor((bundl/active-ProbBundling)*active_ovlps/delta);
        i=1;
        while (i<=x){
            
            sz=1000000;
            for (iseg=1;iseg<=2*MT.number-1;iseg++){
                for(jseg=iseg+1;jseg<=2*MT.number;jseg++) {
                    if (ovlp.type[iseg][jseg]==BUNDLING){
                        
                        foo=soverlap(iseg,jseg);
                        if (foo<sz) {
                            sz=foo;
                            a=iseg;
                            b=jseg;
                        }
                    }
                }
            }
            ovlp.type[a][b]=ZERO;
            ovlp.type[b][a]=ZERO;
            i++;
        }
    }

    
}


void overlapMap()
{
  
  int iMT,jMT,iMTim,jMTim,type,ihigh;
  int iseg,jseg,iseg1,jseg1,iseg2,jseg2,iseg1im,jseg1im,iseg2im,jseg2im;
  float sz,r,r1,r2,r3,r4,tmp,prob1,prob2,prob3;
  float active,notactive,active_ovlps,totvel,bundl,bipol,unipol_dyn, unipol_kin;
  float dx=0.00000001+(MT.xmax_real+MT.length[MT.imax_real]/2-MT.xmin_real+MT.length[MT.imin_real]/2)/nbox;
  float MT_rad;
  float vel0;
  FILE *fp;
  float prob[5];
  
  
  /* prob3=ProbBipolar+ProbTwoMotors+ProbBundling; */

  // If overlaps between parallel MTs are occupied by BIPOLAR motors 
  // or by TWOMOTORS there is no motion between them and the filaments 
  // as well as the motors are stuck. Also, the next thing that is done 
  // in this routine is to calculate the percentage of active overlaps. 
  // Such overlaps will be considered active since they are occupied by motors.
  // However in practice, they are stuck. Thus the estimate of active 
  // overlaps is problematic. 
  // One possibility is to consider parallel overlaps as non-active (ZERO), 
  // and to let it be asigned randomely again in the next iteration. 
  // This possibility treats parallel and anti parallel interactions 
  // differently. I.e., it assumes that motors unbind parallel filaments 
  // with rate equal to the time step while parallel filaments stay connected 
  // as long as the overlap is greater than zero.
  // Another posibility is to shuffle all overlaps at each time step.
 
  /* The following loop shuffles BIPOLAR motors that are stuck 
     in an overlap region of Parallel filaments for more 
     than v0*sz */
 /*
  for (iseg=1;iseg<=2*MT.number-1;iseg++){
      
      for(jseg=iseg+1;jseg<=2*MT.number;jseg++) {
          iMT=segindex(iseg);
          jMT=segindex(jseg);
          sz=soverlap(iseg,jseg);

    
          if ( sz>EPSI && ovlp.type[iseg][jseg]==BIPOLAR && MT.direct[iMT].x*MT.direct[jMT].x>0 &&((float)iter-ovlp.iter0[iseg][jseg])*dt>=sz/vel0_k){
             
                  ovlp.type[iseg][jseg]=ZERO;
                  ovlp.type[jseg][iseg]=ZERO;
          }
      }
  }*/

  /* Calculate the present ovlps */
    active=0.0;
   
    notactive=0.0;
    bipol=0.0;
    unipol_dyn=0.0;
    unipol_kin=0.0;
    bundl=0.0;
    for (iseg=1;iseg<=2*MT.number-1;iseg++)
        for(jseg=iseg+1;jseg<=2*MT.number;jseg++){
            sz=soverlap(iseg,jseg);
            if (sz>EPSI){
                if (ovlp.type[iseg][jseg]==ZERO){
                    notactive++;
                }
                else{
                    active++;
                    if (ovlp.type[iseg][jseg]==BIPOLAR) bipol++;
                    else if ((ovlp.type[iseg][jseg]==LEGUP || ovlp.type[iseg][jseg]==LEGDOWN) && ovlp.motor_direction[iseg][jseg]==1) unipol_dyn++;
                    else if ((ovlp.type[iseg][jseg]==LEGUP || ovlp.type[iseg][jseg]==LEGDOWN) && ovlp.motor_direction[iseg][jseg]==-1) unipol_kin++;
                    else if (ovlp.type[iseg][jseg]==BUNDLING) bundl++;
                    else{
                        printf("critical error in ovlpmap.c!\n");
                        exit(0);
                    }
                }
            }
        }
    active_ovlps=active/(active+notactive);
    /*Delete excess ovlps that do not satisfy user input*/
 
free_ovlps( active, active_ovlps, bipol, unipol_dyn,  unipol_kin, bundl);

  
  /* generate new active crosslinks according to the 
     desired - ProbActive - probability and the designated
     overlap type */


  for (iseg=1;iseg<=2*MT.number-1;iseg++)
    {
      iMT=segindex(iseg);
      for(jseg=iseg+1;jseg<=2*MT.number;jseg++){
          
	    jMT=segindex(jseg);
	/* we add the transient probability that the motors 
	   in the overlap region be active */
	    

	    sz=soverlap(iseg,jseg);
	    type=ovlp.type[iseg][jseg];
	    

 
        if (sz>EPSI && type==ZERO) {
            /*The following calculates a value for the active prob
             and the other probabilities depending on the current ovlp distribution*/
            if (MOTOR==TRUE){
                if (ProbBundling>0) {
                    printf("\n bundling with probdist not implemented!\n");
                    exit(0);
                }
                get_Prob_active( iseg, jseg, prob, dx);
            }
            else{
                prob[0]=ProbActive;
                prob[1]=ProbBipolar;
                prob[2]=ProbKinesin;
                prob[3]=1-ProbBipolar-ProbKinesin-ProbBundling;
                prob[4]=ProbBundling;
                
            }
                /*Check if probabilities are consistent*/
            if ((prob[1]+prob[2]+prob[3]+prob[4])<0.98) {
                printf("iter%i\n",iter);
                printf("act=%f\nbi=%f\nki=%f\ndy=%f bund=%f\n", prob[0], prob[1], prob[2], prob[3],prob[4]);
                exit(0);
            }
            /*Generate some rand*/
            r1=ran1(&idum);
            r2=ran1(&idum);
            r4=ran1(&idum);
            if (r1<prob[0]-active_ovlps){
              if (r2<= prob[1]) /* motors are bipolar */
              {
      
                ovlp.type[iseg][jseg]=BIPOLAR;
                ovlp.type[jseg][iseg]=BIPOLAR;
                ovlp.iter0[iseg][jseg]=(float)iter;
                ovlp.iter0[jseg][iseg]=(float)iter;
                /*only bipolar kinesin*/
                ovlp.motor_direction[iseg][jseg]=-1;
                ovlp.motor_direction[jseg][iseg]=-1;
              }
               
               /*Motors are unipolar plus end directed*/
              else if( r2 <= (prob[2]+prob[1]) && r2>prob[1]){
              /* motors are unipolar plus walkers */
              /* choose their orientation up- or down-wards */
              
                r=ran1(&idum);
                if (r>ProbLegDown){
                    
                    ovlp.type[iseg][jseg]=LEGUP;
                    ovlp.type[jseg][iseg]=LEGUP;
                    ovlp.iter0[iseg][jseg]=(float)iter;
                    ovlp.iter0[jseg][iseg]=(float)iter;
                    ovlp.motor_direction[iseg][jseg]=-1;
                    ovlp.motor_direction[jseg][iseg]=-1;
                }
                else{
                    ovlp.type[iseg][jseg]=LEGDOWN;
                    ovlp.type[jseg][iseg]=LEGDOWN;
                    ovlp.iter0[iseg][jseg]=(float)iter;
                    ovlp.iter0[jseg][iseg]=(float)iter;
                    ovlp.motor_direction[iseg][jseg]=-1;
                    ovlp.motor_direction[jseg][iseg]=-1;
                }
              }
              else if(r2 <= (prob[3]+prob[2]+prob[1]) && r2>(prob[2]+prob[1]))
              /* motors are unipolar minus end walkers */
              /* choose their orientation up- or down-wards */
              {
                  r=ran1(&idum);
                  if (r>ProbLegDown){
                      
                      ovlp.type[iseg][jseg]=LEGUP;
                      ovlp.type[jseg][iseg]=LEGUP;
                      ovlp.iter0[iseg][jseg]=(float)iter;
                      ovlp.iter0[jseg][iseg]=(float)iter;
                      ovlp.motor_direction[iseg][jseg]=1;
                      ovlp.motor_direction[jseg][iseg]=1;
                      
                  }
                  else{
                      ovlp.type[iseg][jseg]=LEGDOWN;
                      ovlp.type[jseg][iseg]=LEGDOWN;
                      ovlp.iter0[iseg][jseg]=(float)iter;
                      ovlp.iter0[jseg][iseg]=(float)iter;
                      ovlp.motor_direction[iseg][jseg]=1;
                      ovlp.motor_direction[jseg][iseg]=1;
                  }
              }
              /*Passive bipolar = Bundling*/
              else if ( r2 <= (prob[4]+prob[3]+prob[2]+prob[1]) &&  r2 > (prob[3]+prob[2]+prob[1])){
                  ovlp.type[iseg][jseg]=BUNDLING;
                  ovlp.type[jseg][iseg]=BUNDLING;
                  ovlp.iter0[iseg][jseg]=(float)iter;
                  ovlp.iter0[jseg][iseg]=(float)iter;
                  /*only bipolar kinesin*/
                  ovlp.motor_direction[iseg][jseg]=-1;
                  ovlp.motor_direction[jseg][iseg]=-1;
              }
              else {
                  /*due to a numerical calc error it might happen that the overall probability 
                   is slightly below 1 so r2 might be above this does ot happen often
                   so for now I just set those ovlps to ZERO as this makes the error hard
                   to find*/
                  char filename[100]="error_log.dat";
                  if ((fp=fopen(filept(filename),"a"))==NULL)
                      printf("COULDNT OPEN FILE: %s",filename);
                  fprintf(fp,"%i active overlap found no motor and therby remains inactive?\n",iter);
                  fprintf(fp,"act=%f\nbi=%f\nki=%f\ndy=%f\n", prob[0], prob[1], prob[2], prob[3]);
                  fprintf(fp,"r2=%f\n\n",r2);
                  fclose(fp);
                  ovlp.type[iseg][jseg]=ZERO;
                  ovlp.type[jseg][iseg]=ZERO;
                  ovlp.motor_direction[iseg][jseg]=0;
                  ovlp.motor_direction[jseg][iseg]=0;
              }
        
            }
		}
          
	    else if (sz<=EPSI)
	      /* non-overlapping segments */
	      {
              ovlp.type[iseg][jseg]=ZERO;
              ovlp.type[jseg][iseg]=ZERO;
              ovlp.motor_direction[iseg][jseg]=0;
              ovlp.motor_direction[jseg][iseg]=0;
	      }

	
      }
    }


  /* Correct for image cases */
  if (REF)
    {
      for (iMT=1;iMT<=MT.number-1;iMT++)
	{
	  iseg1=2*iMT-1;
	  iseg2=2*iMT;
	  for(jMT=iMT+1;jMT<=MT.number;jMT++)
	    {
	      jseg1=2*jMT-1;
	      jseg2=2*jMT;
	      
	      if (seg.length[iseg2]>0 && seg.length[jseg2]>0)
		/* When both MTs cross the reflecting boundary
		   the image motors should be chosen appropriately */
		{
		  if(soverlap(iseg2,jseg2)>0 && soverlap(iseg1,jseg1)>0)
		    {
		      ovlp.type[iseg1][jseg1]=ovlp.type[iseg2][jseg2];
		      ovlp.type[jseg1][iseg1]=ovlp.type[jseg2][iseg2];
                
              ovlp.motor_direction[iseg1][jseg1]=ovlp.motor_direction[iseg2][jseg2];
              ovlp.motor_direction[jseg1][iseg1]=ovlp.motor_direction[jseg2][iseg2];
                
            }
		  
		  if (soverlap(iseg1,jseg2)>0 && soverlap(iseg2,jseg1)>0)
		    {
		      ovlp.type[iseg2][jseg1]=ovlp.type[iseg1][jseg2];
		      ovlp.type[jseg1][iseg2]=ovlp.type[jseg2][iseg1];
                
                ovlp.motor_direction[iseg2][jseg1]=ovlp.motor_direction[iseg1][jseg2];
                ovlp.motor_direction[jseg1][iseg1]=ovlp.motor_direction[jseg2][iseg1];
		    }
		}
	    }
	}
    }
  
  if (twopop==1)
    {
      for (iMT=1;iMT<=MT.number-3;iMT+=2)
	{
	  iseg1=2*iMT-1;
	  iseg2=2*iMT;
	  iMTim=iMT+1;
	  iseg1im=2*iMTim-1;
	  iseg2im=2*iMTim;

	  for(jMT=iMT+2;jMT<=MT.number-1;jMT+=2)
	    {
	      jseg1=2*jMT-1;
	      jseg2=2*jMT;
	      jMTim=jMT+1;
	      jseg1im=2*jMTim-1;
	      jseg2im=2*jMTim;
	      
	      
	      ovlp.type[iseg1][jseg1]=ovlp.type[iseg1im][jseg1im];
	      ovlp.type[jseg1][iseg1]=ovlp.type[jseg1im][iseg1im];
	      ovlp.type[iseg1][jseg1im]=ovlp.type[iseg1im][jseg1];
	      ovlp.type[jseg1im][iseg1]=ovlp.type[jseg1][iseg1im];
            
            
            ovlp.motor_direction[iseg1][jseg1]=ovlp.motor_direction[iseg1im][jseg1im];
            ovlp.motor_direction[jseg1][iseg1]=ovlp.motor_direction[jseg1im][iseg1im];
            ovlp.motor_direction[iseg1][jseg1im]=ovlp.motor_direction[iseg1im][jseg1];
            ovlp.motor_direction[jseg1im][iseg1]=ovlp.motor_direction[jseg1][iseg1im];
	      
	    }
	}
    }



}

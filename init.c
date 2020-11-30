/* #################################################################### */
/* init  --  main Loops                                     */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global_var.h"
#include "nrutil.h"
#include <time.h>

#define EPS 1.0e-10

void init()
{ 
    
   int newMTs,ibox,iMT,N,i;
    clock_t t1, t2, teq=0, ttmp;
    /*If certain features are switched on all averaging calculations
     are useless*/
  int av_flag=1;
  float length,avMTlength;
  long it1,it2,idt,tau;
  DEBUG=1;
    t1=clock();
    
    /*Initialize av structure*/
    initialize_av_variables();
    t2=clock();

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
  /* Ensemble averaging loop */  
  for (iav=1; iav<=nav; iav++)
    {
      //Remove debugging file
      if (remove (filept("Debug_MT_degrade.dat"))==0) printf("Removed old file:%s\n","Debug_MT_degrade.dat");  
        /////////////////////////////////////////////////////////////////
        newMTs=MTfirst; // number of initial MT same every loop
        /*Always reset the index counter to zero! So that at the beginning
         all indexes are "unused"*/
        emptyint(track.index_count,1, niter*mxmts);
        //if (DEBUG==1) printf("1\n"); 
        
        //If we take a diffusing + drifting motor dist into account get array of motor bound density
        //Subsequently fill this array with the initial dist
        dist_dyn=vector(1,nbox);
        dist_kin=vector(1,nbox);
        ProbBi_dist=vector(1,nbox);
        ProbKi_dist=vector(1,nbox);
        ProbDy_dist=vector(1,nbox);
        ProbAct_dist=vector(1,nbox);
        Initial_dist();
        
        /*Initialize iteration structure*/
        /*############################*/
        initialize_it_variables(newMTs);
        /*############################*/
         
        //if (DEBUG==1) printf("2\n"); 
        
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
        //CORE LOOP////
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
    
      /* Start time loop -- trajectory */
      for (iter=1;iter<=niter;iter++)
	{
        
   
	  /* Generate MT array and calculate velocities */
       
        //printf("time between newpos -> mkMT %f \n", (double)(clock()-t)/ CLOCKS_PER_SEC);
        //t=clock();
        ttmp=clock();
        mkMTs(newMTs);  // MT array initialization
        //printf("time between mkMT -> ovlpmap %f \n", (double)(clock()-t)/ CLOCKS_PER_SEC);
        //t=clock();
        teq+=clock()-ttmp;
  
      
        overlapMap();  //get them overlaps

        //printf("time between ovlpmap -> perc %f \n", (double)(clock()-t)/ CLOCKS_PER_SEC);
        //t=clock();
        if (MOTOR==TRUE) {
            prob_update();
        }
        /*If we also simulate actin we need to get an ovlp map for it*/
        if (ACTIN==TRUE) {
            actin_overlapMap();
        }
  
        
        percolation(); //searches for clusters and stores them in MT.cluster array whereas the 0 marks if the system is percolated or not

        //printf("time between perc -> equations %f \n", (double)(clock()-t)/ CLOCKS_PER_SEC);
        
        
        
        if (MT.number>0){
        	equations();   //solve equations
        }
        
        
        //t=clock();
        
        //wdata();       //write the data only debug
      
        update_av();   //update all tracking and av parameters
       // printf("time between update -> newpos %f \n", (double)(clock()-t)/ CLOCKS_PER_SEC);
       // t=clock();

        newpos();// MOVE

        //if (DEBUG==1) printf("3\n"); 
	/* MT flipping if probability arger than 0+EPS*/
	if (pFlp>EPS){
		flipMTs();
	}

        //if (DEBUG==1) printf("4\n"); 
      // The following function applies some potential which forces MTs to move closer to 0 in the yz plane
	  if (MCfreq>0)
	    if ((iter+MCfreq-1)%MCfreq==0) compact();
	  /* if ((iter+MCfreq-1)%MCfreq==0) compactx(); */
	  
	  /* Depending on input, control the entering
	     flux*/
	  if (FixedDensity)
	    {
            printf("FIXED DENS BROKEN!\n");
            exit(0);
	      length=MT.xmax_real-MT.xmin_real;
	      avMTlength=0.5*(MT.maxL+MT.minL);
	      N=(int)floor(0.5+NmaxMTs*(length+avMTlength)/avMTlength);
	      if (MT.number<N)
		newMTs=N-MT.number;
	      else
		newMTs=0;
	      printf("init:> iter=%d, length=%f N=%d MT.number=%d NmaxMTs=%d avMTlength=%f\n",iter,length,N,MT.number,NmaxMTs,avMTlength);
	      printf("init:> iter=%d, adding %d new MTs\n",iter,newMTs);getchar();
	    }
 	  else
            if ((iter+newfreq-1)%newfreq==0){
	      newMTs=newMTs0;
              //MAX 30/10/18 see what happens if newMTs depends on number of MTs
              //printf("%f\n",((float)newMTs0*MT.number)/100.0);
              //newMTs=ceil(((float)newMTs0*MT.number)/100.);
            }
	    else
	      newMTs=0;
        

        
     
        
	}
        
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
        ////////////////////////////////////////////////
        
        /*BEFORE doing the averaging calculations
         Check if the user wants to use the averaging
         iterations for a special plot. Returns 0 or 1 
         PLEASE do not set anything to 1 in this file without knowing exactly what youre doing!*/
      av_flag=special_plots();
        


      

      /* print MT trajectories watch out that it only prints trajectories above neqsteps*/
      //if (iav==1) wav_xtraject();
      
    if (/*iiav==1*/1) {
            if (D3_output==TRUE) wav_traject_3Dmathematica();//mathematica animation output
            else wav_traject_2Dmathematica();
    }
        
    
        /*FREE all MT data*/
        free_MTdat();
        ////////////////////
        
    
      /*Only calculate averages if av_flag gives 1
       There are some features which interpret nav differently
       and therby render averaging useless*/
     if (av_flag){
            /*Normalize all averaged quantities (divide by iav)*/
         
            normalize();
            ///////////////////////////////////////
       
         for (iter=1;iter<=niter;iter++){
             if ((iter+wfreq-1)%wfreq==0){
                 wav_tubulin_dens_it(iter);
                 wav_tubulin_op(iter);
                 wav_avstrands(iter);
                 //wav_motordens(iter);
             }
             
         }
            for (iter=1;iter<=niter;iter++)
             if ((iter+wfreq-1)%wfreq==0) wav_xvelprof(iter);
            for (iter=1;iter<=niter;iter++)
                if ((iter+wfreq-1)%wfreq==0) wav_absveldensity(iter);
            for (iter=1;iter<=niter;iter++)
                if ((iter+wfreq-1)%wfreq==0) wav_veldensity(iter);

         
         
        
            wav_clustervel();
            wav_tubulin_dens();
            wav_end2end();
            wav_avovlpsize_perMT();
            wav_avovlpsize_perOV();
            wav_avstrands_tot();
         
            wav_V();
            wav_totabsveldensity();
            wav_totveldensity();
            wav_xvelprofToT();
        
            
            /*Unnormalize all averaged quantities (multiply by iav)*/
            unnormalize();
        }
    ///////////////////////////////////////
      /*for (iter=1;iter<=niter;iter++)
	for (ibox=0;ibox<nbox;ibox++)
	av.polarcor[iter][ibox] =0.0;  */
        
       free_vector(sum.xl,0,niter);
       free_vector(sum.xl2,0,niter);
       free_vector(sum.xr,0,niter);
       free_vector(sum.xr2,0,niter);
       free_vector(sum.nl,0,niter);
       free_vector(sum.nr,0,niter);
        
    }
  
    /*Free averything for end*/
    free_all();
    printf("process exited in %f sec\n", (float)(clock()-t2)/CLOCKS_PER_SEC);
    printf("process including initialization exited in %f sec\n", (float)(clock()-t1)/CLOCKS_PER_SEC);
     printf("time it took for eq solving %f sec\n", (float)(teq)/CLOCKS_PER_SEC);
    ///////////////////////////
}

  



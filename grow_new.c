/* ######################################################################### */
/* These functions should determine: (i) rate of growth or shrinkage at MT plus ends; (ii) switching between growth and shrinking states*/

#include <stdio.h>
#include "global_var.h"

float rateofgrowth(float growth_rate_at_plus_end_basal){
    int iseg,iMT, MT_growth_direction;
    float growth_rate_at_plus_end_basal; /*some constant - check with Max as to specific value*/
    float axon_tip, soma, axon_length, d_tip, MT_centre, MT_length;

    axon_tip  = rbound0; /*or is it, as in mkMT.c: "rb"?*/
    soma = lbound0; /*or is it, as in mkMT.c: "lb"?*/
    axon_length = axon_tip - soma;

    for (iMT=1; iMT<=MT.number; iMT++) {

      /*Evaluate MT position, distance from axon tip, orientation, and direction of growth*/
      MT_centre = MT.cm[iMT].x; /*in mkMT.c: "xcm"*/
      d_tip = axon_tip - MT_centre;
      MT_orientation = MT.direct[iMT].x;
      MT_growth_direction = MT.grow_direct[iMT];

      /*If MT is growing*/
      if(MT_growth_direction==1){
        /*Growth rate at plus end of plus-end-out microtubules*/
        if(MT_orientation<0.){
          return growth_rate_at_plus_end_basal + exp(d_tip - axon_length);
        }
        /*Growth rate at plus end of minus-end-out microtubules*/
        else {
          return growth_rate_at_plus_end_basal;
        }
      }

      /*If MT is shrinking*/
      else if (MT_growth_direction==-1){
        /*Assume that shrink rate is (a) independent of MT position and (b) equal to basal growth rate*/
          return -growth_rate_at_plus_end_basal;
      }

      /*If MT is idle*/
      else if(MT_growth_direction==0){
        return 0.0;
      }

      else{
          printf("Fatal Error in rateofgrowth function!\n MT growth rate not specified!\n")
        }
    }

}

void switching(){
  int iseg,iMT, MT_growth_direction,rand;
  float MT_centre, MT_plus_end, MT_length, p_switch_to_cat, p_switch_to_res;
  float f_cat = 0.05; /*per second*/
  float f_res = 0.1; /*per second*/
  /*Evaluate probabilties of switching between growth and shrinking states*/
  /*From growth to shrinking*/
  p_switch_to_cat = 1 - exp(-dt*f_cat);
  /*From shrinking to growth*/
  p_switch_to_res = 1 - exp(-dt*f_res);

  /*Now determine switching between growth and shrinking states*/
  for (iMT=1; iMT<=MT.number; iMT++) {
    growth_or_shrink_rate=rateofgrowth(1.0);
    /*Evaluate MT position, length, and direction of growth*/
    MT_centre = MT.cm[iMT].x; /*in mkMT.c: "xcm"*/
    MT_plus_end = MT.pend[iMT].x;
    MT_length = MT.length[iMT].x;
    MT_growth_direction = MT.grow_direct[iMT];
      /*If MT is in growth state*/
      if (MT_growth_direction==1) {
          MT_centre += 0.5*growth_or_shrink_rate*dt;
          MT_plus_end += growth_or_shrink_rate*dt;
          MT_length += growth_or_shrink_rate*dt;
      if (ran1(&idum) < p_switch_to_cat){
              MT_growth_direction=-1;
          }
      }
      /*If MT is in shrinking state*/
      else if (MT_growth_direction==-1) {
          MT_centre += 0.5*growth_or_shrink_rate*dt;
          MT_plus_end += growth_or_shrink_rate*dt;
          MT_length += growth_or_shrink_rate*dt;
          if (ran1(&idum) < p_switch_to_res) {
              MT_growth_direction=1;
          }
      }
      /*If MT is idle*/
      else if (MT_growth_direction==0) {
          if (ran1(&idum) < p_switch_to_cat) {
              MT_growth_direction =-1;
          }
          else{
              MT_growth_direction=1;
          }
      }

      else{
          printf("Fatal Error in switching function!\n MT switching not specified!\n")
      }
  }

}

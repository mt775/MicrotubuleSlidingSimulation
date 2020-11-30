#include <stdio.h>
#include <math.h>
#include "sampleExpDist.h"
#include "ran1.h"
#include "lambert_w.h"
#include "global_var.h"

/*Generate inverse CDF*/
double sampleExpDist(float L){
  float x = ran1(&idum);
  float b = 0.038;
  float a = (1/L)*exp(-b*L);
  double z = (double) (1/a)*(b*exp((-1/a)*b*(x-1)));
  return (1/(a*b))*(-b+(b*x)+(a*LambertW(z)));
}

/* A simple solution to the diffusion advection equation */
/* C. B. Connor
 April 1, 2003
 Department of Geology
 University of South Florida
 Description: The C code calculates concentration as a function of
 time for given boundary conditions and a fixed velocity field using
 an explicit finite difference solution to the diffusion - advection
 equation. This solution is quite prone to instabilities!
 
 Input: the files vel_u.dat and vel_v.dat must be
 present in your directory for this program to function
 
 Usage: Compile the program and run by specifying the total number
 of time steps on the command line
 */

#include <stdio.h>
#include <stdlib.h>

#define KAPPA 0.0001 /*kappa is the diffusion coefficient - m*m/s */
#define DELT 10.0 /* the time step - s */
#define DELX 0.3 /* the distance step  - m */


int main (int argc, char *argv[]) {
    
    /*initialize the variables used in this code */
    FILE *ufile;
    FILE *vfile;
    int i, j, x,y,t; /*counters */
    int imax=50; /*max cells in x-direction */
    int jmax=50; /*max cells in y-direction */
    int tmax; /*max time steps - from command line */
    double u_term, v_term, laplace_term; /*for calculations */
    double part1; /*for calculations */
    
    double u[50][50]; /* the u component of the velocity field - m/s */
    double v[50][50]; /* the v compoent of the velocity field - m/s */
    double c[50][50], c2[50][50]; /* concentration - say in ppm */
    
    /* make sure the correct number of arguments are supplied
     on the command line */
    if (argc != 2) {
        fprintf(stderr, "Usage: <filename> <total number time steps>\n");
        return (0);
    }
    
    
    /*open the files containing the components of the velocity
     vector field  u,v */
    ufile = fopen("vel_u.dat", "r");
    vfile = fopen("vel_v.dat", "r");
    
    /*convert the command line argument (total number
     of time steps) to an integer */
    tmax = atoi(argv[1]);
    
    /* read the velocity field from files into the arrays
     for the arrays u[][] and v[][] - careful with those
     subscripts...*/
    for (i=0; i<imax; i++) {
        for (j=0;j<jmax;j++){
            fscanf(ufile,"%d %d %lf", &x,&y,&u[j][i]);
            fscanf(vfile,"%d %d %lf", &x,&y,&v[j][i]);
            
            c[i][j]=0.0; /*initialize the arrays that track concentration */
            c2[i][j] = 0.0;
        }
    }
    
    /*create the initial plume source */
    c[25][10]=10.0;
    
    /*iterate through time */
    for (t=1;t<tmax; t++) {
        /*iterate across the grid area - not including the
         boundary elements */
        for (i=1; i<imax-1; i++) {
            for (j=1;j<jmax-1;j++){
                
                /*solve the diffusion - advection equation */
                u_term=u[i][j]*(c[i][j]-c[i-1][j])/DELX;
                v_term=v[i][j]*(c[i][j]-c[i][j-1])/DELX;
                laplace_term = KAPPA/(DELX*DELX)*(c[i-1][j]+c[i+1][j]+c[i][j-1]+c[i][j+1]-4.0*c[i][j]);
                
                part1=laplace_term + u_term + v_term;
                c2[i][j] = DELT*part1 + c[i][j];
                
            }
        } /*loop through i and j */
        
        /*update the c[][] array with the new values */
        for (i=1; i<imax-1; i++) {
            for (j=1;j<jmax-1;j++){
                c[i][j]=c2[i][j];
            }
        }
        
    } /* loop through the time step */
    
    /*print the results to standard output */
    for (i=0; i<imax; i++) {
        for (j=0;j<jmax;j++){
            printf("%d %d %lf\n", i,j,c[i][j]);
        }
    }
    
    
}

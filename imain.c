/* #################################################################### */
/* Main Program - Dynamics of Motor Proteins-Microtubule Assembly       */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "global_var.h" 

int main(int argc, char *argv[]){
    
    int i, dir_miss=1;
    lin_density=0; //default
    pFlp=0;
    forcepol=0;
    forcelength=0;
/* Get input/output directory name from command line */
    i=1;
    while (i<argc){
        if ((strcmp(argv[i],"-outdir"))==0) {
            dirname=(char *)malloc((strlen(argv[i+1])+2)*sizeof(char));
            strcpy(dirname,argv[i+1]);
            strcat(dirname,"/");
            dir_miss=0;
            i++;
        }
        
        else if ((strcmp(argv[i],"-lindist"))==0) {
            lin_density=1;
            printf("lindens on!\n");
        }
        else if ((strcmp(argv[i],"-forcepol"))==0) {
            forcepol=1;
            printf("Averaging iterations are interpreted as raising polarity ratio for each calculation\n");
            printf("IMPORTANT: This needs all averaging calculations switched of!\n");
        }
        else if ((strcmp(argv[i],"-forcelength"))==0) {
            forcelength=1;
            printf("Averaging iterations are interpreted as raising externalforce for each calculation\n");
            printf("IMPORTANT: This needs all averaging calculations switched of!\n");
        }
        else usage_err(argv);
        i++; 
    }
    if (dir_miss){
        printf("No output directory specified....\n ");
        dirname=(char *)malloc((8)*sizeof(char));
        strcat(dirname,"OCHECK/");
        printf("Setting %s as output directory\n", dirname);
        
    }

    rdata(argc,argv); // read data from rdat file
    init();
    return(0);

}

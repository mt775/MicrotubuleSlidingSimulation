/* ################################################################# */
/* Functions for memory allocation                                   */

#include <stdlib.h>

/*---------------------------------------------------------------------*/
/* Float 1D */
float *alocF1D(int nx)
{
  float *arrpt;

  arrpt=(float *)calloc(nx,sizeof(float));
  if (arrpt==NULL) 
    return (NULL);
  else
    return(arrpt);
}

/*---------------------------------------------------------------------*/
/* Float 2D */
 
float **alocF2D(int nx, int ny)
{
  int ix;
  float **arrpt;
  
  arrpt=(float **)calloc(nx,sizeof(float *));
  if (arrpt==NULL) 
    return (NULL);
  else
    for (ix=0;ix<nx;ix++)
      {
	arrpt[ix]=(float *)calloc(ny,sizeof(float));
	if (arrpt[ix]==NULL)
	  return (NULL);
      }
  return(arrpt);
}

/*---------------------------------------------------------------------*/
/* Float 3D */
 
float ***alocF3D(int nx, int ny, int nz)
{
  int ix,iy;
  float ***arrpt;
  
  arrpt=(float ***)calloc(nx,sizeof(float **));
  if (arrpt==NULL) 
    return (NULL);
  else
    for (ix=0;ix<nx;ix++)
      {
	arrpt[ix]=(float **)calloc(ny,sizeof(float *));
	if (arrpt[ix]==NULL)
	  return (NULL);
	else
	  for (iy=0;iy<ny;iy++)
	    {
	      arrpt[ix][iy]=(float *)calloc(nz,sizeof(float));
	      if (arrpt[ix][iy]==NULL)
		return (NULL);
	    }
      }
  return(arrpt);
}


/*---------------------------------------------------------------------*/
/* INT 1D */
int *alocI1D(int nx)
{
  int *arrpt;

  arrpt=(int *)calloc(nx,sizeof(int));
  if (arrpt==NULL) 
    return (NULL);
  else
    return(arrpt);
}

/*---------------------------------------------------------------------*/
/* Int 2D */
 
int **alocI2D(int nx, int ny)
{
  int ix;
  int **arrpt;
  
  arrpt=(int **)calloc(nx,sizeof(int *));
  if (arrpt==NULL) 
    return (NULL);
  else
    for (ix=0;ix<nx;ix++)
      {
	arrpt[ix]=(int *)calloc(ny,sizeof(int));
	if (arrpt[ix]==NULL)
	  return (NULL);
      }
  return(arrpt);
}

/*---------------------------------------------------------------------*/
/* Int 3D */
 
int ***alocI3D(int nx, int ny, int nz)
{
  int ix,iy;
  int ***arrpt;
  
  arrpt=(int ***)calloc(nx,sizeof(int **));
  if (arrpt==NULL) 
    return (NULL);
  else
    for (ix=0;ix<nx;ix++)
      {
	arrpt[ix]=(int **)calloc(ny,sizeof(int *));
	if (arrpt[ix]==NULL)
	  return (NULL);
	else
	  for (iy=0;iy<ny;iy++)
	    {
	      arrpt[ix][iy]=(int *)calloc(nz,sizeof(int));
	      if (arrpt[ix][iy]==NULL)
		return (NULL);
	    }
      }
  return(arrpt);
}

/*---------------------------------------------------------------------*/
/* Free 2D matrix */

void freeF2D(float **matrix, long int nx)
{
  int i,j;
  for (i=0;i<nx;i++)
    free(matrix[i]);
  free(matrix);
}

/*---------------------------------------------------------------------*/
void free2D(void **matrix, long int nx)
{
  int i,j;
  for (i=0;i<nx;i++)
    free(matrix[i]);
  free(matrix);
}


/*---------------------------------------------------------------------*/
/* Free 2D matrix */

void freeF3D(float ***matrix, long int nx, long int ny)
{
  int i,j;
  for (j=0;j<nx;j++)
    for (i=0;i<ny;i++)
      free(matrix[j][i]);
  for (j=0;j<nx;j++)
    free(matrix[j]);
  free(matrix);
}


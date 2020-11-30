#include <stdio.h>
#include <math.h>
#include "global_var.h"
#include "nrutil.h"
#include <stdlib.h>

void group_search();
void group_search_n(int *arr_i, int n) ;
int CountConnectedGroups(int *arr_i);
void PrintConnectedGroups(int *arr_i, int display_count);

int **makepercolationArr(int *arr_i);
int Compare_groups(int* group1,int *group2);
int in_same_group(int i1, int i2, int *group);
int group_no(int MT_index, int **percolationGrp);
int* Select_group(int MT_Index,int *group);
void print_group(int *group);

void move_MT(int i,float dx);
void Move_group(int* group,float dx);
int* find_RightMT();
int* find_LeftMT();
float find_center();
void move_to_center();

int MT_to_Seg_ovlp(int i, int j);
int is_in_arr (int n,int *arr_i);
void freepercolationArr(int **perculation_groups);


void percolation()
{
    if (MT.number>0) {

        int *R_MT = find_RightMT(); // a vector contain max x MTs
        int *L_MT = find_LeftMT();  // a vector contain min x MTs
   

        int *group;
        group = ivector(0,2*MT.number);
        //group =  (int*) calloc (2 * MT.number , sizeof(int) );//hint: calloc sets the allocated mem to zero
        emptyint(group,0,2*MT.number); // one should never trust
        group_search(group); //a function returning an array with the indexes of the connected groups in a row and zeros
                            //in between groups
        
        int **percolation_groups = makepercolationArr(group);

        // check if there is any percolation from Left to Right
        int count_r,count_l;
        count_r=count_l=0;
        while(R_MT[count_r]!=0)count_r ++;
        while(L_MT[count_l]!=0) count_l++;
        
        int r,l;
        int perculate = 0;
        for(r=0;r<count_r;r++)
        for(l=0;l<count_l;l++)
          /*checks for each MT in l/r boundary region if in same group
           also passes group vector*/
          if(in_same_group(R_MT[r],L_MT[l],group))
             perculate = group_no(R_MT[r],percolation_groups);
      

        if (perculate)
        percolation_groups [0][1]=perculate;

        free(R_MT);
        free(L_MT);
        free_ivector(group,0,2*MT.number);

        /* if (1)
        { printf("printing percolation_groups\n");
          printpercolationArr(percolation_groups); }

        if (DEBUG_PERCOLATION)
        printf("copying percolation_groups to MT.percolation\n");
        copypercolationArr(percolation_groups,&(MT.percolation) );
        if (DEBUG_PERCOLATION)
        { printf("printing MT.percolation\n");
          printpercolationArr( MT.percolation ); }*/

        /*
         [0][0] - the no of the groups
         [0][1] - is 0 if the bundle is not perculated, and 1 for perculated bundle;
         [1][0] - the no of elements in group #1
         [1][1...n] - the elements in group #1
         [2][0] - the no of elements in group #2
         [2][1...n] - the elements in group #2
         .
         .
         for example:
         the percolated vector  1 2 5 6 0 3 4 8
         will be arranges to arr:
         2 1
         4 1 2 5 6
         3 3 4 8
         */

        // make a vector MT.cluster.
        // every MT will be indicated by the index no of it's cluster.
        // data is saved to MT.cluster[iMT]
        int i,j;
        for (i=1;i<=percolation_groups[0][0];i++){
            for (j=1; j<=percolation_groups[i][0];j++){
                /*gives for each MT the number of total MTs in his group*/
                MT.cluster[ percolation_groups[i][j] ]= percolation_groups[i][0];

            }
        }

        /*for (i=1;i<=percolation_groups[0][0];i++)
            for (j=1; j<=percolation_groups[i][0];j++){

                MT.cluster[ percolation_groups[i][j]  ] = i;
            }*/
        /*Save percolation switch in zero index of MT.cluster*/
        MT.cluster[0]= perculate;
        /*update percolation probability*/
        if (perculate>0)
           perc_prob+=1.;
   
        freepercolationArr(percolation_groups);

    }
}


//###########################################################//

void group_search(int *group)
{ 
  group_search_n(group, 1);
}

//--------------------------------------------------
void group_search_n(int *arr_i, int n)
{
  int i,j;
  static int index;
  
  static int tree;
  // if (n>20)printf("in group search n=%d iter=%d\n",n,iter);
    
  //For the first call sets statics to zero and generates array of lenght 2*MT.number
  if (n == 1 )
    {
      //if (iter ==1) print_ovlp();
      
      index = 0;
      tree = 0;
    }
  ///////////////////////////////////
    
  //If this function went through ALL MTs it returns arr_i to the group_search functions
  // and thereby to the group array in percolation
  if (n > MT.number)
    {
      //PrintConnectedGroups(arr_i, 1);
      //free(arr_i);
      return;
    }
      
  
  // UPDATING GROUPS VECTOR.
    //printf("%i\n",n);
  arr_i[index]=n;
  index++;
  
  // FINDING ALL THE BRANCHES FROM THE ROOT n.
  //
  for (i=1;i<=MT.number;i++)
    {
      
      if (  (i != n)
	    && ( MT_to_Seg_ovlp(i,n) ) 
	    &&  !(is_in_arr(i,arr_i)) )
	
	{
	  tree ++ ;
	  group_search_n(arr_i, i);// in this function tree is NEVER equal to 0 so the cell below is not executed
	  tree --;
	}
    }


  //RECALLING THE FUNCTION IN NEW ROOT n+1.  
  if (tree == 0)
    {
      n++; index++;
      while ( (is_in_arr(n,arr_i))  && (n <= MT.number) ) n++;
      group_search_n( arr_i, n);
    }
}



// ++++++++++++++MT_to_Seg_ovlp+++++++++++++++++
// gets two MT's index, and returns ovlp.type between 
// them in terms of segmants.
//+++++++++++++++++++++++++++++++++++++++++++++
int MT_to_Seg_ovlp(int i, int j)
{
  if ( ovlp.type[2*i-1][2*j-1] ==LEGDOWN || ovlp.type[2*i-1][2*j-1] ==LEGUP)
      return  ( ovlp.type[2*i-1][2*j-1] );
  if ( ovlp.type[2*i-1][2*j] ==LEGDOWN || ovlp.type[2*i-1][2*j] ==LEGUP )
      return  ( ovlp.type[2*i-1][2*j]   );
  if ( ovlp.type[2*i]  [2*j-1] ==LEGDOWN || ovlp.type[2*i]  [2*j-1] ==LEGUP)
      return  ( ovlp.type[2*i]  [2*j-1] );
  if ( ovlp.type[2*i-1][2*j-1] ==LEGDOWN || ovlp.type[2*i-1][2*j-1] ==LEGUP )
      return  ( ovlp.type[2*i-1][2*j-1] );
  
  return 0;
}




void freepercolationArr(int **perculation_groups)
{

    //   printPerculationArr(perculation_groups);
    int i,j;
    for (i=perculation_groups[0][0];i>=0;i--)
    {
        free(perculation_groups[i]);
    }
    free(perculation_groups);
}


//++++++++++++++++++++++a bulean function++++++++++++++++++++++++//
//+++ checks wether n is already located in the groups vector+++//
int is_in_arr (int n,int *arr_i)
{
  int j=0;
  while (arr_i[j]+arr_i[j+1])
    {
      if (arr_i[j]==n)  return 1;
      j++;
    }
  return 0;
  
}





//+++++++++++++++printing the overlap map to the screen++++++++++//
void print_ovlp()
{
  int i,j;
  printf("printing ovlp.map\n iter no %d\n",iter);
  
  // printing the decades headline  
  printf("\t\t\b| "); 
  for (i=1; i<MT.number;i++)
    if ( (i%10) == 0) 
      printf(" |%d \b",10*(int)(i/10));
    else 
      printf("  ");
  printf("\n");
  
  //printing the ones headline
  printf("\t\\jMT no\t\b|");
  for (i=1;i<=MT.number;i++)
    if ( (i%10) != 0) 
      printf("%d ",i%10);
    else  
      printf("10| ");
  printf("\n");

  
  //printing rows values  &  seperators : '|'
  for (i=1 ;i<=MT.number;i++)
    {
      if ((i%10)==1)
	{printf("\t\t\b+");
	  for(j=1;j<=MT.number;j++)
	    {
	      printf("--");
	      if ((j%10)==0) printf("+-");
	    }
	  printf("\n");
	}
      printf("iMT no %d -\t\b|",i);
      for(j=1;j<=MT.number;j++)
	{
	  printf("%d ", MT_to_Seg_ovlp(i,j) );
	  if ((j%10)==0) printf("| ");
	}
      printf("\n");
    }
  

  //getchar();
}


//+++++++++++++++++++++PrintConnectedGroups+++++++++++++++++++
//+++ printing the vector of the conected groups.
//++  if display_count==1
//+   calculates the number of sepearated groups, and prints to the screen.
void PrintConnectedGroups(int *arr_i, int display_count)
{
  int count;
  int i;
  //k will point to the last occupied site in arr_i vector 
  int k = 0;
  while((arr_i[k]+arr_i[k+1])!=0) k++;
  
  // printf("\niter no %d\n",iter);
  // printf ("CONNECTED GROUPS ARE: \n");
  
  printf ("{");
  for(i=0;i<k;i++)
    if (arr_i[i]) 
      printf("%d, ",arr_i[i]);
    else
      {
	printf("\b\b}  ; {");
	//	j++;
      }
  printf("\b\b}\n");
  
  if (display_count)
    {
      count = CountConnectedGroups(arr_i) ;
      if (count==1)
	printf("All MT's arranged in one connected group\n");
      else
	printf("there are %d seperate groups\n",count);
    }
  //getchar();
    
}

int CountConnectedGroups(int *arr_i)
{
  int j=1;
  int i;
  //k will point to the last occupied site in arr_i vector 
  int k = MT.number;
  while((arr_i[k]+arr_i[k+1])!=0) k++;

 for(i=0;i<k;i++)
    if (!arr_i[i]) 
      j++;
      
 return j;      
  
}



//#############################################//
int *find_RightMT() //returns vector or the max x MT's
{
  int i,j,count;
  float MT_right_edge =0, Bundle_Right_edge=lbound0;
  for (i=1;i<=MT.number*2;i++) {
      MT_right_edge = seg.cm[i].x + seg.length[i]/2;
      if  (MT_right_edge>Bundle_Right_edge) {
          Bundle_Right_edge = MT_right_edge;
      }
  }
  
  count = 0;
  for (i=1;i<=MT.number*2;i++)
    {
      MT_right_edge = seg.cm[i].x + seg.length[i]/2;
      if  (MT_right_edge==Bundle_Right_edge)
	count++;
    }

  int* tmp_arr= (int*) calloc (count+1 , sizeof(int) );
  j=0;
  for (i=1;i<=MT.number*2;i++)
    {
      MT_right_edge = seg.cm[i].x + seg.length[i]/2;
      if  ( (MT_right_edge==Bundle_Right_edge) && 
	    ( (j==0) || ((i+1)/2)!=tmp_arr[j-1]) )
	{
	  tmp_arr[j]=(i+1)/2;
	  j++;
	}
    }
  
  return tmp_arr;


}

//--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--//

int *find_LeftMT()//returns vector or the min x MT's
{
  int i,j,count;
  float MT_Left_edge =0, Bundle_Left_edge=rbound0;
  for (i=1;i<=MT.number*2;i++)
    {
      MT_Left_edge = seg.cm[i].x - seg.length[i]/2;
      if  (MT_Left_edge<Bundle_Left_edge) 
	{
	  Bundle_Left_edge = MT_Left_edge;
	}
    }

  count = 0;
  for (i=1;i<=MT.number*2;i++)
    {
      MT_Left_edge = seg.cm[i].x - seg.length[i]/2;
      if  (MT_Left_edge==Bundle_Left_edge) 
	count++;
    }


  int* tmp_arr= (int*) calloc (count+1 , sizeof(int) );
  j=0;
  for (i=1;i<=MT.number*2;i++)
    {
      MT_Left_edge = seg.cm[i].x - seg.length[i]/2;
      if  ( (MT_Left_edge == Bundle_Left_edge) &&
	    ( (j==0) || ((i+1)/2)!=tmp_arr[j-1]) )
	{
	  tmp_arr[j]=(i+1)/2;
	  j++;
	}
    }

  
  return tmp_arr;
  
}



//##############################################
/*Returns the array of indexes contained in MT_indexes group*/
int* Select_group(int MT_Index,int *group)
{
  int i,j,k,l;
  int* tmp_group;
  int is_in_group = 0;
 //tmp_group = (int*) calloc (3 * MT.number , sizeof(int) );
  tmp_group = (int*) calloc (MT.number+1 , sizeof(int) );
  
  i=j=k=l=0;

  while( (k+j) <= 2*MT.number) {
      is_in_group = 0;
      j=0;
      while(group[k+j]){
          if (group[k+j] == MT_Index )
              is_in_group =1;
          j++;
      }
      if (is_in_group)
          for (l=0;l<j;l++){
              tmp_group[i]=group[k+l];
              i++;
          }
      j++;
      k=k+j;
  }

   
  return tmp_group;	

}
	

//####################################################//
void Move_group(int* group,float dx)
{
  int i;

  for (i=1;i<=MT.number;i++)
    if (is_in_arr(i,group))
      move_MT(i,dx);
}

void move_MT(int i,float dx)
{
  MT.cm[i].x = MT.cm[i].x +dx;
  mkseg(i);
}

//##########################################################
/*Takes specific group vector and checks it i1 and i2 are in the same group*/
int in_same_group(int i1, int i2, int *group)
{
  int *group1 = Select_group(i1,group);
  int *group2 = Select_group(i2,group);
  return Compare_groups(group1,group2);
}




//#########################################################
int Compare_groups(int* group1,int *group2)
{
  int i,k1,k2;
  k1=k2=0;

  while(group1[k1]!=0) k1++;
  while(group2[k2]!=0) k2++;
  
  if (k1 != k2){
      free(group1);
      free(group2);
      return 0;
  }
   
  for (i=0;i<k1;i++)
      if (group1[i] != group2[i]){
          free(group1);
          free(group2);
          return 0;
      }
   
  free(group1);
  free(group2);
  return 1;
}


//######################################
void print_group(int *group)
{
  int i = 0;   
  while (group[i])
    {    
      printf("%i ",group[i]);
      i++;
    }
  printf("\n");
}
 
 
//########################################################
  
float find_center()
{
  int i;
  float total_length=0;
  float accumulated_cm=0;
  float average_cm =0;

  for (i=1;i<=MT.number;i++)
    {
      accumulated_cm += MT.cm[i].x * MT.length[i];
      total_length += MT.length[i];
    }
  average_cm =accumulated_cm/ total_length; 
  //printf("CENTER IN FUNCTION =%f\n",average_cm);
  return average_cm;
}

//#################################################//
void move_to_center()
{
  if (!(bound.type[0]== 0 && bound.type[1] == 0)) return;
  float center = find_center();
  int i;
 
 for (i=1;i<=MT.number;i++)
   move_MT(i,-center);
}




 


// change the percolation vector to percolation arr:
/*
[0][0] - the no of the groups
[0][1] - is 0 if the bundle is not perculated, and 1 for perculated bundle;
[1][0] - the no of elements in group #1
[1][1...n] - the element indexes in group #1
[2][0] - the no of elements in group #2
[2][1...n] - the element indexes in group #2
 .
 .
 for example:
 the percolated vector  1 2 5 6 0 3 4 8 
 will be arranges to arr:
 2 1
 4 1 2 5 6 
 3 3 4 8 
*/
int **makepercolationArr(int *arr_i)
{
  int **percolation_groups;
  int count_groups=0;
  int count;
  int i,j,l;
  //k will point to the last occupied site in arr_i vector 
  int k = 0;
  while((arr_i[k]+arr_i[k+1])!=0) k++;

  count = CountConnectedGroups(arr_i) ;
  percolation_groups =        calloc(sizeof(int *),count+1) ;
  percolation_groups[0] =     calloc(sizeof(int),2)         ;
  percolation_groups[0][0] =  count                         ;

  i=0;
  while (i<=k)
    {
      count_groups ++;
      count = 0;
      while (arr_i[i]) {i++;count++;};
      percolation_groups[count_groups]    = calloc(sizeof(int *),count+1);
      percolation_groups[count_groups][0] = count;
      //  printf("%d\t",count);
      l=i;
      // copy form the vector "arr_i" to  matrix "percolation_groups";
      for (i=i-count,j=1;i<l;i++,j++)
	 percolation_groups[count_groups][j]= arr_i[i];
	
      i++; j++;
    }

  return percolation_groups;
  
}


// returns the group index of the MT "MT_index", as is saved in percolationGrp
int group_no(int MT_index, int **percolationGrp)
{
  /*
    [0][0] - the no of the groups
    [0][1] - is 0 if the bundle is not perculated, and 1 for perculated bundle;
    [1][0] - the no of elements in group #1
    [1][1...n] - the elements in group #1
    [2][0] - the no of elements in group #2
    [2][1...n] - the elements in group #2
    .
    .
    for example:
    the percolated vector  1 2 5 6 0 3 4 8 
    will be arranges to arr:
    2 1
    4 1 2 5 6 
    3 3 4 8 
  */
  int i,j;
  for (i=1;i<=percolationGrp[0][0];i++)
    for (j=1;j<=percolationGrp[i][0];j++)
      if (percolationGrp[i][j] == MT_index)
	return i;


  return 0;

}

/*
* Author: * Efthymios Grigorakis - 9694

 * 2021 Aristotle University of Thessaloniki
 * Parallel and Distributed Systems.
*/


#include <stdlib.h>
#include <stdio.h>
//gcc -lm(link math.h)
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "mpi.h"

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
//http://www.stat.cmu.edu/~ryantibs/median/quickselect.c
double quickselect(double *arr, int n, int k) {
  unsigned long i,ir,j,l,mid;
  double a,temp;

  l=0;
  ir=n-1;
  for(;;) {
    if (ir <= l+1) { 
      if (ir == l+1 && arr[ir] < arr[l]) {
	SWAP(arr[l],arr[ir]);
      }
      return arr[k];
    }
    else {
      mid=(l+ir) >> 1; 
      SWAP(arr[mid],arr[l+1]);
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir]);
      }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir]);
      }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1; 
      j=ir;
      a=arr[l+1]; 
      for (;;) { 
	do i++; while (arr[i] < a); 
	do j--; while (arr[j] > a); 
	if (j < i) break; 
	SWAP(arr[i],arr[j]);
      } 
      arr[l+1]=arr[j]; 
      arr[j]=a;
      if (j >= k) ir=j-1; 
      if (j <= k) l=i;
    }
  }
}

typedef struct{
    double *coordinates ;
} point;



void printTable(double **table, int rows,int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			printf(" %f ", table[i][j]);
		}
		printf("\n");
	}

	printf("\n");
}

void printpoints(point *points,int numberofpoints,int dimension,int operation){
    
    for(int i = 0 ; i <numberofpoints ; i++){
        printf("%d Point from operation %d: ",i,operation);
        for(int j = 0 ; j <dimension; j++){
            printf("%f ", points[i].coordinates[j]);
        }
        printf("\n");
    }
    printf("\n");
}

point *makepoints(double **pointsmatrix,int numberofpoints,int dimension){

    point *points = (point *) malloc(numberofpoints * sizeof(point));

    for(int i = 0 ; i < numberofpoints ; i++){
        points[i].coordinates = (double *) malloc(dimension * sizeof(double));
        for(int j = 0 ; j < dimension ; j ++){
            points[i].coordinates[j] = pointsmatrix[j][i] ; 
        } 
    }

    return points;

}

double *distacefrompivot(point pivot,  int numberofpoints, int dimension, point *points){
    double *distances = (double *) malloc(numberofpoints * sizeof(double));

    for(int i = 0 ; i < numberofpoints ; i++){
        double distance = 0 ;
        for(int j = 0  ; j < dimension ; j++){
            distance = pow((points[i].coordinates[j] - pivot.coordinates[j]),2) ;

        }
        distances[i] = distance ;
    }

    for(int i = 0 ; i < numberofpoints ; i++){
        printf("distace from %d is %f \n" , i , distances[i]);
    }

    return distances;

}

int smallerthanmedian(double median,point *points,double *distances,int numberofpoints){
    int smaller = 0 ; 
    // for topically function
    for(int i = 0 ; i < numberofpoints ;i++ ){
       if(distances[i] <median){
           smaller ++;
       }
    }
    /*
    * for mpi function 
    for(int i = 0 ; i < numberofpoints ; i++){
        if(distances[i]>median){
            for(int j = i ; j < numberofpoints ; j++){
                if(distances[i] < median){
                    smaller++;
                    point tmp ; 
                    tmp = points[i];
                    points[i] = points[j];
                    points[j] = tmp ;
                }
            }
        }
    }
    *
    */
    printf("SMALLER THAN MEDIAN %d \n",smaller);
    return smaller;
}

point **partionpoints(point *pointsmatrix,int numberofpoints,int numberofoperations){
    int pointsperoperation = numberofpoints / numberofoperations ;
    point **partpoints = (point **) malloc(numberofoperations * sizeof(point*));
    
    for(int i = 0 ; i < numberofoperations ; i++){
        partpoints[i] = (point *) malloc(pointsperoperation * sizeof(point));
    }

    int currentpoint = 0 ; 
    for(int j = 0 ; j  < numberofoperations ; j++){
        for(int k = 0 ; k < pointsperoperation ; k++){
            partpoints[j][k] = pointsmatrix[currentpoint];
            currentpoint++;
        }
    }


    return partpoints ;
}

double** readbinary(char *filename,int numberofpoints,int dimension){
    FILE* file = fopen(filename,"rb");
    double buffer[numberofpoints];

    double **points = (double **) malloc(dimension * sizeof(double*));
    for(int i = 0 ; i < dimension ; i++){
        points[i] = (double *) malloc(numberofpoints * sizeof(double));
    }

    for(int i = 0 ; i < dimension ; i++){
        fread(buffer,sizeof(double),numberofpoints,file);
        for(int j = 0 ; j < numberofpoints ; j++){
            points[i][j] = buffer[j];
        }
        //need to change line in the here 
        puts("");
    }

    return points ;

}

// 0 -> train_image
// 1 -> test image 
int main(int argc, char **argv){

    //readbinary("mnist.txt",8,4);
    
    int pivot ;
    int pid,nproc;
    int numberofpoints = 8 ;
    int dimension = 4 ;
    double operationpoint[numberofpoints];
    
    MPI_Status status;
   

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int nprocs = nproc;

    
    if(pid == 0 ){
       double **points = (double **)malloc(dimension * sizeof(double*));
        for(int i = 0 ; i < dimension ; i++){
            points[i] = (double *) malloc(numberofpoints * sizeof(double)); 
        }

        for(int i = 0 ; i < dimension ; i++){
            for(int j = 0 ; j < numberofpoints ; j++){
                points[i][j] = (double)rand()/((double)RAND_MAX) ; 
           
            }
        }
        printTable(points,dimension,numberofpoints);
        int pivot = rand() % (numberofpoints/nprocs+1) ;

        for(int i = 1 ; i < nproc ; i ++){
            MPI_Send(&pivot,1,MPI_INT,i,55,MPI_COMM_WORLD);
        }
    
    }
    else{
        MPI_Recv(&pivot, 1, MPI_INT, MPI_ANY_SOURCE, 55, MPI_COMM_WORLD, &status);
        printf("PIVOT %d FROM OPERATION %d ",pivot,pid);
       }
    // printTable(points,dimension,numberofpoints);
    
    
 

    
    // int pivotnumber = rand() % (numberofpoints/nprocs+1);
    //     point pivot = operationpoints[0][pivotnumber] ;
    //     printf("THE PIVOT IS : %d \n" , pivotnumber);
    // point pivot = operationpoints[0][pivotnumber] ;
   /* point *points1 = operationpoints[0];
    point *points2 = operationpoints[1];

    printpoints(points1,numberofpoints/nprocs,dimension,0);
    printpoints(points2,numberofpoints/nprocs,dimension,1);
   

     

     
    
    printf("\n");

    double *distances = distacefrompivot(pivot,numberofpoints,dimension,finalpoints);
    printf("\n");

    double **partdistances = (double **) malloc(nprocs * sizeof(double*));
    for(int i = 0 ; i < nprocs ; i++){
        partdistances[i] = (double *) malloc((numberofpoints/nprocs) * sizeof(double));
        partdistances[i] = distacefrompivot(pivot,numberofpoints/nprocs,dimension,operationpoints[i]);
    }*/


    // double median = (quickselect(distances,numberofpoints,numberofpoints/2)+quickselect(distances,numberofpoints,numberofpoints/2-1))/2;
    // printf("median is %f \n",median);

    // int smaller = smallerthanmedian(median,finalpoints,distances,numberofpoints);
    // printf("\n");

    MPI_Finalize();
   
    return 0 ;
}
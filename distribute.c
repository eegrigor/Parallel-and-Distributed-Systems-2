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

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
//http://www.stat.cmu.edu/~ryantibs/median/quickselect.c
float quickselect(float *arr, int n, int k) {
  unsigned long i,ir,j,l,mid;
  float a,temp;

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
    float *coordinates ;
} point;

void printTable(float **table, int rows,int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			printf(" %f ", table[i][j]);
		}
		printf("\n");
	}

	printf("\n");
}

void printpoints(point *points,int numberofpoints,int dimension){
    
    for(int i = 0 ; i <numberofpoints ; i++){
        printf("%d Point: ",i);
        for(int j = 0 ; j <dimension; j++){
            printf("%f ", points[i].coordinates[j]);
        }
        printf("\n");
    }
    printf("\n");
}

point *makepoints(float **pointsmatrix,int numberofpoints,int dimension){

    point *points = (point *) malloc(numberofpoints * sizeof(point));

    for(int i = 0 ; i < numberofpoints ; i++){
        points[i].coordinates = (float *) malloc(dimension * sizeof(float));
        for(int j = 0 ; j < dimension ; j ++){
            points[i].coordinates[j] = pointsmatrix[i][j] ; 
        } 
    }

    return points;

}

float *distacefrompivot(point pivot,  int numberofpoints, int dimension, point *points){
    float *distances = (float *) malloc(numberofpoints * sizeof(float));

    for(int i = 0 ; i < numberofpoints ; i++){
        float distance = 0 ;
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

int smallerthanmedian(float median,point *points,float *distances,int numberofpoints){
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

    for(int i = 0 ; i < numberofpoints ; i++){
        for(int j = 0 ; j  < numberofoperations ; j++){
            for(int k = 0 ; k < pointsperoperation ; k++){
                partpoints[j][k] = pointsmatrix[i];
            }
        }
    }

    return partpoints ;
}


int main(int argc, char **argv){

    
    int numberofpoints = 8 ;
    int dimension = 4 ;
    int numberofoperations = 2 ;

    float **points = (float **)malloc(numberofpoints * sizeof(float*));
    for(int i = 0 ; i < numberofpoints ; i++){
       points[i] = (float *) malloc(dimension * sizeof(float)); 
    };

    for(int i = 0 ; i < numberofpoints ; i++){
        for(int j = 0 ; j < dimension ; j++){
            points[i][j] = (float)rand()/((float)RAND_MAX) ; 
        }
    }

    printTable(points,numberofpoints,dimension);
    
    
    point *finalpoints = makepoints(points,numberofpoints,dimension);
    point **operationpoints = partionpoints(finalpoints,numberofpoints,numberofoperations);

    printpoints(finalpoints,numberofpoints,dimension);

   

    int pivotnumber = rand() % (numberofpoints+1) ; 

    point pivot = finalpoints[pivotnumber] ; 
    printf("THE PIVOT IS : %d \n" , pivotnumber);
    printf("\n");

    float *distances = distacefrompivot(pivot,numberofpoints,dimension,finalpoints);
    printf("\n");

    // float **partdistances = (float **) malloc(numberofoperations * sizeof(float*));
    // for(int i = 0 ; i < numberofoperations ; i++){
    //     partdistances[i] = (float *) malloc((numberofpoints/numberofoperations) * sizeof(float));
    //     partdistances[i] = distacefrompivot(pivot,numberofpoints/numberofoperations,dimension,operationpoints[i]);
    // }


    float median = (quickselect(distances,numberofpoints,numberofpoints/2)+quickselect(distances,numberofpoints,numberofpoints/2-1))/2;
    printf("median is %f \n",median);

    int smaller = smallerthanmedian(median,finalpoints,distances,numberofpoints);
    printf("\n");

   
    return 0 ;
}
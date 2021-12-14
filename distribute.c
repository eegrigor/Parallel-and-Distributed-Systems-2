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

double **colmajor2rowmajor(double **points,int numberofpoints,int dimension){
    double **retpoints = (double **) malloc(numberofpoints * sizeof(double *));
    for(int i = 0 ; i < numberofpoints ; i++){
        retpoints[i] = (double *) malloc(dimension * sizeof(double));
    } 

    for(int i = 0 ; i < numberofpoints ; i ++){
        for(int j = 0 ; j < dimension ; j++){
            retpoints[i][j] = points[j][i] ;
        }
    }

    return retpoints;
}

// 0 -> train_image
// 1 -> test image 
int main(int argc, char **argv){

    //readbinary("mnist.txt",8,4);
    
    
    int pid,nproc;
    int numberofpoints = 16 ;
    int dimension = 4 ;
    
    
    MPI_Status status;
    

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    double oppoints[numberofpoints/nproc][dimension];
    double pivot[dimension];

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
        
        double **rowpoints = colmajor2rowmajor(points , numberofpoints , dimension);
        printTable(points,dimension,numberofpoints);
        printTable(rowpoints,numberofpoints,dimension);

        for(int i = 0 ; i < numberofpoints/nproc ; i++){
            for(int j = 0 ; j < dimension ;j++){
                oppoints[i][j] = rowpoints[i][j];
            }
        }

        printf("STARTING POINTS OPPERATION %d \n",pid);
        for(int i = 0 ; i < numberofpoints/nproc ; i++){
            for(int j = 0 ; j < dimension ; j++){
                printf("%f ",oppoints[i][j]);
            }
            printf("\n");
        }

        for(int i = 1 ; i < nproc ; i ++){
            for(int j = numberofpoints/nproc * i ; j < (numberofpoints/nproc)*(i+1) ; j++){
                MPI_Send(rowpoints[j], dimension , MPI_DOUBLE,i,55,MPI_COMM_WORLD);
            }
        }

        srand(time(0));
        int pivotnumber = rand() % (numberofpoints/nproc);
        //printf("PIVOT NUMBER IS %d ",pivotnumber);

        for(int i = 0 ; i < dimension ; i++){
            pivot[i] = oppoints[pivotnumber][i];
        //    printf("%f ",pivot[i]);
        }
        
        

    }
    else{
        
        for(int i = 0 ; i < numberofpoints/nproc ; i ++){
            MPI_Recv(oppoints[i],dimension,MPI_DOUBLE,0,55,MPI_COMM_WORLD, &status);
        }
        
        printf("STARTING POINTS OPPERATION %d \n",pid);
        for(int i = 0 ; i < numberofpoints/nproc ; i++){
            for(int j = 0 ; j < dimension ; j++){
                printf("%f ",oppoints[i][j]);
            }
            printf("\n");
        }

       
    }
     
    MPI_Bcast(pivot, dimension , MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double distances[numberofpoints/nproc];
    //if(pid < nproc){
       
        for(int i = 0 ; i < numberofpoints / nproc ; i++){
            double distance = 0;
            for(int j = 0  ; j < dimension ;j++){
                 distance += pow(oppoints[i][j] - pivot[j],2) ;
            }
            distances[i] = distance ;
        }
        
        for(int i = 0 ; i < numberofpoints / nproc ; i++){
            printf("%f from PID : %d ",distances[i],pid);
            
        }
//  }

    printf("\n");
    printf("\n");

    double alldistances[numberofpoints];
    MPI_Gather(distances, numberofpoints/nproc , MPI_DOUBLE , alldistances , numberofpoints/nproc , MPI_DOUBLE , 0 , MPI_COMM_WORLD );

    double median ; 

    if(pid == 0){
        median = (quickselect(alldistances, numberofpoints,numberofpoints/2) + quickselect(alldistances, numberofpoints, numberofpoints/2 - 1))/2;
    }

    MPI_Bcast(&median, 1 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if(pid == 0)
    printf("MEDIAN IS %f \n",median);

    //distributebymedian
    int smaller = 0 ;
    for(int i = 0 ; i < numberofpoints /nproc ;i++){
        if(distances[i] < median)
            smaller++ ;
    }
    double tmp ;
    double distancetmp;
    for(int i = 0 ; i < numberofpoints/nproc - 1 ; i++){
        
        if(distances[i] > median){
            for(int j = i + 1; j < numberofpoints/nproc ; j++){
                if(distances[j] < median){
                    
                    distancetmp = distances[i];
                    distances[i] = distances[j];
                    distances[j] = distancetmp; 
                    //SWAP(distance[i],distances[j]);
                    for(int k = 0 ; k < dimension ; k++){
                        
                        tmp = oppoints[i][k];
                        oppoints[i][k] = oppoints[j][k];
                        oppoints[j][k] = tmp ;
                        //SWAP(oppoints[i][k],oppoints[j][k]);
                        
                    }
                    break ;
                }
                
            }
        }
    }

    if(pid < nproc/2){
        int bigger = numberofpoints/nproc - smaller ;
        smaller = bigger ;
    }
    //  printf("POINTS AFTER SMALLER OPPERATION %d \n",pid);
    // for(int i = 0 ; i < numberofpoints/nproc ; i++){
    //     for(int j = 0 ; j < dimension ; j++){
    //         printf("%f ",oppoints[i][j]);
    //     }
    //     printf("\n");
    // }

    //  for(int i = 0 ; i < numberofpoints / nproc ; i++){
    //         printf("%f from PID : %d ",distances[i],pid);
            
    // }
    // printf("\n");
    //printf("SMALLER  FROM OPPERATION %d : %d \n ",pid,smaller);

    int smalls[nproc];
    MPI_Gather(&smaller, 1 , MPI_INT , smalls , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
    
    int  exchangetable[nproc/2][nproc/2];

    if(pid ==0){
        for(int i = 0 ; i < nproc ; i++){
            printf("%d ",smalls[i]);   
        }
        printf("\n");

        for(int i = 0 ; i < nproc/2 ; i++){
            for(int j = nproc/2 ; j < nproc ; j++){
                exchangetable[i][j-nproc/2] = 0;
                if(smalls[i] > 0 ){
                    if(smalls[j] > 0){
                        if(smalls[i] >= smalls[j]){
                            exchangetable[i][j - nproc/2] = smalls[j];
                            smalls[i] -= smalls[j];
                            smalls[j] =0 ;

                        }
                        else if(smalls[i] < smalls[j]){
                            exchangetable[i][j - nproc/2] = smalls[i] ;
                            smalls[j] -=smalls[i];
                            smalls[i] =0 ;
                        }

                    }
                }
            }
        }

        for(int i = 0 ; i < nproc/2 ; i++){
            for(int j = 0 ; j < nproc/2 ; j++){
                printf("%d ", exchangetable[i][j]);
            }
            printf("\n");
        }
    }

    MPI_Bcast(exchangetable, pow(nproc/2,2) , MPI_INT, 0, MPI_COMM_WORLD);

    if(pid < nproc/2){
        int counter = numberofpoints/nproc - 1 ;
        for(int i = 0 ; i < nproc/2 ; i++){
            if(exchangetable[pid][i] > 0){
                for(int j = 0 ; j < exchangetable[pid][i] ;j++){
                    MPI_Sendrecv(oppoints[counter] , dimension ,MPI_DOUBLE ,i+nproc/2,55,oppoints[counter],dimension,MPI_DOUBLE,i+nproc/2,55,MPI_COMM_WORLD,&status);
                    counter--;
                }
            }
        }
    }
    else{ 
        int counter =0;
        for(int i = 0 ; i < nproc/2 ; i++){
            if(exchangetable[i][pid - nproc/2] > 0){
                for(int j = 0 ; j < exchangetable[i][pid - nproc/2] ;j++){
                    MPI_Sendrecv(oppoints[counter] , dimension ,MPI_DOUBLE , i ,55,oppoints[counter],dimension,MPI_DOUBLE,i,55,MPI_COMM_WORLD,&status);
                    counter++;
                }
            }
        }

    }
   
    printf("FINAL POINTS OPPERATION %d \n",pid);
    for(int i = 0 ; i < numberofpoints/nproc ; i++){
        for(int j = 0 ; j < dimension ; j++){
            printf("%f ",oppoints[i][j]);
        }
        printf("\n");
    }

    // MPI_Bcast(exchangetable, pow(nproc/2 ,2) , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    

    MPI_Finalize();
   
    return 0 ;
}
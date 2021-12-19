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
#include <limits.h>
#include <sys/time.h>

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

//MAKES THE COLUMN-MAJOR MATRIX TO ROW-MAJOR
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

void distributebymedian(int nproc,int numberofpoints,int dimension,int pid,double oppoints[numberofpoints/nproc][dimension],double pivot[]){
     
    MPI_Status status ;
    int lastid = pid ;
    
    //END THE RECURSIVE CALL IF WE HAVE ONLY 1 PROCESS
    if(nproc == 1){
         return ;
    }
    
    // GROUP PROCESSES 
    int color = (pid / nproc) ;

    //MAKE A NEW COMMUNICATOR FOR EACH GROUP 
    MPI_Comm comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, pid, &comm);
    
    //GIVE PROCESSES THE NEW ID IN THE NEW COMMUNICATOR
    MPI_Comm_rank(comm, &pid);

    double distances[numberofpoints/nproc];
    
    //FIND DISTACE FROM EACH POINT    
    for(int i = 0 ; i < numberofpoints / nproc ; i++){
        double distance = 0;
        for(int j = 0  ; j < dimension ;j++){
            distance += pow(oppoints[i][j] - pivot[j],2) ;
        }
        distances[i] = distance ;
    }

    // MASTER PROCESS GATHER ALL DISTANCES   
    double alldistances[numberofpoints];
    MPI_Gather(distances, numberofpoints/nproc , MPI_DOUBLE , alldistances , numberofpoints/nproc , MPI_DOUBLE , 0 , comm );

    double median ; 

    if(pid == 0){
        //FIND MEDIAN WITH QUICKSELECT
        median = (quickselect(alldistances, numberofpoints,numberofpoints/2) + quickselect(alldistances, numberofpoints, numberofpoints/2 - 1))/2;
    }

    //MASTER SENDS MEDIAN TO ALL PROCESSES 
    MPI_Bcast(&median, 1 , MPI_DOUBLE, 0, comm);

    //FIND HOW MANY DISTANCES ARE SMALLER THAN MEDIAN
    int smaller = 0 ;
    for(int i = 0 ; i < numberofpoints /nproc ;i++){
        if(distances[i] < median)
            smaller++ ;
    }
    double tmp ;
    double distancetmp;
    //SPLIT THE "BIG" FROM "SMALL" POINTS 
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
    //IF PROCESS MUST GIVE "BIG" POINTS 
    if(pid < nproc/2){
        int bigger = numberofpoints/nproc - smaller ;
        smaller = bigger ;
    }
    //MASTER GATHER HOW MANY POINTS MUST EACH PROCESS GIVE 
    int smalls[nproc];
    MPI_Gather(&smaller, 1 , MPI_INT , smalls , 1 , MPI_INT , 0 , comm );
    
    int  exchangetable[nproc/2][nproc/2];
    //MASTER MAKE THE EXCHANGE TABLE 
    if(pid ==0){
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
    }
    //SEND TO ALL 
    MPI_Bcast(exchangetable, pow(nproc/2,2) , MPI_INT, 0, comm);

    //TRADE POINTS
    if(pid < nproc/2){
        int counter = numberofpoints/nproc - 1 ;
        for(int i = 0 ; i < nproc/2 ; i++){
            if(exchangetable[pid][i] > 0){
                for(int j = 0 ; j < exchangetable[pid][i] ;j++){
                    MPI_Send(oppoints[counter] , dimension , MPI_DOUBLE,i+nproc/2,55,comm);
                    MPI_Recv(oppoints[counter],dimension,MPI_DOUBLE,i+nproc/2,55,comm, &status);
                    counter--;
                }
            }
        }
    }
    else{ 
        int counter = 0;
        for(int i = 0 ; i < nproc/2 ; i++){
            if(exchangetable[i][pid - nproc/2] > 0){
                for(int j = 0 ; j < exchangetable[i][pid - nproc/2] ;j++){
                    MPI_Send(oppoints[counter] , dimension , MPI_DOUBLE,i,55,comm);
                    MPI_Recv(oppoints[counter],dimension,MPI_DOUBLE,i,55,comm, &status);
                    counter++;
                }
            }
        }

    }
    
    //RECURSIVE CALL
    distributebymedian(nproc/2,numberofpoints/2,dimension,lastid,oppoints,pivot);
    
    return ;
}


int main(int argc, char **argv){

    
    
    
    int pid,nproc;
    int numberofpoints =pow(2,14);
    int dimension = 256 ;
    
    
    MPI_Status status;
    
    srand(time(0));
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    //POINTS FOR EACH PROCESS
    double oppoints[numberofpoints/nproc][dimension];
    double pivot[dimension];
    
    struct timeval stop , start ;

    if(pid == 0 ){
        
        double **points = (double **)malloc(dimension * sizeof(double*));
        for(int i = 0 ; i < dimension ; i++){
            points[i] = (double *) malloc(numberofpoints * sizeof(double)); 
        }
        //MAKE THE POINTS(I KNOW IT IS THE WRONG WAY BECAUSE WE USE MPI FOR MEMORY AND HERE I MAKE POINTS IN ONE PROCESS)
        for(int i = 0 ; i < dimension ; i++){
            for(int j = 0 ; j < numberofpoints ; j++){
                points[i][j] = (double)rand()/((double)RAND_MAX) ; 
        
            }
        }
        //MAKE POINTS ROW-MAJOR
        double **rowpoints = colmajor2rowmajor(points , numberofpoints , dimension);
    
        for(int i = 0 ; i < numberofpoints/nproc ; i++){
            for(int j = 0 ; j < dimension ;j++){
                oppoints[i][j] = rowpoints[i][j];
            }
        }
        gettimeofday(&start, NULL);
        //MASTER SEND POINTS TO EACH PROCESS
        for(int i = 1 ; i < nproc ; i ++){
            for(int j = numberofpoints/nproc * i ; j < (numberofpoints/nproc)*(i+1) ; j++){
                MPI_Send(rowpoints[j], dimension , MPI_DOUBLE,i,55,MPI_COMM_WORLD);
            }
        }

        int pivotnumber = rand() % (numberofpoints/nproc);
       
        //FIND PIVOT FROM MASTER'S POINT 
        for(int i = 0 ; i < dimension ; i++){
            pivot[i] = oppoints[pivotnumber][i];
        }
    }
    else{
        for(int i = 0 ; i < numberofpoints/nproc ; i ++){
            MPI_Recv(oppoints[i],dimension,MPI_DOUBLE,0,55,MPI_COMM_WORLD, &status);
        }
    }

    //SEND PIVOT TO ALL PROCESSES
    MPI_Bcast(pivot, dimension , MPI_DOUBLE, 0, MPI_COMM_WORLD);

    distributebymedian( nproc,numberofpoints,dimension,pid,oppoints,pivot);

    double finaldistances[numberofpoints/nproc];
     
    for(int i = 0 ; i < numberofpoints / nproc ; i++){
        double distance = 0;
        for(int j = 0  ; j < dimension ;j++){
                distance += pow(oppoints[i][j] - pivot[j],2) ;
        }
        finaldistances[i] = distance ;
    }

    double max, min;
    min = max = finaldistances[0];

    //EACH PROCESS FINDS MAX AND MIN DISTANCE
    for(int i = 0; i < numberofpoints/nproc ; i++){
        if(finaldistances[i] < min){
            min = finaldistances[i];
        }
        if(finaldistances[i] > max){
            max = finaldistances[i];
        }
    }       

        double maxarray[nproc];
        double minarray[nproc];

    MPI_Gather(&max, 1 ,MPI_DOUBLE , maxarray , 1 ,MPI_DOUBLE , 0 , MPI_COMM_WORLD );
    MPI_Gather(&min, 1 ,MPI_DOUBLE , minarray , 1 , MPI_DOUBLE , 0 , MPI_COMM_WORLD );

    if(pid ==0){
    
        int i ;
        //MASTER CHECKS IF  ORDER IS RIGHT
        for(i =1; i< nproc; i++){
            if(maxarray[i-1] <= minarray[i]){
                continue;
            }
            else
                break;
        }
        printf("MINS \n");
        for(int i = 0 ;i < nproc ;i++ ){
            printf("%f ",minarray[i]);
        }
        printf("\n");

        printf("MAXS \n");

        for(int i = 0 ;i < nproc ;i++ ){
            printf("%f ",maxarray[i]);
        }
        printf("\n");

        
        if(i == nproc){
            printf("SUCCESS");
        }
        else {
            printf("FAIL");
        }
        gettimeofday(&stop, NULL);
        
        uint timediff = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
        
        printf("TIME TAKEN :%u ",timediff);
    }
    
    MPI_Finalize();
   
    return 0 ;
}
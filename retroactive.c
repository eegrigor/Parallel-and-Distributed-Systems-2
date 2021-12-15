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

void distributebymedian(int nproc,int numberofpoints,int dimension,int pid,double oppoints[numberofpoints/nproc][dimension],double pivot[]){
     
    MPI_Status status ;
    int lastid = pid ;
     if(nproc == 1){
         return ;
     }
     int color = (pid / nproc) ;

     MPI_Comm comm;
     MPI_Comm_split(MPI_COMM_WORLD, color, pid, &comm);
    
     MPI_Comm_rank(comm, &pid);
     //MPI_Comm_size(comm, &nproc);


    double distances[numberofpoints/nproc];
    
       
        for(int i = 0 ; i < numberofpoints / nproc ; i++){
            double distance = 0;
            for(int j = 0  ; j < dimension ;j++){
                 distance += pow(oppoints[i][j] - pivot[j],2) ;
            }
            distances[i] = distance ;
        }
        
    //     for(int i = 0 ; i < numberofpoints / nproc ; i++){
    //         printf("%f from PID : %d ",distances[i],pid);
            
    //     }


    // printf("\n");
    // printf("\n");

    double alldistances[numberofpoints];
    MPI_Gather(distances, numberofpoints/nproc , MPI_DOUBLE , alldistances , numberofpoints/nproc , MPI_DOUBLE , 0 , comm );

    double median ; 

    if(pid == 0){
        median = (quickselect(alldistances, numberofpoints,numberofpoints/2) + quickselect(alldistances, numberofpoints, numberofpoints/2 - 1))/2;
    }

    MPI_Bcast(&median, 1 , MPI_DOUBLE, 0, comm);
    // if(pid == 0)
    // printf("MEDIAN IS %f \n",median);

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
    MPI_Gather(&smaller, 1 , MPI_INT , smalls , 1 , MPI_INT , 0 , comm );
    
    int  exchangetable[nproc/2][nproc/2];

    if(pid ==0){
    //    printf("NUMBER OF PRECESS %d \n" , nproc);
    //     for(int i = 0 ; i < nproc ; i++){
            
    //         printf("%d ",smalls[i]);   
    //     }
    //     printf("\n");

        

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

        // for(int i = 0 ; i < nproc/2 ; i++){
        //     for(int j = 0 ; j < nproc/2 ; j++){
        //         printf("%d ", exchangetable[i][j]);
        //     }
        //     printf("\n");
        // }
    }

    MPI_Bcast(exchangetable, pow(nproc/2,2) , MPI_INT, 0, comm);

    if(pid < nproc/2){
        int counter = numberofpoints/nproc - 1 ;
        for(int i = 0 ; i < nproc/2 ; i++){
            if(exchangetable[pid][i] > 0){
                for(int j = 0 ; j < exchangetable[pid][i] ;j++){
                    //MPI_Sendrecv(oppoints[counter] , dimension ,MPI_DOUBLE ,i+nproc/2,55,oppoints[counter],dimension,MPI_DOUBLE,i+nproc/2,55,comm,&status);
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
                    //MPI_Sendrecv(oppoints[counter] , dimension ,MPI_DOUBLE , i ,55,oppoints[counter],dimension,MPI_DOUBLE,i,55,comm,&status);
                    MPI_Send(oppoints[counter] , dimension , MPI_DOUBLE,i,55,comm);
                    MPI_Recv(oppoints[counter],dimension,MPI_DOUBLE,i,55,comm, &status);
                    counter++;
                }
            }
        }

    }
   
    

    
    distributebymedian(nproc/2,numberofpoints/2,dimension,lastid,oppoints,pivot);
    
    return ;
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

        // printf("STARTING POINTS OPPERATION %d \n",pid);
        // for(int i = 0 ; i < numberofpoints/nproc ; i++){
        //     for(int j = 0 ; j < dimension ; j++){
        //         printf("%f ",oppoints[i][j]);
        //     }
        //     printf("\n");
        // }

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
        
        // printf("STARTING POINTS OPPERATION %d \n",pid);
        // for(int i = 0 ; i < numberofpoints/nproc ; i++){
        //     for(int j = 0 ; j < dimension ; j++){
        //         printf("%f ",oppoints[i][j]);
        //     }
        //     printf("\n");
        // }

       
    }
    MPI_Bcast(pivot, dimension , MPI_DOUBLE, 0, MPI_COMM_WORLD);

    distributebymedian( nproc,numberofpoints,dimension,pid,oppoints,pivot);

    printf("FINAL POINTS OPPERATION %d \n",pid);
    for(int i = 0 ; i < numberofpoints/nproc ; i++){
        for(int j = 0 ; j < dimension ; j++){
            printf("%f ",oppoints[i][j]);
        }
        printf("\n");
    }


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
        for(i = 0 ; i< nproc-1 ; i++){
            if(maxarray[i] <= minarray[i+1]){
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


        if(i == nproc - 1){
            printf("SUCCESS");
        }
        else {
            printf("FAIL");
        }
    }


    
    

    MPI_Finalize();
   
    return 0 ;
}
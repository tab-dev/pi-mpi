/**
 *Author: Timothy Bennett
  Date: 10/03/2018
  Course: CIS 5990 
  Assignment 2 - Approximation of PI using Trapezoid Method with Reduction, Timing, and File I/O
  Description: The program is based off of examples provided by the 
              "Intro to Parallel Programming" textbook by Peter S. Pacheco 
              and from the course material/example code. The program approximates 
              the value of PI using the Trapezoid Method and a series function. 

              10/03/2018 - Modified to uses file I/O instead of command-line input. Updated to use Reduction operation for all calculated areas. Added timing of entire PI calculation for Process 0. Modified algorithm to work with any ratio of processes to intervals.  
               
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define INPUT_FILE "input.txt"
#define ERROR_CODE -1
#define MICRO_PER_SEC 1000000L

double calcTotalTrapArea(int n,double h,int rank,int procSize);
double calcTrapArea(double base1,double base2,double height);
double piFunction(double x);
int readIntervalsFromFile(char *filename);
void createIntevalMPIType(MPI_Datatype *type);
long getElapsedTimeMicro(struct timeval *start,struct timeval *end);

/**
 *The program first retrieves the command line argument for N sub-intervals. By default,
  the program calculates the area over the interval [0,1]. Each height h is equalivalent so 
  its value is calcuated once which is inversely related to N. Once all local varibles and MPI 
  have been initialized, each process calculates N/P sub-intervals, where P is the number of processes.
  Once each process has finished their calculations, they will send their results to process 0, which will
  add all sub-intervals to its sub-interval for a final approximation of PI. 
 *
 */
int main(int argc,char *args[]){
    int rank,size,n=0;
    double a=0.0,b=1.0,h=0.0;
    MPI_Datatype intervalType;
    MPI_Init(&argc,&args);

    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    if(rank == 0){ 
        n = readIntervalsFromFile(INPUT_FILE);
        if(n == ERROR_CODE){
            MPI_Finalize();
	    return ERROR_CODE;
	}
    }
    createIntevalMPIType(&intervalType);

    MPI_Bcast(&n,1,intervalType,0,MPI_COMM_WORLD);

    MPI_Type_free(&intervalType);
     
    struct timeval start,end;
    gettimeofday(&start,NULL);
   
    h=(b-a)/n;
   
    double my_area = calcTotalTrapArea(n,h,rank,size);
    double pi = 0.0;
    pi = my_area;
    MPI_Reduce(&my_area,&pi,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    gettimeofday(&end,NULL);
    
    char result[200];
    sprintf(result,"Area sum for process %d: %7.14f\n",rank,my_area);
    if(rank != 0){
    	MPI_Send(result,strlen(result)+1,MPI_CHAR,0,0,MPI_COMM_WORLD);
    }else{
    	int i = 0;
        char message[200];
	fprintf(stdout,result);
    	for(i = 1;i<size;i++){
	    MPI_Recv(message,200,MPI_CHAR,i,0,MPI_COMM_WORLD,NULL);
	    fprintf(stdout,message);
	}
    }
    if(rank == 0){
	fprintf(stdout,"Estimation of PI is %7.14f\n",pi);
    	fprintf(stdout,"Time to complete calculation is %ld microseconds.\n",getElapsedTimeMicro(&start,&end));
    }

    
    MPI_Finalize();

    return 0;
}

long getElapsedTimeMicro(struct timeval *start,struct timeval *end){
    long elaspedTime = (end->tv_sec*MICRO_PER_SEC+ end->tv_usec) - (start->tv_sec*MICRO_PER_SEC + start->tv_usec);
    return elaspedTime;
} 

/*
 *The function calculates the area of one or more trapezoid areas based on the value of parameter n.
 *The function use cyclic decomposition to partition the intervals based on rank. The function returns the local 
  sum of all of its calculated areas assigned to the process. The  function is based on a combination of 
  the trazpezoid method for approxiating an area under a function  and the example function provided the 
  textbook "An Introduction to Parallel Programming" pg. 99. The f(x) fucntion is based off of a series 
  function for approimating pi.
 */
double calcTotalTrapArea(int n,double h,int rank,int procSize){
    int i = 0;
    double area = 0.0;

    for(i=rank;i<n;i+=procSize){
       area +=calcTrapArea(piFunction(i*h),piFunction(i*h + h),h);
    }

    return area;
}

double calcTrapArea(double base1,double base2,double height){
    double area = 0.0;

    area += (height/2)*(base1+base2);

    return area;
}

/**
 * This function is in part of a series to calulate pi. The function simply takes in a value x and maps to f(x) by returning the 
 * resulting portion of the series. 
 */
double piFunction(double x){
  return 4.0/(1.0+x*x);
}

/**
 * This function reads the number of intervals from a local file. The function takes in a string as the filename. The
 * function will read in the number of intervals from this file. The file should simply have the number of intervals as
 * an integer. Any other data in the file will cause the program to abnormally terminate. If the file does not exist, 
 * the program will alert the user and terminate.
 */
int readIntervalsFromFile(char *filename){
    int n = 0;
    size_t readSize;
    size_t dataSize = 0;
    char *buffer;
    FILE *file = fopen(filename,"r");

    
    if(file == NULL){
	fprintf(stderr,"Unable to open file %s\n",filename);
    	return ERROR_CODE;
    }

    readSize = getline(&buffer,&dataSize,file);

    fclose(file);

    if(readSize == -1){
    	fprintf(stderr,"Error while reading file %s. File is empty!\n",filename);
    	return ERROR_CODE;
    }

    n = atoi(buffer);
    
    if(n == 0){
	fprintf(stderr,"Error, invalid interval data contained in file %s.\n",filename);
    	return ERROR_CODE;
    }

  
    fprintf(stdout,"Interval loaded is %d\n",n); 
    return n;
}

/**
 * This function simplifies the process of creating an MPI type for the interval. The function creates a MPI type
 * for the interval N.
 */
void createIntevalMPIType(MPI_Datatype *type){
    MPI_Aint displacement[1] = {0};
    int lengths[1] = {1};
    MPI_Datatype types[1] = {MPI_INT};

    MPI_Type_create_struct(1,lengths,displacement,types,type);

    MPI_Type_commit(type);
}

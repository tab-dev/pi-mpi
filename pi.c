/**
 *Author: Timothy Bennett
  Date: 09/28/2018
  Course: CIS 5990 
  Assignment 1 - Approximation of PI using Trapezoid Method
  Description: The program is based off of examples provided by the 
              "Intro to Parallel Programming" textbook by Peter S. Pacheco 
              and from the course material/example code. The program approximates 
              the value of PI using the Trapezoid Method and a series function. 
               
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

double calcTrapArea(double x,double x1,int n,double base);
double piFunction(double x);

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
    int my_n = 0;
    double a=0.0,b=1.0,h=0.0;
    double my_a = 0.0,my_b = 0.0;

    MPI_Init(&argc,&args);

    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
   
    if(argc < 2){
    	if(rank == 0){
	    fprintf(stderr,"Invalid arguments. Usage: pi <num_intervals>\n");
	}
        MPI_Finalize();
        return -1;
    } 

    n=atoi(args[1]);

    h=(b-a)/n;
    my_n = n/size;
   
    my_a = a + (double)rank*(double)my_n*h;
    my_b = my_a + my_n*h;

    double my_area = calcTrapArea(my_a,my_b,my_n,h);
    double pi = 0.0,tmp=0.0;
    if(rank != 0){
        MPI_Send(&my_area,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    }else{
	MPI_Status status;
	printf("Process 0 calculated an area of %7.14f\n",my_area);
        pi = my_area;
        int i=1;
        for(i=1;i<size;i++){
	    MPI_Recv(&tmp,1,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	    pi+=tmp;
	    printf("Process %d calculated an area of %7.14f\n",status.MPI_SOURCE,tmp);
	}
	printf("Result pi approximation: %7.14f\n",pi);
    }

    MPI_Finalize();

    return 0;
}

/*
 *The function calculates the area of one or more trapezoid areas based on the value of parameter n.
 *The function iterates from 1 to n-1 (inclusive) and returns the sum of all of its calculated areas. The function 
  is based on a combination of the trazpezoid method for approxiating an area under a function and the example function
  provided the textbook "An Introduction to Parallel Programming" pg. 99. The f(x) fucntion is based off of a series 
  function for approimating pi.
 */
double calcTrapArea(double x,double x1,int n,double base){
    int i = 0;
    double area = 0.0;

    area = (piFunction(x) + piFunction(x1))/2.0;
    for(i=1;i<=n-1;i++){
       area +=piFunction(x+(double)(i*base));
    }
    area *=base;

    return area;
}

/**
 * This function is in part of a series to calulate pi. The function simply takes in a value x and maps to f(x) by returning the 
 * resulting portion of the series. 
 */
double piFunction(double x){
  return 4.0/(1.0+x*x);
}

#include <stdio.h>
#include <math.h>


#define N 2187
#define REPS 200000

void sum(void); 
double second(void); 

float a[N], b[N], c[N], d=5.0;

int main(int argv, char **argc){

  int i,j;
  double time1, time2; 

  
  for (j=0;j<N;j++) {
    a[j] = 7.0;
    b[j] = 6.0;
    c[j] = 23.0;
  }  

  time1=second();
  for (i=0;i<REPS;i++){
     sum(); 
  }
  if (a[0] < 0) printf("Help!\n"); 
  time2=second();

  printf("a[N-1] = %f\n",a[N-1]); 
  printf("Time = %f\n",(time2-time1));   
  printf("Mflops = %f\n",((float)(2*N*REPS))/(1000000.0*(time2-time1)));
  
	   
}

void sum(){

   int j; 

   for (j=0;j<N;j++) {
     a[j] = b[j] + d * c[j];
   }

}


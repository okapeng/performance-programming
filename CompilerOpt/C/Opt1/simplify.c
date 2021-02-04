#include <stdio.h>
#include <math.h>


#define N 2187
#define REPS 50000

void fred(int*); 
double second(void); 


int main(int argv, char **argc){

  int i,a;
  double time1, time2; 



  a=7; 

  time1=second();
  for (i=0;i<REPS;i++){
     fred(&a); 
  }
  if (a < 0) printf("Help!\n"); 
  time2=second();

  printf("a = %d\n",a); 
  printf("Time = %f\n",(time2-time1));   
  printf("Iterations per second %f\n",((float)(N*REPS))/(time2-time1));
  
	   
}

void fred(int* pa){

   int i,j; 
   int a,b,c,d,e;

   a = *pa;
   i=0;
   e=4;
   c=23;

   for (j=0;j<=N;j++) {
     e = c - j - j;
     b = c - j + c - j + c - j;
     d = b + e;
     i += 2;
     d = i;
     a += d + d + d;
     
   }

   *pa = a; 
}


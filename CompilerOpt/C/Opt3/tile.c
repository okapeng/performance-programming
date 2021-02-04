#include <stdio.h>
#include <math.h>


#define N 729
#define REPS 1

void matmul(void); 
double second(void); 

float a[N][N], b[N][N], c[N][N];

int main(int argv, char **argc){

  int i,j;
  double time1, time2; 

  for (i=0;i<N;i++) { 
    for (j=0;j<N;j++) {
      a[i][j] = 7.0;
      b[i][j] = 6.0;
      c[i][j] = 23.0;
    }  
  }

  time1=second();
  for (i=0;i<REPS;i++){
     matmul(); 
  }
  if (c[0][0] < 0) printf("Help!\n"); 
  time2=second();

  printf("c[0][0] = %f\n",c[0][0]); 
  printf("Time = %f\n",(time2-time1));   
  printf("Mflops = %f\n",((float)(2*N*N) * (float)(N*REPS))/(1000000.0*(time2-time1)));
  
	   
}

void matmul(){

   int i,j,k; 

   for (k=0;k<N;k++) {
     for (i=0;i<N;i++) {
       for (j=0;j<N;j++) {
	 c[j][i] += a[k][i] * b[j][k];
       }
     }
   }

}


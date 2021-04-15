#include <math.h>

void vis_forces(int N,double *f, double* restrict vis, double* restrict velo)
{
  int i;
          for(i=0;i<N;i++){
            f[i] = -vis[i] * velo[i];
          }
}
void wind_force(int N,double *f, double* restrict vis, double velo)
{
  int i;
          for(i=0;i<N;i++){
            f[i] = f[i] -vis[i] * velo;
          }
}
void add_norms(int N,double *r, double* restrict delta)
{
  int k;
        for(k=0;k<N;k++){
          r[k] += (delta[k] * delta[k]);
        }
}

__inline double force(double W, double delta, double r){
  return W*delta/(pow(r,3.0));
}




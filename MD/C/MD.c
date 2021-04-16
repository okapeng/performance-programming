/*
 *  Simple molecular dynamics code.
 *  2021
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "coord.h"
#include <omp.h>
// #pragma omp parallel num_threads(8)

void vis_forces(int N,double *f, double *vis, double *vel);
void add_norms(int N,double *r, double *delta);
double force(double W, double delta, double r);
void wind_force(int N,double *f, double *vis, double vel);


void evolve(int count,double dt){
int step;
int i,j,k,l;
int have_collided;
double Size;
/*
 * Loop over timesteps.
 */
      for(step = 1;step<=count;step++){
        printf("timestep %d\n",step);
        printf("collisions %d\n",collisions);

/* set the viscosity term in the force calculation */
        #pragma ivdep
        for(j=0;j<Ndim;j++){
          vis_forces(Nbody,f[j],vis,velo[j]);
        }
/* add the wind term in the force calculation */
        #pragma ivdep
        for(j=0;j<Ndim;j++){
          wind_force(Nbody,f[j],vis,wind[j]);
        }
/* calculate distance from central mass */
        __assume_aligned(r, 16);
        for(k=0;k<Nbody;k++){
          r[k] = 0.0;
        }
        #pragma ivdep
        for(i=0;i<Ndim;i++){
	  add_norms(Nbody,r,pos[i]);
        }
        // #pragma ivdep
        #pragma omp simd aligned(r:16)
        for(k=0;k<Nbody;k++){
          r[k] = sqrt(r[k]);
        }
       /* calculate central force */
        for(l=0;l<Ndim;l++){
	        for(i=0;i<Nbody;i++){
                f[l][i] = f[l][i] - 
                   force(G*mass[i]*M_central,pos[l][i],r[i]);
	  }
	}
/* calculate pairwise separation of particles */
        for(l=0;l<Ndim;l++){
          for(i=0;i<Nbody;i++){
            k = 0;
            #pragma vector aligned
            #pragma omp simd
            for(j=i+1;j<Nbody;j++){
              delta_pos[l][k] = pos[l][i] - pos[l][j];
              k = k + 1;
            }
          }
        }

/* calculate norm of separation vector */
        __assume_aligned(delta_r, 16);
        for(k=0;k<Npair;k++){
          delta_r[k] = 0.0;
        }
        for(i=0;i<Ndim;i++){
	  add_norms(Npair,delta_r,delta_pos[i]);
        }
        // #pragma omp parallel for schedule(static, 4)
        #pragma ivdep
        for(k=0;k<Npair;k++){
          delta_r[k] = sqrt(delta_r[k]);
        }

/*
 * add pairwise forces.
 */
        k = 0;
        int* have_collided[Nbody];
        have_collided[0] = _mm_malloc(Nbody*Nbody*sizeof(int),64);
        #pragma vector aligned
        for(l=0;l<Ndim;l++){
          k=0;
          for(i=0;i<Nbody;i++){
            for(j=i+1;j<Nbody;j++){
              if(l==0) have_collided[i][j]=0;
              Size = radius[i] + radius[j];
/*  flip force if close in */
              if( delta_r[k] >= Size ){
                f[l][i] = f[l][i] - 
                   force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
                f[l][j] = f[l][j] + 
                   force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
              }else{
                f[l][i] = f[l][i] + 
                   force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
                f[l][j] = f[l][j] - 
                   force(G*mass[i]*mass[j],delta_pos[l][k],delta_r[k]);
		have_collided[i][j]+=1;
              }
            k = k + 1;
	    if( have_collided[i][j] == 1 ){
	      collisions++;
	    }
            }
          }
        }

/* update positions */
#pragma ivdep
        for(j=0;j<Ndim;j++){
          #pragma vector aligned
          #pragma omp simd
          for(i=0;i<Nbody;i++){
            pos[j][i] = pos[j][i] + dt * velo[j][i];
          }
        }

/* update velocities */
        for(j=0;j<Ndim;j++){
          #pragma vector aligned
          #pragma ivdep
          #pragma omp simd aligned(velo:64)
           for(i=0;i<Nbody;i++){
            velo[j][i] = velo[j][i] + dt * (f[j][i]/mass[i]);
          }
        }


      }

}





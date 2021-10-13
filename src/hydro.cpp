#include "sd3d.hpp"

double sound_speed(double rho, double p){
    double c2 = gmma*p/rho;
    c2 = max(c2,min_c2);
    return sqrt(c2);
}

void prim_to_cons(double *W, double *U, int N){
    double kin;
    for(int i=0; i<N; i++){
        U[i] = W[i];
        kin=0.0;
        #ifdef X
        U[i+_vx_*N] = W[i]*W[i+_vx_*N];
        kin += U[i+_vx_*N]*W[i+_vx_*N];
        #endif
        #ifdef Y
        U[i+_vy_*N] = W[i]*W[i+_vy_*N];
        kin += U[i+_vy_*N]*W[i+_vy_*N];
        #endif
        #ifdef Z
        U[i+_vz_*N] = W[i]*W[i+_vz_*N];
        kin += U[i+_vz_*N]*W[i+_vz_*N];
        #endif
        U[i+_e_*N] = W[i+_p_*N]/(gmma-1.0)+0.5*kin;
    }   
}

void cons_to_prim(double *U, double *W, int N){
    double kin,dim;
    for(int i=0; i<N; i++){
        W[i] = U[i];
        kin=0.0;
        #ifdef X
        W[i+_vx_*N] = U[i+_vx_*N]/W[i];
        kin += U[i+_vx_*N]*W[i+_vx_*N];
        #endif
        #ifdef Y
        W[i+_vy_*N] = U[i+_vy_*N]/W[i];
        kin += U[i+_vy_*N]*W[i+_vy_*N];
        #endif
        #ifdef Z
        W[i+_vz_*N] = U[i+_vz_*N]/W[i];
        kin += U[i+_vz_*N]*W[i+_vz_*N];
        #endif
        W[i+_p_*N] = (gmma-1.0)*(U[i+_e_*N]-0.5*kin);
    }   
}

void fluxes(double *U, double *W, double *F, int N, int _v1_, int _v2_, int _v3_){
    double vel;
    for(int i=0; i<N; i++){
        vel = W[i+_v1_*N];
        F[i]        = U[i       ]*vel;
        F[i+_v1_*N] = U[i+_v1_*N]*vel + W[i+_p_*N];
        #ifdef _2D_
        F[i+_v2_*N] = U[i+_v2_*N]*vel;
        #endif
        #ifdef _3D_
        F[i+_v3_*N] = U[i+_v3_*N]*vel;
        #endif
        F[i+ _p_*N] = (U[i+_e_*N]+W[i+_p_*N])*vel;
    }   
}

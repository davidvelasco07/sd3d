#include "sd3d.hpp"

void compute_dt(){
    double c_s,c_max,dt_max;
    double dt_min=1E3;//Just a big number
    double dx_min=1E3;
    double *u = malloc_host<double>(nvar);
    double *w = malloc_host<double>(nvar);
    //transform_sp_to_cv(U_sp,U_cv);    
    for(int i=0;i<size_cv;i++){
        for(int var=0;var<nvar;var++){
            u[var] = U_sp[i+var*size_cv];
        }
        cons_to_prim(u,w,1);
        c_s = sound_speed(w[0],w[_p_]);
        c_max = 0.0;
        #ifdef X
        c_max += (abs(w[_vx_]) + c_s);
        dx_min = min(dx_min,dx);
        #endif
        #ifdef Y
        c_max += (abs(w[_vy_]) + c_s);
        dx_min = min(dx_min,dy);
        #endif
        #ifdef Z
        c_max += (abs(w[_vz_]) + c_s);
        dx_min = min(dx_min,dz);
        #endif
        dt_max=cfl*dx_min/c_max/n_cv;
        dt_min = min(dt_min,dt_max);
    }
    #ifdef MPI
    MPI_Allreduce(&dt_min, &dt, 1, MPI_DOUBLE, MPI_MIN, Comm);
    #else
    dt=dt_min;
    #endif
}

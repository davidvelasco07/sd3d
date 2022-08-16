#include "sd3d.hpp"

void riemann_llf(double *F, double *U_L, double *U_R, int faces, int _v1_, int _v2_, int _v3_){
    double *u_L = malloc_host<double>(nvar);
    double *u_R = malloc_host<double>(nvar);
    double *w_L = malloc_host<double>(nvar);
    double *w_R = malloc_host<double>(nvar);
    double *f_L = malloc_host<double>(nvar);
    double *f_R = malloc_host<double>(nvar);
    double c_L,c_R,c_max;
    for(int i=0; i<faces; i++){
        for(int var=0; var<nvar; var++){
            u_L[var] =  U_L[i+var*faces];
            u_R[var] =  U_R[i+var*faces];
        }
        cons_to_prim(u_L, w_L, 1);
        cons_to_prim(u_R, w_R, 1);  
        fluxes(u_L,w_L,f_L,1,_v1_,_v2_,_v3_);
        fluxes(u_R,w_R,f_R,1,_v1_,_v2_,_v3_);
        c_L = sound_speed(w_L[0],w_L[_p_])+abs(w_L[_v1_]);
        c_R = sound_speed(w_R[0],w_R[_p_])+abs(w_R[_v1_]);
        c_max = max(c_L,c_R);
        for(int var=0; var<nvar; var++){
            F[i+var*faces] = 0.5*(f_R[var]+f_L[var])-0.5*c_max*(u_R[var]-u_L[var]);
        }
    }
}

void riemann_hllc(double *F, double *U_L, double *U_R, int faces, int _v1_, int _v2_, int _v3_, bool conservatives){
    double *u_L = malloc_host<double>(nvar);
    double *u_R = malloc_host<double>(nvar);
    double *w_L = malloc_host<double>(nvar);
    double *w_R = malloc_host<double>(nvar);
    double c_L,c_R,c_max;
    double v_L,v_R,s_L,s_R,rc_L,rc_R;
    double v_star,p_star;
    double r_starL,r_starR,e_starL,e_starR;
    double r_gdv,v_gdv,p_gdv,e_gdv;
    for(int i=0; i<faces; i++){
        for(int var=0; var<nvar; var++){
            if(conservatives){
                u_L[var] =  U_L[i+var*faces];
                u_R[var] =  U_R[i+var*faces];
            }
            else{
                w_L[var] =  U_L[i+var*faces];
                w_R[var] =  U_R[i+var*faces];
            }
        }
        if(conservatives){
            cons_to_prim(u_L, w_L, 1);
            cons_to_prim(u_R, w_R, 1);
        }
        else{
            prim_to_cons(w_L, u_L, 1);
            prim_to_cons(w_R, u_R, 1);
        }

        c_L = sound_speed(w_L[0],w_L[_p_])+abs(w_L[_v1_]);
        c_R = sound_speed(w_R[0],w_R[_p_])+abs(w_R[_v1_]);
        c_max = max(c_L,c_R);

        v_L = w_L[_v1_];
        v_R = w_R[_v1_];

        //Compute HLL wave speed
        s_L = min(v_L,v_R)-c_max;
        s_R = max(v_L,v_R)+c_max;
        //Compute lagrangian sound speed
        rc_L = w_L[0]*(v_L-s_L);
        rc_R = w_R[0]*(s_R-v_R);
        //Compute acoustic star state
        v_star = (rc_R*v_R + rc_L*v_L + (w_L[_p_]-w_R[_p_]))/(rc_R+rc_L);
        p_star = (rc_R*w_L[_p_] + rc_L*w_R[_p_] + rc_L*rc_R*(v_L-v_R))/(rc_R+rc_L);
        //Left star region variables
        r_starL = w_L[0]*(s_L-v_L)/(s_L-v_star);
        e_starL = ((s_L-v_L)*u_L[_e_]-w_L[_p_]*v_L+p_star*v_star)/(s_L-v_star);
        //Right star region variables
        r_starR = w_R[0]*(s_R-v_R)/(s_R-v_star);
        e_starR = ((s_R-v_R)*u_R[_e_]-w_R[_p_]*v_R+p_star*v_star)/(s_R-v_star);

        //If   s_L>0 -> U_gdv = U_L
        //elif v*>0  -> U_gdv = U*_L
        //eilf s_R>0 -> U_gdv = U*_R
        //else       -> U_gnv = U_R
        if(s_L>0){
            r_gdv = w_L[0];
            v_gdv = w_L[_v1_];
            p_gdv = w_L[_p_];
            e_gdv = u_L[_e_];
        }
        else if(v_star>0){
            r_gdv = r_starL;
            v_gdv = v_star;
            p_gdv = p_star;
            e_gdv = e_starL;
        }
        else if(s_R>0){
            r_gdv = r_starR;
            v_gdv = v_star;
            p_gdv = p_star;
            e_gdv = e_starR;
        }
        else{
            r_gdv = w_R[0];
            v_gdv = w_R[_v1_];
            p_gdv = w_R[_p_];
            e_gdv = u_R[_e_];
        }
 
        F[i           ] = r_gdv*v_gdv;
        F[i+_v1_*faces] = r_gdv*v_gdv*v_gdv + p_gdv;
        #ifdef _2D_
        F[i+_v2_*faces] = r_gdv*v_gdv*(v_star>0 ? w_L[_v2_] : w_R[_v2_]); 
        #endif
        #ifdef _3D_
        F[i+_v3_*faces] = r_gdv*v_gdv*(v_star>0 ? w_L[_v3_] : w_R[_v3_]);
        #endif
        F[i+_p_*faces]  = v_gdv*(e_gdv + p_gdv);
    }
}

void riemann_solver_x(){
    //Store fluxes at the interface between elements to then solve the unique flux
    int cell,face;
    int fp_x=n_fp*cells_x;
    for(int var=0;var<nvar;var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<cells_x-1;i++){
                        cell = (i*n_fp) + j*fp_x         + k*fp_x*cv_y         + i_ader*fp_x*cv_y*cv_z         + var*fp_x*cv_y*cv_z*n_cv;
                        face = i        + j*(cells_x-1)  + k*(cells_x-1)*cv_y  + i_ader*(cells_x-1)*cv_y*cv_z  + var*(cells_x-1)*cv_y*cv_z*n_cv;
                        U_R[face] = U_ader_fp_x[cell+n_fp];
                        U_L[face] = U_ader_fp_x[cell+n_cv];
                    }
                }
            }
        }
    }
    //This should be a pointer towards whatever riemann solver was defined
    riemann_llf(F_x,U_L,U_R,(cells_x-1)*cv_y*cv_z*n_cv,_vx_,_vy_,_vz_);

    for(int var=0;var<nvar;var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<cells_x-1;i++){
                        cell = (i*n_fp) + j*fp_x         + k*fp_x*cv_y         + i_ader*fp_x*cv_y*cv_z         + var*fp_x*cv_y*cv_z*n_cv;
                        face = i        + j*(cells_x-1)  + k*(cells_x-1)*cv_y  + i_ader*(cells_x-1)*cv_y*cv_z  + var*(cells_x-1)*cv_y*cv_z*n_cv;
                        F_ader_fp_x[cell+n_fp] = F_x[face];
                        F_ader_fp_x[cell+n_cv] = F_x[face];
                    }
                }
            }
        }
    }
}

void riemann_solver_y(){
    //Store fluxes at the interface between elements to then solve the unique flux
    int cell,face;
    int fp_y=n_fp*cells_y;
    for(int var=0;var<nvar;var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<cells_y-1;j++){
                    for(int i=0;i<cv_x;i++){
                        cell = i + (j*n_fp)*cv_x + k*cv_x*fp_y        + i_ader*cv_x*fp_y*cv_z        + var*cv_x*fp_y*cv_z*n_cv;
                        face = i + j*cv_x        + k*cv_x*(cells_y-1) + i_ader*cv_x*(cells_y-1)*cv_z + var*cv_x*(cells_y-1)*cv_z*n_cv;
                        U_R[face] = U_ader_fp_y[cell+n_fp*cv_x];
                        U_L[face] = U_ader_fp_y[cell+n_cv*cv_x];
                    }
                }
            }
        }
    }
    //This should be a pointer towards whatever riemann solver was defined
    riemann_llf(F_y,U_L,U_R,cv_x*(cells_y-1)*cv_z*n_cv,_vy_,_vx_,_vz_);

    for(int var=0;var<nvar;var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<cells_y-1;j++){
                    for(int i=0;i<cv_x;i++){
                        cell = i + (j*n_fp)*cv_x + k*cv_x*fp_y        + i_ader*cv_x*fp_y*cv_z        + var*cv_x*fp_y*cv_z*n_cv;
                        face = i + j*cv_x        + k*cv_x*(cells_y-1) + i_ader*cv_x*(cells_y-1)*cv_z + var*cv_x*(cells_y-1)*cv_z*n_cv;
                        F_ader_fp_y[cell+n_fp*cv_x] = F_y[face];
                        F_ader_fp_y[cell+n_cv*cv_x] = F_y[face];
                    }
                }
            }
        }
    }
}

void riemann_solver_z(){
    //Store fluxes at the interface between elements to then solve the unique flux
    int cell,face;
    int fp_z=n_fp*cells_z;
    for(int var=0;var<nvar;var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int k=0;k<cells_z-1;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<cv_x;i++){
                        cell = i + j*cv_x + (k*n_fp)*cv_x*cv_y + i_ader*cv_x*cv_y*fp_z        + var*cv_x*cv_y*fp_z*n_cv;
                        face = i + j*cv_x +  k*cv_x*cv_y       + i_ader*cv_x*cv_y*(cells_z-1) + var*cv_x*cv_y*(cells_z-1)*n_cv;
                        U_R[face] = U_ader_fp_z[cell+n_fp*(cv_x*cv_y)];
                        U_L[face] = U_ader_fp_z[cell+n_cv*(cv_x*cv_y)];
                    }
                }
            }
        }
    }
    //This should be a pointer towards whatever riemann solver was defined
    riemann_llf(F_z,U_L,U_R,cv_x*cv_y*(cells_z-1)*n_cv,_vz_,_vy_,_vx_);

    for(int var=0;var<nvar;var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int k=0;k<cells_z-1;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<cv_x;i++){
                        cell = i + j*cv_x + (k*n_fp)*cv_x*cv_y + i_ader*cv_x*cv_y*fp_z        + var*cv_x*cv_y*fp_z*n_cv;
                        face = i + j*cv_x +  k*cv_x*cv_y       + i_ader*cv_x*cv_y*(cells_z-1) + var*cv_x*cv_y*(cells_z-1)*n_cv;
                        F_ader_fp_z[cell+n_fp*(cv_x*cv_y)] = F_z[face];
                        F_ader_fp_z[cell+n_cv*(cv_x*cv_y)] = F_z[face];
                    }
                }
            }
        }
    }
}
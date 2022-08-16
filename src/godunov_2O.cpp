#include "sd3d.hpp"

double minmod(double S_L, double S_R, double x_L, double x_R){
    double ratio,slope;
    //Compute ratio between slopes SlopeR/SlopeL
    ratio = S_R/S_L;
    //limit the ratio to in (0,1)
    ratio = max(0.0,min(ratio,1.0));
    //multiply by SlopeL to get the limited slope at the cell center
    slope = ratio*S_L;
    //Compute SlopeC*dx/2                             
    return 0.5*slope*(x_R-x_L);
}

void compute_prediction(double *dWt, double *W, double *dWx, double *dWy, double *dWz){
    int cv, size_cv;
    int ivx,ivy,ivz,ip;
    for(int k=0;k<cv_z;k++){
        for(int j=0;j<cv_y;j++){
            for(int i=0;i<cv_x;i++){
                cv = i + j*cv_x + k*cv_x*cv_y;
                ivx = cv+_vx_*size_cv;
                ivy = cv+_vy_*size_cv;
                ivz = cv+_vz_*size_cv;
                ip  = cv+ _p_*size_cv;
                dWt[cv] = 0;
                dWt[ip] = 0;
                #ifdef X
                dWt[cv] -= W[ivx]*dWx[cv] + W[cv]*dWx[ivx];
                dWt[ivx] = -W[ivx]*dWx[ivx] - dWx[ip]/W[cv];
                dWt[ip] -= W[ivx]*dWx[ip] + dWx[ivx]*gmma*W[ip];
                #ifdef Y
                dWt[ivx] -= W[ivy]*dWy[ivx];
                #endif
                #ifdef Z
                dWt[ivx] -= W[ivz]*dWz[ivx];
                #endif
                #endif

                #ifdef Y
                dWt[cv] -= W[ivy]*dWx[cv] + W[cv]*dWy[ivy];
                dWt[ivy] = -W[ivy]*dWy[ivy] - dWy[ip]/W[cv];
                dWt[ip] -= W[ivy]*dWy[ip] + dWy[ivy]*gmma*W[ip];
                #ifdef X
                dWt[ivy] -= W[ivx]*dWx[ivy];
                #endif
                #ifdef Z
                dWt[ivy] -= W[ivz]*dWz[ivy];
                #endif
                #endif
                
                
                #ifdef Z
                dWt[cv] -= W[ivz]*dWz[cv] + W[cv]*dWz[ivz];
                dWt[ivz] = -W[ivz]*dWz[ivz] - dWz[ip]/W[cv];
                dWt[ip] -= W[ivz]*dWz[ip] + dWz[ivz]*gamma*W[ip];
                #ifdef X
                dWt[ivz] -= W[ivx]*dWx[ivz];
                #endif
                #ifdef Y
                dWt[ivz] -= W[ivy]*dUy[ivz];
                #endif
                #endif
            }
        }
    }
}

void derivative_x(double *U, double *dUdx){
    int cv, nx=cv_x-1, face;
    double dx;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=0;i<nx;i++){
                    cv   = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                    face = i + j*nx   + k*nx*cv_y   + var*nx*cv_y*cv_z;
                    dx   = (X_centers[i+1]-X_centers[i]);
                    dUdx[face] = (U[cv+1] - U[cv])/dx;
                }
            }
        }  
    }
}

void derivative_y(double *U, double *dUdy){
    int cv, ny=cv_y-1, face;
    double dy;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=0;j<ny;j++){
                for(int i=0;i<cv_x;i++){
                    cv   = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                    face = i + j*cv_x + k*cv_x*ny   + var*cv_x*ny*cv_z;
                    dy   = (Y_centers[j+1]-Y_centers[j]);
                    dUdy[face] = (U[cv+cv_x] - U[cv])/dy;
                }
            }
        }  
    }
}

void derivative_z(double *U, double *dUdz){
    int cv, nz=cv_z-1, face;
    double dz;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<nz;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=0;i<cv_x;i++){
                    cv   = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                    face = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*nz;
                    dz   = (Z_centers[k+1]-Z_centers[k]);
                    dUdz[face] = (U[cv+cv_x*cv_y] - U[cv])/dz;
                }
            }
        }  
    }
}

void limit_slope_x(double *U, double* dU, double* dUc){
    double du,dx;
    int id, cv, face, nx=cv_x-1, faces=(Nx*n_cv)+1;
    int i_min = (NGHx*n_cv)-1;
    int i_max = (Nx+NGHx)*n_cv;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=i_min; i<=i_max;i++){
                    cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                    id = i-1 + j*nx + k*nx*cv_y + var*nx*cv_y*cv_z;
                    dUc[cv] = minmod(dU[id],dU[id+1],X_faces[i],X_faces[i+1]);
                }
            }
        }  
    }
}

void limit_slope_y(double *U, double* dU, double* dUc){
    double du, dy;
    int id, cv, face, ny=cv_y-1, faces=(Ny*n_cv)+1;
    int j_min = (NGHy*n_cv)-1;
    int j_max = (Ny+NGHy)*n_cv;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=j_min;j<=j_max;j++){
                for(int i=0;i<cv_x;i++){
                    cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                    id = i + (j-1)*cv_x + k*cv_x*ny + var*cv_x*ny*cv_z;
                    dUc[cv] = minmod(dU[id],dU[id+cv_x],Y_faces[j],Y_faces[j+1]);
                }
            }
        }  
    }
}

void limit_slope_z(double *U, double* dU, double* dUc){
    double du, dz;
    int id, cv, face, nz=cv_z-1, faces=(Nz*n_cv)+1;
    int k_min = (NGHz*n_cv)-1;
    int k_max = (Nz+NGHz)*n_cv;
    for(int var=0; var<nvar; var++){
        for(int k=k_min;k<=k_max;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=0;i<cv_x;i++){
                    cv = i + j*cv_x + k*cv_x*cv_y     + var*cv_x*cv_y*cv_z;
                    id = i + j*cv_x + (k-1)*cv_x*cv_y + var*cv_x*cv_y*nz;
                    dUc[cv] = minmod(dU[id],dU[id+cv_x*cv_y],Z_faces[k],Z_faces[k+1]);
                }
            }
        }  
    }
}
    

void interpolate_x(double *U, double* dUx, double* dUt, double* U_L, double* U_R){
    double du,dx;
    int id, cv, face, nx=cv_x-1, faces=(Nx*n_cv)+1;
    int i_min = (NGHx*n_cv)-1;
    int i_max = (Nx+NGHx)*n_cv;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=i_min; i<=i_max;i++){
                    cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                    face = (i-i_min) + j*faces + k*faces*cv_y + var*faces*cv_y*cv_z;
                    dx = X_faces[i+1]-X_faces[i];
                    if(i>i_min)
                        U_R[face-1] = U[cv] - dUx[cv] + dUt[cv]*dt/dx;
                    if(i<i_max)
                        U_L[face  ] = U[cv] + dUx[cv] + dUt[cv]*dt/dx;
                }
            }
        }  
    }
}

void interpolate_y(double *U, double* dUy, double* dUt, double* U_L, double* U_R){
    double du, dy;
    int id, cv, face, ny=cv_y-1, faces=(Ny*n_cv)+1;
    int j_min = (NGHy*n_cv)-1;
    int j_max = (Ny+NGHy)*n_cv;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=j_min;j<=j_max;j++){
                for(int i=0;i<cv_x;i++){
                    cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                    face = i + (j-j_min)*cv_x + k*cv_x*faces + var*cv_x*faces*cv_z;
                    dy = Y_faces[j+1]-Y_faces[j];
                    if(j>j_min)
                        U_R[face-cv_x] = U[cv] - dUy[cv] + dUt[cv]*dt/dy;
                    if(j<j_max)
                        U_L[face     ] = U[cv] + dUy[cv] + dUt[cv]*dt/dy;
                }
            }
        }  
    }
}

void interpolate_z(double *U, double* dUz, double* dUt, double* U_L, double* U_R){
    double du, dz;
    int id, cv, face, nz=cv_z-1, faces=(Nz*n_cv)+1;
    int k_min = (NGHz*n_cv)-1;
    int k_max = (Nz+NGHz)*n_cv;
    for(int var=0; var<nvar; var++){
        for(int k=k_min;k<=k_max;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=0;i<cv_x;i++){
                    cv = i + j*cv_x + k*cv_x*cv_y     + var*cv_x*cv_y*cv_z;
                    face = i + j*cv_x + (k-k_min)*cv_x*cv_y + var*cv_x*cv_y*faces;
                    dz = Z_faces[k+1] - Z_faces[k];
                    if(k>k_min)
                        U_R[face-cv_x*cv_y] = U[cv] - dUz[cv] + dUt[cv]*dt/dz;
                    if(k<k_max)
                        U_L[face          ] = U[cv] + dUz[cv] + dUt[cv]*dt/dz;
                }
            }
        }  
    }
}

void fv_godunov2O_update(double *U_new, double *U_cv, int i_ader){
    int cell, cv, fv;
    int face_x, fp_x=(Nx*n_cv+1);
    int face_y, fp_y=(Ny*n_cv+1);
    int face_z, fp_z=(Nz*n_cv+1);
    double dx,dy,dz;
    double F;
   for(int var=0; var<nvar; var++){
        for(int kk=NGHz;kk<Nz+NGHz;kk++){
            for(int jj=NGHy;jj<Ny+NGHy;jj++){
                for(int ii=NGHx;ii<Nx+NGHx;ii++){
                    cell   = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                    face_x = (ii-NGHx)*nx_cv + (jj*ny_cv)*fp_x        + (kk*nz_cv)*fp_x*cv_y        + var*fp_x*cv_y*cv_z;
                    face_y = (ii*nx_cv)      + ((jj-NGHy)*ny_cv)*cv_x + (kk*nz_cv)*cv_x*fp_y        + var*cv_x*fp_y*cv_z;
                    face_z = (ii*nx_cv)      + (jj*ny_cv)*cv_x        + ((kk-NGHz)*nz_cv)*cv_x*cv_y + var*cv_x*cv_y*fp_z;
                    for(int k=0;k<nz_cv;k++){
                        for(int j=0;j<ny_cv;j++){
                            for(int i=0;i<n_cv;i++){
                                cv=i + j*cv_x + k*cv_x*cv_y + cell;
                                F=0.0;
                                #ifdef X
                                fv = i + j*fp_x + k*fp_x*cv_y + face_x;
                                dx  = (X_faces[i+1 + ii*n_cv]-X_faces[i + ii*n_cv]);
                                F += (F_x[fv+1]-F_x[fv])/dx;
                                #endif
                                #ifdef Y
                                fv = i + j*cv_x + k*cv_x*fp_y + face_y;
                                dy  = (Y_faces[j+1 + jj*n_cv]-Y_faces[j+ jj*n_cv]);
                                F += (F_y[fv+cv_x]-F_y[fv])/dy;
                                #endif
                                #ifdef Z
                                fv = i + j*cv_x + k*cv_x*cv_y + face_z;
                                dz  = (Z_faces[k+1 + kk*n_cv]-Z_faces[k+ kk*n_cv]);
                                F += (F_z[fv+cv_x*cv_y]- F_z[fv])/dz;
                                #endif
                                U_new[cv] = U_cv[cv] - w_t[i_ader]*dt*F;
                            }
                        }
                    }
                }
            }
        }
    }
}

void godunov_2O(){
    #ifdef X
    derivative_x(W_cv,U_L);
    limit_slope_x(W_cv,U_L,dUdx);
    #endif
    #ifdef Y
    derivative_y(W_cv,U_L);
    limit_slope_y(W_cv,U_L,dUdy);
    #endif
    #ifdef Z
    derivative_z(W_cv,U_L);
    limit_slope_z(W_cv,U_L,dUdz);
    #endif
    compute_prediction(dUdt,W_cv,dUdx,dUdy,dUdz);
    #ifdef X
    interpolate_x(W_cv, dUdx, dUdt, U_L, U_R);
    //riemann_llf(F_x,U_L,U_R,(Nx*nx_cv+1)*cv_y*cv_z,_vx_,_vy_,_vz_);
    riemann_hllc(F_x,U_L,U_R,(Nx*nx_cv+1)*cv_y*cv_z,_vx_,_vy_,_vz_,0);
    #endif
    #ifdef Y
    interpolate_y(W_cv, dUdy, dUdt, U_L, U_R);
    //riemann_llf(F_y,U_L,U_R,cv_x*(Ny*ny_cv+1)*cv_z,_vy_,_vx_,_vz_);
    riemann_hllc(F_y,U_L,U_R,cv_x*(Ny*ny_cv+1)*cv_z,_vy_,_vx_,_vz_,0);
    #endif
    #ifdef Z
    interpolate_z(W_cv, dUdz, dUdt, U_L, U_R);
    //riemann_llf(F_z,U_L,U_R,cv_x*cv_y*(Nz*nz_cv+1),_vz_,_vx_,_vy_);
    riemann_hllc(F_z,U_L,U_R,cv_x*cv_y*(Nz*nz_cv+1),_vz_,_vx_,_vy_,0);
    #endif
}

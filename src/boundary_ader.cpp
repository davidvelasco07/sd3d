#include "sd3d.hpp"

void Boundary_Conditions_ader(double* U){
    #ifdef X
    #ifdef Periodic_X
    store_boundaries_ader_x(U);
    #endif
    if(cpu_x>1)Exec_comms_ader_x(U);
    apply_boundaries_ader_x(U);
    #endif
    #ifdef Y
    #ifdef Periodic_Y
    store_boundaries_ader_y(U);
    #endif
    if(cpu_y>1)Exec_comms_ader_y(U);
    apply_boundaries_ader_y(U);
    #endif
    #ifdef Z
    #ifdef Periodic_Z
    store_boundaries_ader_z(U);
    #endif
    if(cpu_z>1)Exec_comms_ader_z(U);
    apply_boundaries_ader_z(U);
    #endif
}

void Store_boundaries_ader(double *U){
    #ifdef X
    store_boundaries_ader_x(U);
    #endif
    #ifdef Y
    store_boundaries_ader_y(U);
    #endif
    #ifdef Z
    store_boundaries_ader_z(U);
    #endif
}

void store_boundaries_ader_x(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_cv; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<NGH*n_cv;i++){
                        cv = i + j*NGH*n_cv + k*NGH*n_cv*cv_y + i_a*NGH*n_cv*cv_y*cv_z + var*NGH*n_cv*cv_y*cv_z*n_cv;
                        #ifdef Periodic_X
                        BC_ader_x[0][cv] = U[i+Nx*n_cv  + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv]; 
                        BC_ader_x[1][cv] = U[i+NGH*n_cv + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv]; 
                        #else
                        BC_ader_x[0][cv] = U[i               + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv]; 
                        BC_ader_x[1][cv] = U[i+(Nx+NGH)*n_cv + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv];
                        #endif
                    }
                }
            }
        }
    }
}

void store_boundaries_ader_y(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_cv; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<NGH*n_cv;j++){
                    for(int i=0;i<cv_x;i++){
                        cv = i + j*cv_x + k*cv_x*NGH*n_cv + i_a*cv_x*NGH*n_cv*cv_z + var*cv_x*NGH*n_cv*cv_z*n_cv;
                        #ifdef Periodic_Y
                        BC_ader_y[0][cv] = U[i + (j+Ny*n_cv)*cv_x  + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv]; 
                        BC_ader_y[1][cv] = U[i + (j+NGH*n_cv)*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv]; 
                        #else
                        BC_ader_y[0][cv] = U[i +  j*cv_x                + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv]; 
                        BC_ader_y[1][cv] = U[i + (j+(NGH+Ny)*n_cv)*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv]; 
                        #endif
                    }  
                }
            }
        }
    }
}

void store_boundaries_ader_z(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_cv; i_a++){
            for(int k=0;k<NGH*n_cv;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<cv_x;i++){
                        cv = i + j*cv_x + k*cv_x*cv_y + i_a*cv_x*cv_y*NGH*n_cv + var*cv_x*cv_y*NGH*n_cv*n_cv;
                        #ifdef Periodic_Z
                        BC_ader_z[0][cv] = U[i + j*cv_x + (k+Nz*n_cv)*cv_x*cv_y  + i_a*size_cv + var*size_cv*n_cv]; 
                        BC_ader_z[1][cv] = U[i + j*cv_x + (k+NGH*n_cv)*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv]; 
                        #else
                        BC_ader_z[0][cv] = U[i + j*cv_x +  k*cv_x*cv_y                + i_a*size_cv + var*size_cv*n_cv]; 
                        BC_ader_z[1][cv] = U[i + j*cv_x + (k+(NGH+Nz)*n_cv)*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv]; 
                        #endif  
                    }
                }
            }
        }
    }
}

void apply_boundaries_ader_x(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_cv; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<NGH*n_cv;i++){
                        cv = i + j*NGH*n_cv + k*NGH*n_cv*cv_y + i_a*NGH*n_cv*cv_y*cv_z + var*NGH*n_cv*cv_y*cv_z*n_cv;
                        U[i               + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv] = BC_ader_x[0][cv]; 
                        U[i+(Nx+NGH)*n_cv + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv] = BC_ader_x[1][cv]; 
                    }
                }
            }
        }
    }
}

void apply_boundaries_ader_y(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_cv; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<NGH*n_cv;j++){
                    for(int i=0;i<cv_x;i++){
                        cv = i + j*cv_x + k*cv_x*NGH*n_cv + i_a*cv_x*NGH*n_cv*cv_z + var*cv_x*NGH*n_cv*cv_z*n_cv;
                        U[i +  j*cv_x                + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv] = BC_ader_y[0][cv];
                        U[i + (j+(Ny+NGH)*n_cv)*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv] = BC_ader_y[1][cv];
                    }
                }
            }
        }
    }
}

void apply_boundaries_ader_z(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_cv; i_a++){
            for(int k=0;k<NGH*n_cv;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<cv_x;i++){
                        cv = i + j*cv_x + k*cv_x*cv_y + i_a*cv_x*cv_y*NGH*n_cv + var*cv_x*cv_y*NGH*n_cv*n_cv;
                        U[i + j*cv_x + k*cv_x*cv_y                 + i_a*size_cv + var*size_cv*n_cv] = BC_ader_z[0][cv];
                        U[i + j*cv_x + (k+(Nz+NGH)*n_cv)*cv_x*cv_y + i_a*size_cv + var*size_cv*n_cv] = BC_ader_z[1][cv];
                    }
                }
            }
        }
    }
}
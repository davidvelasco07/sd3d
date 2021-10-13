#include "sd3d.hpp"

void Boundary_Conditions(double* U){
    //Boundaries
    #if defined(X) && defined(Periodic_X)
    store_boundaries_x(U);
    #endif
    #if defined(Y) && defined(Periodic_Y)
    store_boundaries_y(U);
    #endif
    #if defined(Z) && defined(Periodic_Z)
    store_boundaries_z(U);
    #endif

    if(cpu_size>1)Exec_comms(U);
        
    #ifdef X
    apply_boundaries_x(U);
    #endif
    #ifdef Y
    apply_boundaries_y(U);
    #endif
    #ifdef Z
    apply_boundaries_z(U);
    #endif
}

void Store_boundaries(double *U){
    #ifdef X
    store_boundaries_x(U);
    #endif
    #ifdef Y
    store_boundaries_y(U);
    #endif
    #ifdef Z
    store_boundaries_z(U);
    #endif
}

void store_boundaries_x(double *U){
    int cell;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=0;i<NGH*n_cv;i++){
                    cell = i + j*NGH*n_cv + k*NGH*n_cv*cv_y + var*NGH*n_cv*cv_y*cv_z;
                    #ifdef Periodic_X
                    BC_x[0][cell] = U[i+Nx*n_cv  + j*cv_x + k*cv_x*cv_y + var*size_cv]; 
                    BC_x[1][cell] = U[i+NGH*n_cv + j*cv_x + k*cv_x*cv_y + var*size_cv]; 
                    #else
                    BC_x[0][cell] = U[i               + j*cv_x + k*cv_x*cv_y + var*size_cv]; 
                    BC_x[1][cell] = U[i+(Nx+NGH)*n_cv + j*cv_x + k*cv_x*cv_y + var*size_cv];
                    #endif
                }
            }
        }
    }
}

void store_boundaries_y(double *U){
    int cell;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=0;j<NGH*n_cv;j++){
                for(int i=0;i<cv_x;i++){
                    cell = i + j*cv_x + k*cv_x*NGH*n_cv + var*cv_x*NGH*n_cv*cv_z;
                    #ifdef Periodic_Y
                    BC_y[0][cell] = U[i + (j+Ny*n_cv)*cv_x  + k*cv_x*cv_y + var*size_cv]; 
                    BC_y[1][cell] = U[i + (j+NGH*n_cv)*cv_x + k*cv_x*cv_y + var*size_cv]; 
                    #else
                    BC_y[0][cell] = U[i +  j*cv_x               + k*cv_x*cv_y + var*size_cv]; 
                    BC_y[1][cell] = U[i + (j+(NGH+Ny)*n_cv)*cv_x + k*cv_x*cv_y + var*size_cv]; 
                    #endif  
                }
            }
        }
    }
}

void store_boundaries_z(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<NGH*n_cv;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=0;i<cv_x;i++){
                    cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*NGH*n_cv;
                    #ifdef Periodic_Z
                    BC_z[0][cv] = U[i + j*cv_x + (k+Nz*n_cv)*cv_x*cv_y  + var*size_cv]; 
                    BC_z[1][cv] = U[i + j*cv_x + (k+NGH*n_cv)*cv_x*cv_y + var*size_cv]; 
                    #else
                    BC_z[0][cv] = U[i + j*cv_x +  k*cv_x*cv_y                + var*size_cv]; 
                    BC_z[1][cv] = U[i + j*cv_x + (k+(NGH+Nz)*n_cv)*cv_x*cv_y + var*size_cv]; 
                    #endif  
                }
            }
        }
    }
}

void apply_boundaries_x(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=0;i<NGH*n_cv;i++){
                    cv = i + j*NGH*n_cv + k*NGH*n_cv*cv_y + var*NGH*n_cv*cv_y*cv_z;
                    U[i               + j*cv_x + k*cv_x*cv_y + var*size_cv] = BC_x[0][cv]; 
                    U[i+(Nx+NGH)*n_cv + j*cv_x + k*cv_x*cv_y + var*size_cv] = BC_x[1][cv]; 
                }
            }
        }
    }
}

void apply_boundaries_y(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=0;j<NGH*n_cv;j++){
                for(int i=0;i<cv_x;i++){
                    cv = i + j*cv_x + k*cv_x*NGH*n_cv + var*cv_x*NGH*n_cv*cv_z;
                    U[i + j*cv_x                 + k*cv_x*cv_y + var*size_cv] = BC_y[0][cv];
                    U[i + (j+(Ny+NGH)*n_cv)*cv_x + k*cv_x*cv_y + var*size_cv] = BC_y[1][cv];
                }
            }
        }
    }
}

void apply_boundaries_z(double *U){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<NGH*n_cv;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=0;i<cv_x;i++){
                    cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*NGH*n_cv;
                    U[i + j*cv_x + k*cv_x*cv_y                 + var*size_cv] = BC_z[0][cv];
                    U[i + j*cv_x + (k+(Nz+NGH)*n_cv)*cv_x*cv_y + var*size_cv] = BC_z[1][cv];
                }
            }
        }
    }
}
#include "sd3d.hpp"

void Boundary_Conditions(double* U, int n_a){
    //Boundaries
    
    #if defined(X) && defined(Periodic_X)
    store_boundaries_x(U,n_a);
    #endif
    #if defined(Y) && defined(Periodic_Y)
    store_boundaries_y(U,n_a);
    #endif
    #if defined(Z) && defined(Periodic_Z)
    store_boundaries_z(U,n_a);
    #endif
   
    if(cpu_size>1)Exec_comms(U,n_a);
        
    #ifdef X
    apply_boundaries_x(U,n_a);
    #endif
    #ifdef Y
    apply_boundaries_y(U,n_a);
    #endif
    #ifdef Z
    apply_boundaries_z(U,n_a);
    #endif
}

void Store_boundaries(double *U){
    #ifdef X
    store_boundaries_x(U,1);
    #endif
    #ifdef Y
    store_boundaries_y(U,1);
    #endif
    #ifdef Z
    store_boundaries_z(U,1);
    #endif
}

void store_boundaries_x(double *U, int n_a){
    int cell;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_a; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<NGH*n_cv;i++){
                        cell = i + j*NGH*n_cv + k*NGH*n_cv*cv_y + i_a*NGH*n_cv*cv_y*cv_z + var*NGH*n_cv*cv_y*cv_z*n_a;
                        #ifdef Periodic_X
                        BC_x[0][cell] = U[i+Nx*n_cv  + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                        BC_x[1][cell] = U[i+NGH*n_cv + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                        #else
                        BC_x[0][cell] = U[i               + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                        BC_x[1][cell] = U[i+(Nx+NGH)*n_cv + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a];
                        #endif
                    }
                }
            }
        }
    }
}

void store_boundaries_y(double *U,int n_a){
    int cell;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_a; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<NGH*n_cv;j++){
                    for(int i=0;i<cv_x;i++){
                        cell = i + j*cv_x + k*cv_x*NGH*n_cv + i_a*cv_x*NGH*n_cv*cv_z + var*cv_x*NGH*n_cv*cv_z*n_a;
                        #ifdef Periodic_Y
                        BC_y[0][cell] = U[i + (j+Ny*n_cv)*cv_x  + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                        BC_y[1][cell] = U[i + (j+NGH*n_cv)*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                        #else
                        BC_y[0][cell] = U[i +  j*cv_x                + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                        BC_y[1][cell] = U[i + (j+(NGH+Ny)*n_cv)*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                        #endif
                    }  
                }
            }
        }
    }
}

void store_boundaries_z(double *U, int n_a){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_a; i_a++){
            for(int k=0;k<NGH*n_cv;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<cv_x;i++){
                        cv = i + j*cv_x + k*cv_x*cv_y + i_a*cv_x*cv_y*NGH*n_cv + var*cv_x*cv_y*NGH*n_cv*n_a;
                        #ifdef Periodic_Z
                        BC_z[0][cv] = U[i + j*cv_x + (k+Nz*n_cv)*cv_x*cv_y  + i_a*size_cv + var*size_cv*n_a]; 
                        BC_z[1][cv] = U[i + j*cv_x + (k+NGH*n_cv)*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                        #else
                        BC_z[0][cv] = U[i + j*cv_x +  k*cv_x*cv_y                + i_a*size_cv + var*size_cv*n_a]; 
                        BC_z[1][cv] = U[i + j*cv_x + (k+(NGH+Nz)*n_cv)*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                        #endif  
                    }
                }
            }
        }
    }
}

void apply_boundaries_x(double *U, int n_a){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_a; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<NGH*n_cv;i++){
                        cv = i + j*NGH*n_cv + k*NGH*n_cv*cv_y + i_a*NGH*n_cv*cv_y*cv_z + var*NGH*n_cv*cv_y*cv_z*n_a;
                        U[i               + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a] = BC_x[0][cv]; 
                        U[i+(Nx+NGH)*n_cv + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a] = BC_x[1][cv]; 
                    }
                }
            }
        }
    }
}

void apply_boundaries_y(double *U, int n_a){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_a; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<NGH*n_cv;j++){
                    for(int i=0;i<cv_x;i++){
                        cv = i + j*cv_x + k*cv_x*NGH*n_cv + i_a*cv_x*NGH*n_cv*cv_z + var*cv_x*NGH*n_cv*cv_z*n_a;
                        U[i + j*cv_x                 + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a] = BC_y[0][cv];
                        U[i + (j+(Ny+NGH)*n_cv)*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a] = BC_y[1][cv];
                    }
                }
            }
        }
    }
}

void apply_boundaries_z(double *U, int n_a){
    int cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_a; i_a++){
            for(int k=0;k<NGH*n_cv;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<cv_x;i++){
                        cv = i + j*cv_x + k*cv_x*cv_y + i_a*cv_x*cv_y*NGH*n_cv + var*cv_x*cv_y*NGH*n_cv*n_a;
                        U[i + j*cv_x + k*cv_x*cv_y                 + i_a*size_cv + var*size_cv*n_a] = BC_z[0][cv];
                        U[i + j*cv_x + (k+(Nz+NGH)*n_cv)*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a] = BC_z[1][cv];
                    }
                }
            }
        }
    }
}
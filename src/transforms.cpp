#include "sd3d.hpp"

void transform_cv_to_sp(double *U_cv, double *U_sp){
    int cell;
    double s;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int kk=0;kk<cells_z;kk++){
            for(int jj=0;jj<cells_y;jj++){
                for(int ii=0;ii<cells_x;ii++){
                    cell = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + var*size_cv;
                    for(int k=0;k<nz_cv;k++){
                        for(int j=0;j<ny_cv;j++){
                            for(int i=0;i<nx_cv;i++){
                                U_sp[i+j*cv_x+k*cv_x*cv_y+cell]=0.0;
                                for(int n=0;n<nz_cv;n++){
                                    for(int m=0;m<ny_cv;m++){
                                        for(int l=0;l<nx_cv;l++){
                                            s = U_cv[l+m*cv_x+n*cv_x*cv_y+cell];
                                            #ifdef X
                                            s*=cv_to_sp[l+i*n_cv];
                                            #endif
                                            #ifdef Y
                                            s*=cv_to_sp[m+j*n_cv];
                                            #endif
                                            #ifdef Z
                                            s*=cv_to_sp[n+k*n_cv];
                                            #endif
                                            U_sp[i+j*cv_x+k*cv_x*cv_y+cell] += s;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void transform_sp_to_cv(double *U_sp, double *U_cv){
    int cell;
    double s;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int kk=0;kk<cells_z;kk++){
            for(int jj=0;jj<cells_y;jj++){
                for(int ii=0;ii<cells_x;ii++){
                    cell = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + var*size_cv;
                    for(int k=0;k<nz_cv;k++){
                        for(int j=0;j<ny_cv;j++){
                            for(int i=0;i<nx_cv;i++){
                                U_cv[i+j*cv_x+k*cv_x*cv_y+cell]=0.0;
                                for(int n=0;n<nz_cv;n++){
                                    for(int m=0;m<ny_cv;m++){
                                        for(int l=0;l<nx_cv;l++){
                                            s = U_sp[l+m*cv_x+n*cv_x*cv_y+cell];
                                            #ifdef X
                                            s*=sp_to_cv[l+i*n_cv];
                                            #endif
                                            #ifdef Y
                                            s*=sp_to_cv[m+j*n_cv];
                                            #endif
                                            #ifdef Z
                                            s*=sp_to_cv[n+k*n_cv];
                                            #endif
                                            U_cv[i+j*cv_x+k*cv_x*cv_y+cell] += s;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void transform_sp_to_fp_x(double *U_sp, double *U_fp_x, int n_ader){
    int cell,face,Nfp=n_fp*cells_x;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int i_ader = 0; i_ader<n_ader; i_ader++){
            for(int kk=0;kk<cells_z;kk++){
                for(int jj=0;jj<cells_y;jj++){
                    for(int ii=0;ii<cells_x;ii++){
                        cell = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + i_ader*size_cv       + var*size_cv*n_ader;
                        face = (ii*n_fp) + (jj*ny_cv)*Nfp  + (kk*nz_cv)*Nfp*cv_y  + i_ader*Nfp*cv_y*cv_z + var*Nfp*cv_y*cv_z*n_ader;
                        for(int k=0;k<nz_cv;k++){
                            for(int j=0;j<ny_cv;j++){
                                for(int i=0;i<n_fp;i++){
                                    U_fp_x[i+j*Nfp+k*Nfp*cv_y+face]=0.0;
                                    for(int l=0;l<n_cv;l++){
                                        U_fp_x[i+j*Nfp+k*Nfp*cv_y+face] += U_sp[l+j*cv_x+k*cv_x*cv_y+cell]*sp_to_fp[l+i*n_cv];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void transform_sp_to_fp_y(double *U_sp, double *U_fp_y, int n_ader){
    int cell,face,fp_y=n_fp*cells_y;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int i_ader = 0; i_ader<n_ader; i_ader++){
            for(int kk=0;kk<cells_z;kk++){
                for(int jj=0;jj<cells_y;jj++){
                    for(int ii=0;ii<cells_x;ii++){
                        cell = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + i_ader*size_cv        + var*size_cv*n_ader;
                        face = (ii*nx_cv) + (jj*ny_fp)*cv_x + (kk*nz_cv)*cv_x*fp_y + i_ader*cv_x*fp_y*cv_z + var*cv_x*fp_y*cv_z*n_ader;
                        for(int k=0;k<nz_cv;k++){
                            for(int j=0;j<n_fp;j++){
                                for(int i=0;i<nx_cv;i++){
                                    U_fp_y[i+j*cv_x+k*cv_x*fp_y+face]=0.0;
                                    for(int m=0;m<n_cv;m++){
                                        U_fp_y[i+j*cv_x+k*cv_x*fp_y+face] += U_sp[i+m*cv_x+k*cv_x*cv_y+cell]*sp_to_fp[m+j*n_cv]; 
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void transform_sp_to_fp_z(double *U_sp, double *U_fp_z, int n_ader){
    int cell,face,Nfp=n_fp*cells_z;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int i_ader = 0; i_ader<n_ader; i_ader++){
            for(int kk=0;kk<cells_z;kk++){
                for(int jj=0;jj<cells_y;jj++){
                    for(int ii=0;ii<cells_x;ii++){
                        cell = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + i_ader*size_cv       + var*size_cv*n_ader ;
                        face = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*n_fp)*cv_x*cv_y  + i_ader*cv_x*cv_y*Nfp + var*cv_x*cv_y*Nfp*n_ader ;
                        for(int k=0;k<n_fp;k++){
                            for(int j=0;j<ny_cv;j++){
                                for(int i=0;i<nx_cv;i++){
                                    U_fp_z[i+j*cv_x+k*cv_x*cv_y+face]=0.0;
                                    for(int l=0;l<n_cv;l++){
                                        U_fp_z[i+j*cv_x+k*cv_x*cv_y+face] += U_sp[i+j*cv_x+l*cv_x*cv_y+cell]*sp_to_fp[l+k*n_cv];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }   
}

void derive_fp_x_to_sp(double *F_fp_x, double *U_sp){
    int cell,face,Nfp=n_fp*cells_x;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int kk=0;kk<cells_z;kk++){
                for(int jj=0;jj<cells_y;jj++){
                    for(int ii=0;ii<cells_x;ii++){
                        cell = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + i_ader*size_cv       + var*size_cv*n_cv;
                        face = (ii*n_fp)  + (jj*ny_cv)*Nfp  + (kk*nz_cv)*Nfp*cv_y  + i_ader*Nfp*cv_y*cv_z + var*Nfp*cv_y*cv_z*n_cv;
                        for(int k=0;k<nz_cv;k++){
                            for(int j=0;j<ny_cv;j++){
                                for(int i=0;i<nx_cv;i++){
                                    U_sp[i+j*cv_x+k*cv_x*cv_y+cell]=0.0;
                                    for(int l=0;l<n_fp;l++){
                                        U_sp[i+j*cv_x+k*cv_x*cv_y+cell] += (F_fp_x[l+j*Nfp+k*Nfp*cv_y+face]*dfp_to_sp[l+i*n_fp])/dx;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void derive_fp_y_to_sp(double *F_fp_y, double *U_sp){
    int cell,face,fp_y=n_fp*cells_y;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int kk=0;kk<cells_z;kk++){
                for(int jj=0;jj<cells_y;jj++){
                    for(int ii=0;ii<cells_x;ii++){
                        cell = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + i_ader*size_cv       + var*size_cv*n_cv;
                        face = (ii*nx_cv) + (jj*ny_fp)*cv_x + (kk*nz_cv)*cv_x*fp_y + i_ader*cv_x*fp_y*cv_z + var*cv_x*fp_y*cv_z*n_cv;
                        for(int k=0;k<nz_cv;k++){
                            for(int j=0;j<ny_cv;j++){
                                for(int i=0;i<nx_cv;i++){
                                    #ifndef X
                                    U_sp[i+j*cv_x+k*cv_x*cv_y+cell]=0.0;
                                    #endif
                                    for(int m=0;m<n_fp;m++){
                                        U_sp[i+j*cv_x+k*cv_x*cv_y+cell] += (F_fp_y[i+m*cv_x+k*cv_x*fp_y+face]*dfp_to_sp[m+j*n_fp])/dy;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void derive_fp_z_to_sp(double *F_fp_z, double *U_sp){
    int cell,face,Nfp=n_fp*cells_z;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int kk=0;kk<cells_z;kk++){
                for(int jj=0;jj<cells_y;jj++){
                    for(int ii=0;ii<cells_x;ii++){
                        cell = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + i_ader*size_cv       + var*size_cv*n_cv;
                        face = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*n_fp)*cv_x*cv_y   + i_ader*cv_x*cv_y*Nfp + var*cv_x*cv_y*Nfp*n_cv;
                        for(int k=0;k<nz_cv;k++){
                            for(int j=0;j<ny_cv;j++){
                                for(int i=0;i<nx_cv;i++){
                                    #if !defined(X) && !defined(Y)
                                    U_sp[i+j*cv_x+k*cv_x*cv_y+cell]=0.0;
                                    #endif
                                    for(int l=0;l<n_fp;l++){
                                        U_sp[i+j*cv_x+k*cv_x*cv_y+cell] += (F_fp_z[i+j*cv_x+l*cv_x*cv_y+face]*dfp_to_sp[l+k*n_fp])/dz;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void face_integral_x(double *F_fp_x, double *F_fv_x, int i_a){
    int cell,face,face_ader,fp_x=nx_fp*cells_x,fv_x=cv_faces_x;
    double s;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int kk=0;kk<cells_z;kk++){
            for(int jj=0;jj<cells_y;jj++){
                for(int ii=0;ii<cells_x;ii++){
                    face      = (ii*nx_cv)  + (jj*ny_cv)*fv_x  + (kk*nz_cv)*fv_x*cv_y  + var*fv_x*cv_y*cv_z;
                    face_ader = (ii*nx_fp)  + (jj*ny_cv)*fp_x  + (kk*nz_cv)*fp_x*cv_y  + i_a*fp_x*cv_y*cv_z + var*fp_x*cv_y*cv_z*n_cv;
                    for(int k=0;k<nz_cv;k++){
                        for(int j=0;j<ny_cv;j++){
                            for(int i=0;i<nx_cv;i++){
                                F_fv_x[i + j*fv_x + k*fv_x*cv_y + face]=0.0;
                                for(int n=0;n<nz_cv;n++){
                                    for(int m=0;m<ny_cv;m++){
                                        s = F_fp_x[i + m*fp_x + n*fp_x*cv_y + face_ader];
                                        #ifdef Y
                                        s*=sp_to_cv[m+j*n_cv];
                                        #endif
                                        #ifdef Z
                                        s*=sp_to_cv[n+k*n_cv];
                                        #endif
                                        F_fv_x[i + j*fv_x + k*fv_x*cv_y + face] += s;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void face_integral_y(double *F_fp_y, double *F_fv_y, int i_a){
    int cell,face,face_ader,fp_y=ny_fp*cells_y,fv_y=cv_faces_y;
    double s;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int kk=0;kk<cells_z;kk++){
            for(int jj=0;jj<cells_y;jj++){
                for(int ii=0;ii<cells_x;ii++){
                    face      = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*fv_y + var*cv_x*fv_y*cv_z;
                    face_ader = (ii*nx_cv) + (jj*ny_fp)*cv_x + (kk*nz_cv)*cv_x*fp_y + i_a*cv_x*fp_y*cv_z + var*cv_x*fp_y*cv_z*n_cv;
                    for(int k=0;k<nz_cv;k++){
                        for(int j=0;j<ny_cv;j++){
                            for(int i=0;i<nx_cv;i++){
                                F_fv_y[i+j*cv_x+k*cv_x*fv_y+face]=0.0;
                                for(int n=0;n<nz_cv;n++){
                                    for(int l=0;l<nx_cv;l++){
                                        s = F_fp_y[l+j*cv_x+n*cv_x*fp_y+face_ader];
                                        #ifdef X
                                        s*=sp_to_cv[l+i*n_cv];
                                        #endif
                                        #ifdef Z
                                        s*=sp_to_cv[n+k*n_cv];
                                        #endif
                                        F_fv_y[i+j*cv_x+k*cv_x*fv_y+face] += s;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void face_integral_z(double *F_fp_z, double *F_fv_z, int i_a){
    int cell,face,face_ader,fp_z=nz_fp*cells_z,fv_z=cv_faces_z;
    double s;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int kk=0;kk<cells_z;kk++){
            for(int jj=0;jj<cells_y;jj++){
                for(int ii=0;ii<cells_x;ii++){
                    face      = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_cv)*cv_x*cv_y + var*cv_x*cv_y*fv_z;
                    face_ader = (ii*nx_cv) + (jj*ny_cv)*cv_x + (kk*nz_fp)*cv_x*cv_y + i_a*cv_x*cv_y*fp_z + var*cv_x*cv_y*fp_z*n_cv;
                    for(int k=0;k<nz_cv;k++){
                        for(int j=0;j<ny_cv;j++){
                            for(int i=0;i<nx_cv;i++){
                                F_fv_z[i+j*cv_x+k*cv_x*cv_y+face]=0.0;
                                for(int m=0;m<ny_cv;m++){
                                    for(int l=0;l<nx_cv;l++){
                                        s = F_fp_z[l+m*cv_x+k*cv_x*cv_y+face_ader];
                                        #ifdef X
                                        s*=sp_to_cv[l+i*n_cv];
                                        #endif
                                        #ifdef Y
                                        s*=sp_to_cv[m+j*n_cv];
                                        #endif
                                        F_fv_z[i+j*cv_x+k*cv_x*cv_y+face] += s;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void volume_integral(double *U_fp, double *U_sp){
    int cv,corner;
    int fp_x=cells_x*nx_fp;
    int fp_y=cells_y*ny_fp;
    int fp_z=cells_z*nz_fp;
    double s;
    //Get values at solution points from control volume averages
    for(int var=0; var<nvar; var++){
        for(int kk=0;kk<cells_z;kk++){
            for(int jj=0;jj<cells_y;jj++){
                for(int ii=0;ii<cells_x;ii++){
                    corner = (ii*nx_fp)  + (jj*ny_fp)*fp_x  + (kk*nz_fp)*fp_x*fp_y  + var*fp_x*fp_y*fp_z;
                    cv     = (ii*nx_cv)  + (jj*ny_cv)*cv_x  + (kk*nz_cv)*cv_x*cv_y  + var*cv_x*cv_y*cv_z;
                    for(int k=0;k<nz_cv;k++){
                        for(int j=0;j<ny_cv;j++){
                            for(int i=0;i<n_cv;i++){
                                U_sp[i + j*cv_x + k*cv_x*cv_y + cv]=0.0;
                                for(int n=0;n<nz_fp;n++){
                                    for(int m=0;m<ny_fp;m++){
                                        for(int l=0;l<nx_fp;l++){
                                            s = U_fp[l + m*fp_x + n*fp_x*fp_y + corner];
                                            #ifdef X
                                            s*=fp_to_sp[l+i*n_fp];
                                            #endif
                                            #ifdef Y
                                            s*=fp_to_sp[m+j*n_fp];
                                            #endif
                                            #ifdef Z
                                            s*=fp_to_sp[n+k*n_fp];
                                            #endif
                                            U_sp[i + j*cv_x + k*cv_x*cv_y + cv] += s;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
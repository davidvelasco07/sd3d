#include "sd3d.hpp"

void ader_subupdate(double *U_ader_sp, double *U_sp, double *dU_ader_sp){
    int cell;
    for(int var=0; var<nvar; var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int i=0;i<size_cv;i++){
                cell=i+var*size_cv*n_cv;
                U_ader_sp[cell+i_ader*size_cv] = U_sp[i+var*size_cv];
                for(int l=0;l<n_cv;l++){
                    U_ader_sp[cell+i_ader*size_cv] -= (dU_ader_sp[cell+l*size_cv]*invader[l+i_ader*n_cv]*w_t[l])*dt;
                }
            }
        }
    }
}

void ader_update(double *U_sp, double *dU_sp){
    int cell;
    for(int var=0; var<nvar; var++){
        for(int i_ader = 0; i_ader<n_cv; i_ader++){
            for(int i=0;i<size_cv;i++){
                U_sp[i+var*size_cv] -= dU_sp[i + i_ader*size_cv + var*size_cv*n_cv]*w_t[i_ader]*dt;
            }
        }
    }
}              

void fv_update(double *U_new, double *U_cv, int i_ader){
    int cell, cv;
    int face_x, fv_x=cv_faces_x;
    int face_y, fv_y=cv_faces_y;
    int face_z, fv_z=cv_faces_z;
    double dx,dy,dz;
    double f;
    for(int var=0; var<nvar; var++){
        for(int k=0;k<cv_z;k++){
            for(int j=0;j<cv_y;j++){
                for(int i=0;i<cv_x;i++){
                    f=0.0;
                    #ifdef X
                    dx  = (X_faces[i+1]-X_faces[i]);
                    face_x = i + j*fv_x + k*fv_x*cv_y + var*fv_x*cv_y*cv_z;
                    f += (F_fv_x[face_x+1] - F_fv_x[face_x])/dx;
                    #endif
                    #ifdef Y
                    dy  = (Y_faces[j+1]-Y_faces[j]);
                    face_y = i + j*cv_x + k*cv_x*fv_y + var*cv_x*fv_y*cv_z;
                    f += (F_fv_y[face_y+cv_x]- F_fv_y[face_y])/dy;
                    #endif
                    #ifdef Z
                    dz  = (Z_faces[k+1]-Z_faces[k]);
                    face_z = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*fv_z;
                    f += (F_fv_z[face_z+cv_x*cv_y]- F_fv_z[face_z])/dz;
                    #endif
                    cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                    U_new[cv] = U_cv[cv] - w_t[i_ader]*dt*f;
                }
            }
        }
    }
}
               
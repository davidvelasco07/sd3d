#include "sd3d.hpp"

void fv_corrected_update(double *U_new, double *U_cv, int i_ader){
    int cv, cv_U, fv, fc, ngh=NGH*n_cv;
    int face_x, fv_x=cv_faces_x, fc_x=(Nx*nx_cv)+1;
    int face_y, fv_y=cv_faces_y, fc_y=(Ny*ny_cv)+1;
    int face_z, fv_z=cv_faces_z, fc_z=(Nz*nz_cv)+1;
    double dx,dy,dz;
    double F, Fp, Fm;
    double affected=0;
    for(int var=0; var<nvar; var++){
        for(int k=NGHz*nz_cv;k<(Nz+NGHz)*nz_cv;k++){
            for(int j=NGHy*ny_cv;j<(Ny+NGHy)*ny_cv;j++){
                for(int i=NGHx*nx_cv;i<(Nx+NGHx)*nx_cv;i++){
                    cv   = i + j*cv_x + k*cv_x*cv_y;
                    cv_U = cv + var*cv_x*cv_y*cv_z;
                    F=0.0;
                    #ifdef X
                    dx  = (X_faces[i+1]-X_faces[i]);
                    fv = i              + j*fv_x + k*fv_x*cv_y + var*fv_x*cv_y*cv_z;
                    fc = i-ngh + j*fc_x + k*fc_x*cv_y + var*fc_x*cv_y*cv_z;
                    affected = max(troubles[cv],troubles[cv+1]);
                    Fp = F_fv_x[fv+1]*(1-affected) + F_x[fc+1]*(affected);
                    affected = max(troubles[cv],troubles[cv-1]); 
                    Fm = F_fv_x[fv  ]*(1-affected) + F_x[fc  ]*(affected);
                    F += (Fp-Fm)/dx;
                    #endif
                    #ifdef Y
                    dy  = (Y_faces[j+1]-Y_faces[j]);
                    fv = i +  j*cv_x      + k*cv_x*fv_y + var*cv_x*fv_y*cv_z;
                    fc = i + (j-ngh)*cv_x + k*cv_x*fc_y + var*cv_x*fc_y*cv_z;
                    affected = max(troubles[cv],troubles[cv+cv_x]); 
                    Fp = F_fv_y[fv+cv_x]*(1-affected) + F_y[fc+cv_x]*(affected);
                    affected = max(troubles[cv],troubles[cv-cv_x]); 
                    Fm = F_fv_y[fv     ]*(1-affected) + F_y[fc     ]*(affected);
                    F += (Fp-Fm)/dy;
                    #endif
                    #ifdef Z
                    dz  = (Z_faces[k+1]-Z_faces[k]);
                    fv = i + j*cv_x +  k*cv_x*cv_y      + var*cv_x*cv_y*fv_z;
                    fc = i + j*cv_x + (k-ngh)*cv_x*cv_y + var*cv_x*cv_y*fc_z;
                    affected = max(troubles[cv],troubles[cv+cv_x*cv_y]); 
                    Fp = F_fv_z[fv+cv_x*cv_y]*(1-affected) + F_z[fc+cv_x*cv_y]*(affected);
                    affected = max(troubles[cv],troubles[cv-cv_x*cv_y]); 
                    Fm = F_fv_z[fv          ]*(1-affected) + F_z[fc          ]*(affected);
                    F += (Fp-Fm)/dz;
                    #endif
                    U_new[cv_U] =  U_cv[cv_U] - w_t[i_ader]*dt*F;
                }
            }
        }
    }
}
   
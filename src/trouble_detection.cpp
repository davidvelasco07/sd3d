#include "sd3d.hpp"

double max3(double a, double b, double c){
    return max(max(a,b),c);
}
double min3(double a, double b, double c){
    return min(min(a,b),c);
}

void extrema(double *U_new, double *U, double *extrema, int var, double tolerance){
    int cv;
    double maximum, minimum;
    int nm=0,np=0,mm=0,mp=0,lm=0,lp=0;
    int i=0,j=0,k=0;
    #ifdef Z
    nm=-1;np=1;
    for(k=1;k<cv_z-1;k++){
    #endif
        #ifdef Y
        mm=-1;mp=1;
        for(j=1;j<cv_y-1;j++){  
        #endif
            #ifdef X
            lm=-1;lp=1;
            for(i=1;i<cv_x-1;i++){  
            #endif
                cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                maximum = U[cv];
                minimum = U[cv];
                //for(int n=nm;n<=np;n++){
                //    for(int m=mm;m<=mp;m++){
                //        for(int l=lm;l<=lp;l++){
                //            maximum = max(maximum,U[cv +l + m*cv_x + n*cv_x*cv_y]);
                //            minimum = min(minimum,U[cv +l + m*cv_x + n*cv_x*cv_y]);
                //        }
                //    }
                //}
                #ifdef X
                maximum = max3(U[cv-1],maximum,U[cv+1]);
                minimum = min3(U[cv-1],minimum,U[cv+1]);
                #endif
                #ifdef Y
                maximum = max3(U[cv-cv_x],maximum,U[cv+cv_x]);
                minimum = min3(U[cv-cv_x],minimum,U[cv+cv_x]);
                #endif
                #ifdef Z
                maximum = max3(U[cv-cv_x*cv_y],maximum,U[cv+cv_x*cv_y]);
                minimum = min3(U[cv-cv_x*cv_y],minimum,U[cv+cv_x*cv_y]);
                #endif
                minimum -= abs(minimum)*tolerance;
                maximum += abs(maximum)*tolerance;
                if( U_new[cv] > maximum || U_new[cv] < minimum)
                    extrema[i + j*cv_x + k*cv_x*cv_y] = 1;
            #ifdef X
            }
            #endif
        #ifdef Y
        }
        #endif
    #ifdef Z
    }
    #endif
}

double Alpha(double dv, double dUm, double dU, double dUp){
    double vL,vR;
    double alpha_m,alpha_p,alphaL,alphaR;
    //vL = dU(i-1)-dU(i)
    vL = dU-dUm;
    //alphaL = min(1,-max(vL,0)/dv),1,min(1,-min(vL,0)/dv) for dv<0,dv=0,dv>0
    alpha_p = min(1.0,max(vL,0.0)/dv);
    alpha_m = min(1.0,min(vL,0.0)/dv);
    alphaL = (dv  <  0.0 ? alpha_m : alpha_p);
    alphaL = (dv  == 0.0 ? 1.0 : alphaL);
    //vR = dU(i+1)-dU(i)
    vR = dUp-dU;
    //alphaR = min(1,max(vR,0)/dv),1,min(1,min(vR,0)/dv) for dv>0,dv=0,dv<0
    alpha_p = min(1.0,max(vR,0.0)/dv);
    alpha_m = min(1.0,min(vR,0.0)/dv);
    alphaR = (dv  <  0.0 ? alpha_m : alpha_p);
    alphaR = (dv  == 0.0 ? 1.0 : alphaR);
    return min(alphaL,alphaR);
}

void smooth_extrema_x(double *U, double *alpha, int var){
    int cv;
    double dU,dUm,dUp,d2U,dv;
    for(int k=0;k<cv_z;k++){
        for(int j=0;j<cv_y;j++){
            for(int i=2;i<cv_x-2;i++){
                cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                //First derivative dUdx[i] = (U[i+1]-U[i-1])/(x[i+1]-x[i-1])
                dU  = (U[cv+1]-U[cv-1])/(X_centers[i+1]-X_centers[i-1]);
                dUp = (U[cv+2]-U[cv  ])/(X_centers[i+2]-X_centers[i  ]);
                dUm = (U[cv  ]-U[cv-2])/(X_centers[i  ]-X_centers[i-2]);
                //Second derivative d2Udx2[i] = (dU[i+1]-dU[i-1])/(x[i+1]-x[i-1])
                d2U = (dUp-dUm)/(X_centers[i+1]-X_centers[i-1]);

                dv = 0.5*d2U*(X_faces[i+1]-X_faces[i]);
                cv = i + j*cv_x + k*cv_x*cv_y;
                alpha[cv] = Alpha(dv,dUm,dU,dUp);
            }
        }
    }
}

void smooth_extrema_y(double *U, double *alpha, int var){
     int cv;
    double dU,dUm,dUp,d2U,dv;
    for(int k=0;k<cv_z;k++){
        for(int j=2;j<cv_y-2;j++){
            for(int i=0;i<cv_x;i++){
                cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                //First derivative dUdy[j] = (U[j+1]-U[j-1])/(y[j+1]-y[j-1])
                dU  = (U[cv+cv_x  ]-U[cv-cv_x  ])/(Y_centers[j+1]-Y_centers[j-1]);
                dUp = (U[cv+2*cv_x]-U[cv       ])/(Y_centers[j+2]-Y_centers[j  ]);
                dUm = (U[cv       ]-U[cv-2*cv_x])/(Y_centers[j  ]-Y_centers[j-2]);
                //Second derivative d2Udy2[j] = (dU[j+1]-dU[j-1])/(y[j+1]-y[j-1])
                d2U = (dUp-dUm)/(Y_centers[j+1]-Y_centers[j-1]);

                dv = 0.5*d2U*(Y_faces[j+1]-Y_faces[j]);
                cv = i + j*cv_x + k*cv_x*cv_y;
                alpha[cv] = Alpha(dv,dUm,dU,dUp);
            }
        }
    }
}

void smooth_extrema_z(double *U, double *alpha, int var){
    int cv,stride=cv_x*cv_y;
    double dU,dUm,dUp,d2U,dv;
    for(int k=2;k<cv_z-2;k++){
        for(int j=0;j<cv_y;j++){
            for(int i=0;i<cv_x;i++){
                cv = i + j*cv_x + k*cv_x*cv_y + var*cv_x*cv_y*cv_z;
                //First derivative dUdy[j] = (U[j+1]-U[j-1])/(y[j+1]-y[j-1])
                dU  = (U[cv+stride  ]-U[cv-stride  ])/(Z_centers[k+1]-Z_centers[k-1]);
                dUp = (U[cv+2*stride]-U[cv         ])/(Z_centers[k+2]-Z_centers[k  ]);
                dUm = (U[cv         ]-U[cv-2*stride])/(Z_centers[k  ]-Z_centers[k-2]);
                //Second derivative d2Udy2[j] = (dU[j+1]-dU[j-1])/(y[j+1]-y[j-1])
                d2U = (dUp-dUm)/(Z_centers[k+1]-Z_centers[k-1]);

                dv = 0.5*d2U*(Z_faces[k+1]-Z_faces[k]);
                cv = i + j*cv_x + k*cv_x*cv_y;
                alpha[cv] = Alpha(dv,dUm,dU,dUp);
            }
        }
    }
}

void apply_alpha(double *extrema, double *possible, double *alpha_x, double *alpha_y, double *alpha_z){
    int cv;
    double alpha;
    #ifdef X
    int i_min=NGHx*nx_cv-1;
    int i_max=(NGHx+Nx)*nx_cv+1;
    #else
    int i_min=0;
    int i_max=1;
    #endif
    #ifdef Y
    int j_min=NGHy*ny_cv-1;
    int j_max=(NGHy+Ny)*ny_cv+1;
    #else
    int j_min=0;
    int j_max=1;
    #endif
    #ifdef Z
    int k_min=NGHz*nz_cv-1;
    int k_max=(NGHz+Nz)*nz_cv+1;
    #else
    int k_min=0;
    int k_max=1;
    #endif

    for(int k=k_min;k<k_max;k++){
        for(int j=j_min;j<j_max;j++){
            for(int i=i_min;i<i_max;i++){
                alpha=0.0;
                cv = i + j*cv_x + k*cv_x*cv_y;
                #ifdef X
                alpha = max(alpha,min3(alpha_x[cv-1],alpha_x[cv],alpha_x[cv+1]));
                #endif
                #ifdef Y
                alpha = max(alpha,min3(alpha_y[cv-cv_x],alpha_y[cv],alpha_y[cv+cv_x]));
                #endif
                #ifdef Z
                alpha = max(alpha,min3(alpha_z[cv-cv_x*cv_y],alpha_z[cv],alpha_z[cv+cv_x*cv_y]));
                #endif
                //alpha==1 -> smooth extrema
                //It's the max(...) so that it throws a trouble if any tested variable exhibits a trouble
                extrema[cv] = max(extrema[cv],possible[cv]*(alpha<1)); 
            }
        }
    }
}

void PAD_criteria(double *W, double *troubles, int var, double vmin, double vmax){
    int cv;
    double value;
    for(int k=NGHz*nz_cv;k<(NGHz+Nz)*nz_cv;k++){
        for(int j=NGHy*ny_cv;j<(NGHy+Ny)*ny_cv;j++){
            for(int i=NGHx*nx_cv;i<(NGHx+Nx)*nx_cv;i++){
                cv = i + j*cv_x + k*cv_x*cv_y;
                value = W[cv+var*cv_x*cv_y*cv_z];
                if(value<vmin || value>vmax)
                    troubles[cv] = 1;
            }
        }
    }
}

void detect_field_troubles(int var){
    memset(possible_troubles, 0, size_cv*sizeof(double));
    //First we check the DMP criteria
    extrema(W_new, W_cv, possible_troubles, var, 1E-14);
    //Then relax the DMP criteria for smooth extrema
    #ifdef X
    smooth_extrema_x(W_new,dUdx,var);
    #endif
    #ifdef Y
    smooth_extrema_y(W_new,dUdy,var);
    #endif
    #ifdef Z
    smooth_extrema_z(W_new,dUdz,var);
    #endif
    apply_alpha(troubles, possible_troubles, dUdx, dUdy, dUdz);
}

void detect_troubles(){
    memset(troubles, 0, size_cv*sizeof(double));
    Boundary_Conditions(U_new,1);
    cons_to_prim(U_new,W_new,size_cv);
    cons_to_prim(U_cv ,W_cv ,size_cv);
    detect_field_troubles(0);
    #ifndef ADVECTION
    detect_field_troubles(_p_);
    #endif
    //Finally we check the PAD criteria
    PAD_criteria(W_new, troubles, 0, rho_min, rho_max);
    #ifndef ADVECTION
    PAD_criteria(W_new, troubles, _p_, p_min, p_max);
    #endif
}
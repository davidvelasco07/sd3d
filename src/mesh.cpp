#include "sd3d.hpp"

void build_faces(double *faces, double *x_fp, int N, int n, int x_i, double boxlen, double dx){
    for(int j=0;j<N;j++){
        for(int i=0;i<n+1;i++){
            faces[i+j*(n+1)]= ((j-NGH + x_i) + x_fp[i])*dx;
        }
    }
    faces[N*(n+1)] = boxlen+(NGH*dx);
}

void build_sp(double *sp, double *x_sp, int N, int n, int x_i, double boxlen, double dx){
    for(int j=0;j<N;j++){
        for(int i=0;i<n_cv;i++){
            sp[i+j*n_cv]= ((j-NGH + x_i) + x_sp[i])*dx;
        }
    }
}

void build_centers(double *centers, double *faces, int N){
    for(int i=0;i<N;i++){
        centers[i] = 0.5*(faces[i]+faces[i+1]);
    }
}

void Build_mesh(){
    int i,j,k;    
    cells =  cells_x * cells_y * cells_z;
    size_cv = cv_x * cv_y * cv_z;
    total_size_cv = size_cv*nvar;
    //Build Mesh
    #ifdef X
    X_faces = malloc_host<double>(cv_faces_x);
    X_centers = malloc_host<double>(cv_x);
    X_sp = malloc_host<double>(cv_x);
    build_faces(X_faces, x_fp, cells_x, n, x_i, boxlen_x/cpu_x, dx);
    build_centers(X_centers, X_faces, cv_x);
    build_sp(X_sp, x_sp, cells_x, n, x_i, boxlen_x/cpu_x, dx);
    #endif
    #ifdef Y
    Y_faces = malloc_host<double>(cv_faces_y);
    Y_centers = malloc_host<double>(cv_y);
    Y_sp = malloc_host<double>(cv_y);
    build_faces(Y_faces, x_fp, cells_y, n, y_i, boxlen_y/cpu_y, dy);
    build_centers(Y_centers, Y_faces, cv_y);
    build_sp(Y_sp, x_sp, cells_y, n, y_i, boxlen_y/cpu_y, dy);
    #endif
    #ifdef Z
    Z_faces = malloc_host<double>(cv_faces_z);
    Z_centers = malloc_host<double>(cv_z);
    build_faces(Z_faces, x_fp, cells_z, n, z_i, boxlen_z/cpu_z, dz);
    build_centers(Z_centers, Z_faces, cv_z);
    #endif

    U_cv  = malloc_host<double>(total_size_cv);
    U_new = malloc_host<double>(total_size_cv);
    W_cv  = malloc_host<double>(total_size_cv);
    W_new = malloc_host<double>(total_size_cv);

    U_sp  = malloc_host<double>(total_size_cv);
   
    //We might only need to store space for the active zones in this case
    U_ader_sp = malloc_host<double>(total_size_cv*(n+1));
    dU_ader_sp = malloc_host<double>(total_size_cv*(n+1));
    #ifdef X
    U_ader_fp_x = malloc_host<double>(n_fp*cells_x*cv_y*cv_z*n_cv*nvar);
    W_ader_fp_x = malloc_host<double>(n_fp*cells_x*cv_y*cv_z*n_cv*nvar);
    F_ader_fp_x = malloc_host<double>(n_fp*cells_x*cv_y*cv_z*n_cv*nvar);
    #endif
    #ifdef Y
    U_ader_fp_y = malloc_host<double>(cv_x*n_fp*cells_y*cv_z*n_cv*nvar);
    W_ader_fp_y = malloc_host<double>(cv_x*n_fp*cells_y*cv_z*n_cv*nvar);
    F_ader_fp_y = malloc_host<double>(cv_x*n_fp*cells_y*cv_z*n_cv*nvar);
    #endif
    #ifdef Z
    U_ader_fp_z = malloc_host<double>(cv_x*cv_y*n_fp*cells_z*n_cv*nvar);
    W_ader_fp_z = malloc_host<double>(cv_x*cv_y*n_fp*cells_z*n_cv*nvar);
    F_ader_fp_z = malloc_host<double>(cv_x*cv_y*n_fp*cells_z*n_cv*nvar);
    #endif
    int flux_size = max(max((cv_faces_x)*cv_y*cv_z,cv_x*(cv_faces_y)*cv_z),cv_x*cv_y*(cv_faces_z));
    U_L = malloc_host<double>(flux_size*nvar); 
    U_R = malloc_host<double>(flux_size*nvar);
    #ifdef X
    F_x = malloc_host<double>(flux_size*nvar); 
    #endif
    #ifdef Y
    F_y = malloc_host<double>(flux_size*nvar); 
    #endif
    #ifdef Z
    F_z = malloc_host<double>(flux_size*nvar); 
    #endif

    #ifdef X
    BC_x_size = NGH*n_cv*cv_y*cv_z*nvar;
    BC_x[0] = malloc_host<double>(BC_x_size); 
    BC_x[1] = malloc_host<double>(BC_x_size); 
    #endif
    #ifdef Y
    BC_y_size = cv_x*NGH*n_cv*cv_z*nvar;
    BC_y[0] = malloc_host<double>(BC_y_size);
    BC_y[1] = malloc_host<double>(BC_y_size);
    #endif
    #ifdef Z
    BC_z_size = cv_x*cv_y*NGH*n_cv*nvar;
    BC_z[0] = malloc_host<double>(BC_z_size);
    BC_z[1] = malloc_host<double>(BC_z_size);
    #endif

    #ifndef SD
    #ifdef X
    F_fv_x = malloc_host<double>(n_fp*cells_x*cv_y*cv_z*nvar);
    #endif
    #ifdef Y
    F_fv_y = malloc_host<double>(cv_x*n_fp*cells_y*cv_z*nvar);
    #endif
    #ifdef Z
    F_fv_z = malloc_host<double>(cv_x*cv_y*n_fp*cells_z*nvar);
    #endif
    possible_troubles = malloc_host<double>(cv_x*cv_y*cv_z);
    troubles = malloc_host<double>(cv_x*cv_y*cv_z);
    //FallBack Scheme
    dUdx = malloc_host<double>(cv_x*cv_y*cv_z*nvar);
    dUdy = malloc_host<double>(cv_x*cv_y*cv_z*nvar);
    dUdz = malloc_host<double>(cv_x*cv_y*cv_z*nvar);
    #endif
    Write_edges();
}

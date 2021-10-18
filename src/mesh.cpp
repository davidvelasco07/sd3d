#include "sd3d.hpp"

void build_faces(double *faces, double *x_fp, int N, int x_i, double dx){
    for(int j=0;j<N;j++){
        for(int i=0;i<n_cv;i++){
            faces[i+j*n_cv]= ((j-NGH + x_i) + x_fp[i])*dx;
        }
    }
    faces[N*n_cv] = (N-NGH + x_i)*dx;
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
    build_faces(X_faces, x_fp, cells_x, x_i, dx);
    build_centers(X_centers, X_faces, cv_x);
    #endif
    #ifdef Y
    Y_faces = malloc_host<double>(cv_faces_y);
    Y_centers = malloc_host<double>(cv_y);
    build_faces(Y_faces, x_fp, cells_y, y_i, dy);
    build_centers(Y_centers, Y_faces, cv_y);
    #endif
    #ifdef Z
    Z_faces = malloc_host<double>(cv_faces_z);
    Z_centers = malloc_host<double>(cv_z);
    build_faces(Z_faces, x_fp, cells_z, z_i, dz);
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
    int flux_size = max(max((cv_faces_x-1)*cv_y*cv_z,cv_x*(cv_faces_y-1)*cv_z),cv_x*cv_y*(cv_faces_z-1));
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
    BC_x[0] = malloc_host<double>(BC_x_size*n_cv); 
    BC_x[1] = malloc_host<double>(BC_x_size*n_cv); 
    #endif
    #ifdef Y
    BC_y_size = cv_x*NGH*n_cv*cv_z*nvar;
    BC_y[0] = malloc_host<double>(BC_y_size*n_cv);
    BC_y[1] = malloc_host<double>(BC_y_size*n_cv);
    #endif
    #ifdef Z
    BC_z_size = cv_x*cv_y*NGH*n_cv*nvar;
    BC_z[0] = malloc_host<double>(BC_z_size*n_cv);
    BC_z[1] = malloc_host<double>(BC_z_size*n_cv);
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

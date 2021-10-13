#include "sd3d.hpp"

void Init_variables(){
    int i,j,k;
    n_cv = n+1;
    n_fp = n+2;
    nvar = 1; //Density
    #ifdef X
    nx = n;
    cells_x = (Nx+2*NGH);
    faces_x = cells_x+1;
    cv_x = cells_x*(n+1);
    cv_faces_x = cv_x+1;
    _vx_ = nvar;
    nvar++;
    #else
    nx = 1;
    cells_x = 1;
    faces_x = 1;
    cv_x = 1;
    cv_faces_x = 1;
    #endif
    #ifdef Y
    ny = n;
    cells_y = (Ny+2*NGH);
    faces_y = cells_y+1;
    cv_y = cells_y*(n+1);
    cv_faces_y = cv_y+1;
    _vy_ = nvar;
    nvar++;
    #else
    ny = 1;
    cells_y = 1;
    faces_y = 1;
    cv_y = 1;
    cv_faces_y=1;
    #endif
    #ifdef Z
    nz = n;
    cells_z = (Nz+2*NGH);
    faces_z = cells_z+1;
    cv_z = cells_z*(n+1);
    cv_faces_z = cv_z+1;
    _vz_ = nvar;
    nvar++;
    #else
    nz = 1;
    cells_z = 1;
    faces_z = 1;
    cv_z = 1;
    cv_faces_z = 1;
    #endif
    _e_ = _p_ = nvar;
    nvar++; //Energy

    dx=boxlen_x/NX;
    dy=boxlen_y/NY;
    dz=boxlen_z/NZ;

    if(Master)
        cout<< "n = "<<n<<", N = ("<<NX<<","<<NY<<","<<NZ<<")"<<endl;

    double *x = malloc_host<double>(n);
    double *w = malloc_host<double>(n);

    x_t = malloc_host<double>(n_cv);
    w_t = malloc_host<double>(n_cv);

    gauss_legendre(0.0, 1.0, n, x, w);
    gauss_legendre(0.0, 1.0, n_cv, x_t, w_t);

    x_sp = malloc_host<double>(n_cv);   
    x_fp = malloc_host<double>(n_fp);

    flux_points(x_fp,x,n);
    solution_points(x_sp,n);
    
    if(Master){
        cout<<"x_fp = [";
        for(i=0;i<n_cv;i++){
            cout<<x_fp[i]<<",";
        }
        cout<<x_fp[n_cv]<<"]"<<endl; 

        cout<<"x_sp = [";
        for(i=0;i<n;i++){
            cout<<x_sp[i]<<",";
        }
        cout<<x_sp[n]<<"]"<<endl;

        cout<<"w_t = [";
        for(i=0;i<n;i++){
            cout<<w_t[i]<<",";
        }
        cout<<w_t[n]<<"]"<<std::endl;
    }
    sp_to_fp  = malloc_host<double>(n_cv*n_fp);
    fp_to_sp  = malloc_host<double>(n_fp*n_cv);
    dfp_to_sp = malloc_host<double>(n_fp*n_cv);
    sp_to_cv  = malloc_host<double>(n_cv*n_cv);
    cv_to_sp  = malloc_host<double>(n_cv*n_cv);
    ader      = malloc_host<double>(n_cv*n_cv);
    invader   = malloc_host<double>(n_cv*n_cv);

    lagrange_matrix(sp_to_fp, x_sp, x_fp, n_cv, n_fp);
    lagrange_matrix(fp_to_sp, x_fp, x_sp, n_fp, n_cv);
    lagrange_prime_matrix(dfp_to_sp, x_fp, x_sp, n_fp, n_cv);
    integral_matrix(sp_to_cv, x_fp, x_sp, x, w, n);
    inverse(sp_to_cv, cv_to_sp, n_cv);
    ader_matrix(ader,x_t,w_t,n_cv);
    inverse(ader,invader,n_cv);
    
}
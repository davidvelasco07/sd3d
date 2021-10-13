
//MPI global variables
extern int cpu_rank;
extern int cpu_size;
extern int Comm;
extern int cpu_x;
extern int cpu_y;
extern int cpu_z;
extern int rank_x;
extern int rank_y;
extern int rank_z;
extern int x_i;
extern int y_i;
extern int z_i;

extern int Master;

extern int N_comms;

//Global variables
extern int _x_;
extern int _y_;
extern int _z_;
extern int _vx_;
extern int _vy_;
extern int _vz_;
extern int _e_;
extern int _p_;

extern int N;
extern int NX;
extern int NY;
extern int NZ;
extern int Nx;
extern int Ny;
extern int Nz;
extern int n;
extern int nx;
extern int ny;
extern int nz;
extern int nvar;
extern double boxlen_x;
extern double boxlen_y;
extern double boxlen_z;
extern double dx;
extern double dy;
extern double dz;
extern double cfl;
extern double dt;
extern double t;
extern double t_end;
extern double min_c2;

extern double gmma;

extern int n_step;
extern int n_output;

extern int n_fp;
extern int n_cv;

extern int cells_x;
extern int cells_y;
extern int cells_z;
extern int cells;

extern int size_cv;
extern int total_size_cv;

extern int faces_x;
extern int faces_y;
extern int faces_z;

extern int cv_x;
extern int cv_y;
extern int cv_z;
extern int cv_faces_x;
extern int cv_faces_y;
extern int cv_faces_z;

extern double *X_centers;
extern double *Y_centers;
extern double *Z_centers;
extern double *X_faces;
extern double *Y_faces;
extern double *Z_faces;
extern double *X_sp;
extern double *Y_sp;
extern double *Z_sp;

//Matrices for interchangeability
extern double *x_sp;
extern double *x_fp;
extern double *x_t;
extern double *w_t;
extern double *sp_to_fp;
extern double *fp_to_sp;
extern double *dfp_to_sp;
extern double *sp_to_cv;
extern double *cv_to_sp;
extern double *ader;
extern double *invader;

extern double *U_cv;
extern double *U_new;
extern double *W_cv; 
extern double *W_new;
extern double *U_sp;

//ADER fields
extern double *U_ader_sp;
extern double *dU_ader_sp;
extern double *U_ader_fp_x;
extern double *U_ader_fp_y;
extern double *U_ader_fp_z;
extern double *W_ader_fp_x;
extern double *W_ader_fp_y;
extern double *W_ader_fp_z;
extern double *F_ader_fp_x;
extern double *F_ader_fp_y;
extern double *F_ader_fp_z;

//FV fields
extern double *F_fv_x;
extern double *F_fv_y;
extern double *F_fv_z;

extern double *U_R;
extern double *U_L;
extern double *F_x;
extern double *F_y;
extern double *F_z;

extern double *BC_x[2];
extern double *BC_y[2];
extern double *BC_z[2];
extern int BC_x_size;
extern int BC_y_size;
extern int BC_z_size;

extern int cpu_neighbour[3][2];
extern int cpu_boundary[3][2];
extern comm_cpu comms[3][2];

//Trouble
extern double *possible_troubles;
extern double *troubles;

//FallBack Scheme
extern double *dUdx;
extern double *dUdy;
extern double *dUdz;

extern double rho_min;
extern double rho_max;
extern double p_min;
extern double p_max;

extern void (*riemann_solver)(double*, double*, int, int, int, int);
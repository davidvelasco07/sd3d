//#include "sd3d.hpp"
//MPI global variables
int cpu_rank;
int cpu_size;
int Comm;
int cpu_x=1;
int cpu_y=1;
int cpu_z=1;
int rank_x;
int rank_y;
int rank_z;
int x_i=0;
int y_i=0;
int z_i=0;

int Master=1;

int N_comms;

int _x_=0;
int _y_=1;
int _z_=2;
int _vx_;
int _vy_;
int _vz_;
int _e_;
int _p_;

int N;
int NX;
int NY;
int NZ;
int Nx;
int Ny;
int Nz;
int n;
int nx;
int ny;
int nz;
int nvar;
double boxlen_x;
double boxlen_y;
double boxlen_z;
double dx;
double dy;
double dz;
double cfl;
double dt;
double t;
double t_end;
double gmma;
double min_c2=1E-14;
int n_step;
int n_output;

int n_fp;
int n_cv;

int cells_x;
int cells_y;
int cells_z;
int cells;

int size_cv;
int total_size_cv;

int faces_x;
int faces_y;
int faces_z;

int cv_x;
int cv_y;
int cv_z;
int cv_faces_x;
int cv_faces_y;int cv_faces_z;

double *X_centers;
double *Y_centers;
double *Z_centers;
double *X_faces;
double *Y_faces;
double *Z_faces;
double *X_sp;
double *Y_sp;
double *Z_sp;

double *x_sp;
double *x_fp;
double *x_t;
double *w_t;
double *sp_to_fp;
double *fp_to_sp;
double *dfp_to_sp;
double *sp_to_cv;
double *cv_to_sp;
double *ader;
double *invader;

double *U_cv;
double *U_new;
double *W_cv; 
double *W_new;
double *U_sp;

//ADER fields
double *U_ader_sp;
double *dU_ader_sp;
double *U_ader_fp_x;
double *U_ader_fp_y;
double *U_ader_fp_z;
double *W_ader_fp_x;
double *W_ader_fp_y;
double *W_ader_fp_z;
double *F_ader_fp_x;
double *F_ader_fp_y;
double *F_ader_fp_z;

//FV fields
double *F_fv_x;
double *F_fv_y;
double *F_fv_z;

double *U_R;
double *U_L;
double *F_x;
double *F_y;
double *F_z;

double *BC_x[2];
double *BC_y[2];
double *BC_z[2];
int BC_x_size;
int BC_y_size;
int BC_z_size;

int cpu_neighbour[3][2];
int cpu_boundary[3][2];

//Trouble
double *possible_troubles;
double *troubles;

//FallBack Scheme
double *dUdx;
double *dUdy;
double *dUdz;

double rho_min=1E-5;
double rho_max=1E5;
double p_min=1E-5;
double p_max=1E5;

void (*riemann_solver)(double*, double*, int, int, int, int);
#define MPI

#define PI 3.141592653589793

//#define SD
#define NGH 1
#define X
#define Y
//#define Z

//#define ADVECTION

#define Periodic_X
#define Periodic_Y
#define Periodic_Z

#if (defined(X) && defined(Y)) || (defined(X) && defined(Z)) || (defined(Z) && defined(Y))
#define _2D_
#endif
#if defined(X) && defined(Y) && defined(Z)
#define _3D_
#endif
#define LLF

#ifdef X
#define nx_cv (nx+1)
#define nx_fp (nx+2)
#define NGHx NGH
#else
#define nx_cv 1
#define nx_fp 1
#define NGHx 0
#endif
#ifdef Y
#define ny_cv (ny+1)
#define ny_fp (ny+2)
#define NGHy NGH
#else
#define ny_cv 1
#define ny_fp 1
#define NGHy 0
#endif
#ifdef Z
#define nz_cv (nz+1)
#define nz_fp (nz+2)
#define NGHz NGH
#else
#define nz_cv 1
#define nz_fp 1
#define NGHz 0
#endif
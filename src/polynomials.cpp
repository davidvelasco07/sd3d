#include "sd3d.hpp"

double lagrange(double *x, double y, int i, int n){
    int j;
    double lag=1.;
    for(j=0;j<n;j++){
        if(j!=i)
            lag *= (y-x[j])/(x[i]-x[j]);
    }
    return lag;
}

void lagrange_matrix(double *sp_to_fp, double *x_sp, double *x_fp, int n_sp, int n_fp){
  for(int j=0;j<n_fp;j++){
        for(int i=0;i<n_sp;i++){
            sp_to_fp[i+j*n_sp]=lagrange(x_sp,x_fp[j],i,n_sp);
        }
    }
}

double lagrange_prime(double *x, double y, int i, int n){
    int j,k;
    double lag,lagp=0;
    for(k=0;k<n;k++){
      if(k!=i){
        lag=1;
        for(j=0;j<n;j++){
            if(j!=i && j!=k){
              lag *= (y-x[j])/(x[i]-x[j]);
            }
        }
        lagp+=lag/(x[i]-x[k]);
      }
    }
    return lagp;
}

void lagrange_prime_matrix(double *fp_to_sp, double *x_fp, double *x_sp, int n_fp, int n_sp){
  for(int j=0;j<n_sp;j++){
        for(int i=0;i<n_fp;i++){
            fp_to_sp[i+j*n_fp]=lagrange_prime(x_fp,x_sp[j],i,n_fp);
        }
    }
}

void gauss_legendre(double xi, double xf, int n, double *x, double *w){
    int i,j,m;
    double p1, p2, p3, pp, xl, xm, z, z1;
    double eps=3E-14;
    m = (n+1)/2;
    xm = 0.5*(xf+xi);
    xl = 0.5*(xf-xi);
    if(n>0){
      for(i=0;i<m;i++){
          z = cos(3.141592654*(i+0.75)/(n+0.5));
          z1 = 0.0;
          while(abs(z-z1) > eps){
              p1 = 1.0;
              p2 = 0.0;
              for(j=0;j<n;j++){
                  p3 = p2;
                  p2 = p1;
                  p1 = ((2.0*j+1)*z*p2-j*p3)/(j+1);
              }
              pp = n*(z*p1-p2)/(z*z-1.0);
              z1 = z;
              z = z1 - p1/pp;   
          }
          x[i] = xm - xl*z;
          x[n-i-1] = xm + xl*z;
          w[i] = (2.0*xl)/((1.0-z*z)*pp*pp);
          w[n-i-1] = w[i];
      }
    }
    else{
      x[0]=0.5;
      w[0]=1;
    }
}

void flux_points(double *x_fp, double *x, int n){
  x_fp[0]=0.0;
  for(int i=1;i<n+1;i++){
    x_fp[i] = x[i-1];
  }
  x_fp[n+1]=1.0;
}

void solution_points(double *x_sp, int n){
  for(int i=0;i<=n;i++){
    x_sp[i] = 0.5*(1.0 - cos((2*i + 1)/(2.0*(n+1))*PI));   
  }
}

void ader_matrix(double *ader, double *x_t, double *w_t, int n_cv){
    for(int j=0;j<n_cv;j++){
        for(int i=0;i<n_cv;i++){
          ader[i+j*n_cv]=lagrange(x_t,1,i,n_cv)*lagrange(x_t,1,j,n_cv)-lagrange_prime(x_t,x_t[i],j,n_cv)*w_t[i];
        }
    }
}

void integral_matrix(double *sp_to_cv, double *x_fp, double *x_sp, double *x, double *w, int n){
    double integral;
    int n_cv=n+1;
    for(int k=0;k<n_cv;k++){
        if(n>0)
            gauss_legendre(x_fp[k], x_fp[k+1], n, x, w);
        for(int j=0;j<n_cv;j++){
            if(n>0){
                integral=0.0;
                for(int i=0;i<n;i++){
                    integral+=lagrange(x_sp,x[i],j,n+1)*w[i];
                }
            }
            else
                integral = 1.0;
            sp_to_cv[j+k*n_cv]=integral/(x_fp[k+1]-x_fp[k]);
        }
    }
    gauss_legendre(0.0, 1.0, n, x, w);
}

void inverse(double *A, double *C, int n){
  double coeff;
  double *b = malloc_host<double>(n);
  double *d = malloc_host<double>(n);
  double *x = malloc_host<double>(n);

  double *B = malloc_host<double>(n*n);
  double *L = malloc_host<double>(n*n);
  double *U = malloc_host<double>(n*n);
  int i,j,k;
  for(j=0; j<n; j++){
    b[j]=0.0;
    for(i=0; i<n; i++){
      B[i+j*n]=A[i+j*n];
      U[i+j*n]=0;
      L[i+j*n]=0;
    }
  }
  //Step 1: forward elimination
  for(k=0; k<(n-1); k++){
    for(i=k+1; i<n; i++){
      coeff = B[k+i*n]/B[k+k*n];
      L[k+i*n] = coeff;
      for(j=k+1; j<n; j++){
        B[j+i*n] -= coeff*B[j+k*n];
      }
    }
  }
  //Step 2: prepare L and U matrices 
  //L matrix is a matrix of the elimination coefficient
  // + the diagonal elements are 1.0
  for(i=0;i<n;i++){
    L[i+i*n]=1.0;
  }
  //U matrix is the upper triangular part of A
  for(j=0; j<n; j++){
    for(i=0; i<=j; i++){
      U[j+i*n] = B[j+i*n];
    }
  }
  //Step 3: compute columns of the inverse matrix C
  for(k=0;k<n;k++){
    b[k]=1.0;
    d[0]=b[0];
    //Step 3a: Solve Ld=b using the forward substitution
    for(i=1;i<n;i++){
      d[i]=b[i];
      for(j=0;j<i;j++){
        d[i] -= L[j+i*n]*d[j];
      }
    }
    //Step 3b: Solve Ux=d using the back substitution
    x[n-1]=d[n-1]/U[n*n-1];
    for(i=n-2;i>=0;i--){
      x[i]=d[i];
      for(j=n-1;j>=i+1;j--){
        x[i]=x[i]-U[j+i*n]*x[j];
      }
      x[i] = x[i]/U[i+i*n];
    }
    //Step 3c: fill the solutions x(n) into column k of C
    for(i=0;i<n;i++){
      C[k+i*n] = x[i];
    }
    b[k]=0.0;
  }
  free(B);
  free(L);
  free(U);
  free(x);
  free(b);
  free(d);
}

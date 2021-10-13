#include "sd3d.hpp"

double square_signal(int,double,double,double);
double sine_wave(int,double,double,double);
double kh_instability(int,double,double,double);
double gaussian_plus_square_signal(int,double,double,double);
double sedov_blast(int,double,double,double);
double sod_3D(int,double,double,double);
double spherical_blast(int,double,double,double);

void Initial_Conditions(){
    double x,y,z,dx,dy,dz;
    double s,value;
    for(int var=0; var<nvar; var++){
        for(int k=0; k<cv_z; k++){
            for(int j=0; j<cv_y; j++){
                dy=Y_faces[j+1]-Y_faces[j];
                for(int i=0; i<cv_x; i++){
                    dx=X_faces[i+1]-X_faces[i];
                    value=0.0;
                    for( int m=0; m<ny_cv; m++){
                        y = Y_faces[j] + x_t[m]*dy;
                        for( int l=0; l<nx_cv; l++){
                            x = X_faces[i] + x_t[l]*dx;
                            //s=square_signal(var,x,y,z);
                            //s=sine_wave(var,x,y,z);
                            //s=kh_instability(var,x,y,z);
                            //s=sod_3D(var,x,y,z);
                            s=spherical_blast(var,x,y,z);
                            s*=w_t[m];
                            s*=w_t[l];
                            value+=s;
                        }
                    }
                    W_cv[i + j*cv_x + k*cv_x*cv_y + var*size_cv] = value;
                }
            }
        }
    }
    prim_to_cons(W_cv,U_cv,size_cv);
    transform_cv_to_sp(U_cv,U_sp);    
}

double square_signal(int var, double x, double y, double z){
    rho_min=1;
    rho_max=2;
    if(var==0){
        if( abs(x-0.5)<0.25 && abs(y-0.5)<0.25)
            return 2;
        else
            return 1;
    }
    else
        return 1;
}

double gaussian_plus_square_signal(int var, double x, double y, double z){
    double gaussian;
    if(var==0){
        gaussian = 1+exp(-0.5*pow((x-0.25)/0.05,2));
        if( abs(x-0.7)<0.1 && abs(y-0.7)<0.1)
            return gaussian+1;
        else
            return gaussian;
    }
    else
        return 0;
}

double sine_wave(int var, double x, double y, double z){
    if(var==0)
       return 1.0+0.125*(sin(2*PI*(x+y)));
    else
        return 1;
}

double kh_instability(int var, double x, double y, double z){
    double w0=0.1;
    double sigma2 = pow(0.05,2);
    double rho,vx,vy;
    if( y>0.25 && y<0.75){
        rho = 2;
        vx = 0.5;
    }
    else{
        rho = 1;
        vx = -0.5;
    }
    if(var==0)
        return rho;
    else if(var==_vx_)
        return vx;
    else if(var==_vy_)
        return w0*sin(4*PI*x)*(exp(-pow(y-0.25,2)/sigma2)+exp(-pow(y-0.75,2)/sigma2));
    else if(var==_p_)
        return 2.5;
    else
        return 0;
}

double sedov_blast(int var, double x, double y, double z){
    gmma=5./3.;
    double r,R=0.01;
    #ifndef _3D_
    r=sqrt(pow(x-0.5,2) + pow(y-0.5,2));
    #else
    r=sqrt(pow(x-0.5,2) + pow(y-0.5,2) + pow(z-0.5,2));
    #endif
    if(var==0){
        return 1.0;
    }
    else if(var==_p_){
        if(r<R)
            return 1E-6+(gmma-1);
        else
            return 1E-6;
    }
    else
        return 0;
}

double sod_3D(int var, double x, double y, double z){
    gmma=1.4;
    double r,R=0.25;
    #ifndef _3D_
    r=sqrt(pow(x-0.5,2.0) + pow(y-0.5,2.0));
    #else
    r=sqrt(pow(x-0.5,2) + pow(y-0.5,2) + pow(z-0.5,2));
    #endif
    if(var==0){
        if(r<R)
            return 1;
        else
            return 0.125;
    }
    else if(var==_p_){
        if(r<R)
            return 1;
        else
            return 0.1;
    }
    else
        return 0;
}

double spherical_blast(int var, double x, double y, double z){
    gmma=5./3.;
    double r,R=0.1;
    #ifndef _3D_
    r=sqrt(pow(x-0.5,2) + pow(y-0.5,2));
    #else
    r=sqrt(pow(x-0.5,2) + pow(y-0.75,2) + pow(z-0.5,2));
    #endif
    if(var==0){
        return 1;
    }
    else if(var==_p_){
        if(r<R)
            return 10;
        else
            return 0.1;
    }
    else
        return 0;
}
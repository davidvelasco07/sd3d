#include "sd3d.hpp"

double square_signal(int,double,double,double);
double sine_wave(int,double,double,double);
double kh_instability(int,double,double,double);
double gaussian_plus_square_signal(int,double,double,double);
double sedov_blast(int,double,double,double);
double sod_3D(int,double,double,double);
double spherical_blast(int,double,double,double);
double isentropic_vortex(int,double,double,double);
double sound_wave(int,double,double,double);

void Initial_Conditions(){
    double x,y,z,dx,dy,dz;
    double s,value;
    for(int var=0; var<nvar; var++){
        for(int k=0; k<cv_z; k++){
            #ifdef Z
            dz=Z_faces[k+1]-Z_faces[k];
            #endif
            for(int j=0; j<cv_y; j++){
                #ifdef Y
                dy=Y_faces[j+1]-Y_faces[j];
                #endif
                for(int i=0; i<cv_x; i++){
                    #ifdef X
                    dx=X_faces[i+1]-X_faces[i];
                    #endif
                    value=0.0;
                    for( int n=0; n<nz_cv; n++){
                        #ifdef Z
                        z = Z_faces[k] + x_t[n]*dz;
                        #endif
                        for( int m=0; m<ny_cv; m++){
                            #ifdef Y
                            y = Y_faces[j] + x_t[m]*dy;
                            #endif
                            for( int l=0; l<nx_cv; l++){
                                #ifdef X
                                x = X_faces[i] + x_t[l]*dx;
                                #endif
                                //s=square_signal(var,x,y,z);
                                //s=sine_wave(var,x,y,z);
                                //s=kh_instability(var,x,y,z);
                                //s=sod_3D(var,x,y,z);
                                //s=spherical_blast(var,x,y,z);
                                //s=sedov_blast(var,x,y,z);
                                //s=isentropic_vortex(var,x,y,z);
                                s=sound_wave(var,x,y,z);
                                #ifdef X
                                s*=w_t[l];
                                #endif
                                #ifdef Y
                                s*=w_t[m];
                                #endif
                                #ifdef Z
                                s*=w_t[n];
                                #endif
                                value+=s;
                            }
                        }
                    }
                    W_cv[i + j*cv_x + k*cv_x*cv_y + var*size_cv] = value;
                }
            }
        }
    }
    prim_to_cons(W_cv,U_cv,size_cv);
    transform_cv_to_sp(U_cv,U_sp); 
    for(int var=0;var<nvar;var++){
        for(int i=0;i<n_cv;i++)
            memcpy(&U_ader_sp[i*size_cv+var*size_cv*n_cv], &U_sp[var*size_cv], size_cv*sizeof(double));
    }   
}

double square_signal(int var, double x, double y, double z){
    rho_min=1;
    rho_max=2;
    if(var==0){
        #ifdef _3D_
        if( abs(x-0.5)<0.25 && abs(y-0.5)<0.25 && abs(z-0.5)<0.25 )
        #else
        if( abs(x-0.5)<0.25 && abs(y-0.5)<0.25 )
        #endif
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
    double xc=boxlen_x/2;
    double yc=boxlen_y/2;
    double zc=boxlen_z/2;
    #ifndef _3D_
    r=sqrt(pow(x-xc,2) + pow(y-yc,2));
    #else
    r=sqrt(pow(x-xc,2) + pow(y-yc,2) + pow(z-zc,2));
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

double isentropic_vortex(int var, double x, double y, double z){
    gmma=1.4;
    double r,r2;
    double xc=boxlen_x/2;
    double yc=boxlen_y/2;
    double zc=boxlen_z/2;
    double vphi,rho,beta = 5.;
    #ifndef _3D_
    r=sqrt(pow(x-xc,2) + pow(y-yc,2));
    #else
    r=sqrt(pow(x-xc,2) + pow(y-yc,2) + pow(z-zc,2));
    #endif
    r2   = pow(r,2.);
    vphi = beta/(2.*PI)*exp(0.5*(1-r2));
    rho  = pow(1.-(gmma-1.0)*pow(beta,2.)/(8*gmma*pow(PI,2.))*exp(1-r2),1./(gmma-1.));
    if(var==0){
        return rho;
    }
    else if(var==_p_){
        return pow(rho,gmma);
    }
    else if(var==_vx_){
        return -vphi*(y-yc)+0.;
    }
     else if(var==_vy_){
        return vphi*(x-xc)+0.;
    }
    else
        return 0;
}

double sound_wave(int var, double x, double y, double z){
    gmma=1.4;
    double rho=1,p=1,cs0=1.;
    double kmode  = 2.*PI/boxlen_x;
    double A = 1E-4;
    double lambda_r = 0.0, lambda_i = 2.*PI;
    double vel_r = lambda_i/kmode;
    double vel_i = lambda_r/kmode;
    double dvel = A*(sin(kmode*x));//A*(vel_r*cos(kmode*x) - vel_i*sin(kmode*x));

    rho *= (1+dvel/cs0);
    p   *= (1+gmma*dvel/cs0);

    if(var==0){
        return rho; 
    }
    else if(var==_p_){
        return p;
    }
    else if(var==_vx_){
        return dvel;
    }
     else if(var==_vy_){
        return 0;
    }
    else
        return 0;
}

//double sound_wave(int var, double x, double y, double z){
//    double epsilon = 0.0001;
//    double cs=1;
//    double projx = 1./sqrt(1.+pow(boxlen_y/boxlen_x,2.0));
//    double projy = 1./sqrt(1.+pow(boxlen_x/boxlen_y,2.0));
//    double qx = x/boxlen_x;
//    double qy = y/boxlen_y;
//    double kernel = exp(-(((qx-.5)*(qx-.5)+(qy-.5)*(qy-.5))*25.));
//    double rho = (1.+kernel*epsilon*cos(20.*PI*(qx+qy)));
//    double vx  = cs*epsilon*cos(20.*PI*(qx+qy))*kernel*projx;
//    double vy  = cs*epsilon*cos(20.*PI*(qx+qy))*kernel*projy;
//    double p   = (1.+gmma*kernel*epsilon*cos(20.*PI*(qx+qy)))*cs*cs/gmma;
//    if(var==0)
//        return rho;
//    else if(var==_vx_)
//        return vx;
//    else if(var==_vy_)
//        return vy;
//    else if(var==_p_)
//        return p;
//    else return 0;
//}
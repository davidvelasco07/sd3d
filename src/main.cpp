#include "sd3d.hpp"

void finish(){
    #ifdef MPI
    MPI_Finalize();
    #endif
    exit(0);
}

int main(int argc, char** argv){
    int i,j,k,m;
#ifdef MPI
    MPI_Init (&argc, &argv);
    Comm = MPI_COMM_WORLD;
    MPI_Comm_size (Comm, &cpu_size);
    MPI_Comm_rank (Comm, &cpu_rank);
    Master=(cpu_rank==0);
#endif
    if(Master){
        cout<< "#####################"<< std::endl;
        cout<< "### SD3D          ###"<< std::endl;
        cout<< "#####################"<< std::endl;
    #ifdef MPI
        cout<< "Parallel version with "<< cpu_size<< " processes"<<endl<<endl;
    #else
        cout<< "Serial version -- ";
    #endif
    #ifdef _3D_
        cout<<"3D"<<endl<<endl;
    #else
        #ifdef _2D_
        cout<<"2D"<<endl<<endl;
        #else
        cout<<"1D"<<endl<<endl;
        #endif
    #endif
    }
    //Parameters
    n=6;
    N = 8;
    #ifdef X
    NX = N;
    #else
    NX = Nx = 1;
    #endif
    #ifdef Y
    NY = int(N*1.);
    #else
    NY = 1;
    #endif
    #ifdef Z
    NZ = N;
    #else
    NZ = 1;
    #endif
    Nx=NX;
    Ny=NY;
    Nz=NZ;
    boxlen_x = 10.;
    boxlen_y = 10.;
    boxlen_z = 1.;
    cfl=0.6;
    gmma=7./5.;
    t=0;
    t_end=8.451542547285166;
    double dt_output=8.451542547285166;
    double t_output=dt_output;
    double dt_min;
    Split_Domain();
    Init_variables();
    Build_mesh();
    Build_comms();
    Initial_Conditions();
 
    Store_boundaries_ader(U_ader_sp);
    #ifndef SD
    Store_boundaries(U_cv);
    #endif

    n_step=0;
    n_output=0;

    compute_dt();
    Write(n_output++);

    while(t<t_end){
        t+=dt;
        n_step++;
        //Initialize ADER time slices
        for(int var=0;var<nvar;var++){
            for(i=0;i<n_cv;i++)
                memcpy(&U_ader_sp[i*size_cv+var*size_cv*n_cv], &U_sp[var*size_cv], size_cv*sizeof(double));
        }
        //Picard iteration
        for(i=0;i<n_cv;i++){   
            #ifdef X
            transform_sp_to_fp_x(U_ader_sp,U_ader_fp_x,n_cv);
            cons_to_prim(U_ader_fp_x,W_ader_fp_x,(n_fp*cells_x)*cv_y*cv_z*n_cv);
            fluxes(U_ader_fp_x,W_ader_fp_x,F_ader_fp_x,(n_fp*cells_x)*cv_y*cv_z*n_cv,_vx_,_vy_,_vz_);
            riemann_solver_x();
            derive_fp_x_to_sp(F_ader_fp_x, dU_ader_sp);
            #endif
            #ifdef Y
            transform_sp_to_fp_y(U_ader_sp,U_ader_fp_y,n_cv);
            cons_to_prim(U_ader_fp_y,W_ader_fp_y,cv_x*(n_fp*cells_y)*cv_z*n_cv);
            fluxes(U_ader_fp_y,W_ader_fp_y,F_ader_fp_y,cv_x*(n_fp*cells_y)*cv_z*n_cv,_vy_,_vx_,_vz_);
            riemann_solver_y();
            derive_fp_y_to_sp(F_ader_fp_y, dU_ader_sp);
            #endif
            #ifdef Z
            transform_sp_to_fp_z(U_ader_sp,U_ader_fp_z,n_cv);
            cons_to_prim(U_ader_fp_z,W_ader_fp_z,cv_x*cv_y*(n_fp*cells_z)*n_cv);
            fluxes(U_ader_fp_z,W_ader_fp_z,F_ader_fp_z,cv_x*cv_y*(n_fp*cells_z)*n_cv,_vz_,_vx_,_vy_);
            riemann_solver_z();
            derive_fp_z_to_sp(F_ader_fp_z, dU_ader_sp);
            #endif
            if(i<n){
                ader_subupdate(U_ader_sp,U_sp,dU_ader_sp);
                Boundary_Conditions_ader(U_ader_sp);
            }
            
        }
        #ifdef SD
        ader_update(U_sp,dU_ader_sp);
        Boundary_Conditions(U_sp,1);
        #else
        //Change to Finite Volume scheme
        transform_sp_to_cv(U_sp,U_cv);
        for(i=0;i<n_cv;i++){
            #ifdef X
            face_integral_x(F_ader_fp_x,F_fv_x,i);
            #endif
            #ifdef Y
            face_integral_y(F_ader_fp_y,F_fv_y,i);
            #endif
            #ifdef Z
            face_integral_z(F_ader_fp_z,F_fv_z,i);
            #endif
            //Updates both active and ghost cells of U_new
            fv_update(U_new,U_cv,i);
            Boundary_Conditions(U_new,1);
            //FallBack scheme and trouble detection are performed
            //over primitive variables
            cons_to_prim(U_cv ,W_cv ,size_cv);
            cons_to_prim(U_new ,W_new ,size_cv);
            //Trouble Detection
            detect_troubles();
            //Godunov 2nd 
            godunov_2O();
            //Flux correction//
            //fv_update(U_cv,U_cv,i);
            //fv_godunov2O_update(U_cv,U_cv,i);
            //Updates only active cells of U_cv
            fv_corrected_update(U_cv,U_cv,i);
            //Now we update the ghost cells of U_cv
            Boundary_Conditions(U_cv,1);//Write(n_output++);
        }//Write(n_output++);if(n_step>8)finish();
        transform_cv_to_sp(U_cv,U_sp);
        #endif

        compute_dt();
        
        //Outputs
        if(t==t_output){
            t_output=t+dt_output;
            Write(n_output++);
        }
        if(Master){
            cout<<".";
        }
        if(t+dt>t_output){
            dt=t_output-t;
        }
    }
#ifdef MPI
    MPI_Finalize();
#endif
    return 0;
}
#include "sd3d.hpp"

void fill_comms_x(double *U, double *buffer, int n_a, int side){
    int cv, stride=((1-side)*NGH + side*Nx)*n_cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_a; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<NGH*n_cv;i++){
                        cv = i + j*NGH*n_cv + k*NGH*n_cv*cv_y  + i_a*NGH*n_cv*cv_y*cv_z + var*NGH*n_cv*cv_y*cv_z*n_a;
                        buffer[cv] = U[i+stride  + j*cv_x + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                    }
                }
            }
        }
    }
}

void fill_comms_y(double *U, double *buffer, int n_a, int side){
    int cv, stride=((1-side)*NGH + side*Ny)*n_cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_a; i_a++){
            for(int k=0;k<cv_z;k++){
                for(int j=0;j<NGH*n_cv;j++){
                    for(int i=0;i<cv_x;i++){
                        cv = i + j*cv_x + k*cv_x*NGH*n_cv + i_a*cv_x*NGH*n_cv*cv_z + var*cv_x*NGH*n_cv*cv_z*n_a;
                        buffer[cv] = U[i + (j+stride)*cv_x  + k*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                    }
                }
            }
        }
    }
}

void fill_comms_z(double *U, double *buffer, int n_a, int side){
    int cv, stride=((1-side)*NGH + side*Nz)*n_cv;
    for(int var=0; var<nvar; var++){
        for(int i_a=0; i_a<n_a; i_a++){
            for(int k=0;k<NGH*n_cv;k++){
                for(int j=0;j<cv_y;j++){
                    for(int i=0;i<cv_x;i++){
                        cv = i + j*cv_x + k*cv_x*cv_y + i_a*cv_x*cv_y*NGH*n_cv + var*cv_x*cv_y*NGH*n_cv*n_a;
                        buffer[cv] = U[i + j*cv_x + (k+stride)*cv_x*cv_y + i_a*size_cv + var*size_cv*n_a]; 
                    }
                }
            }
        }
    }
}

void Exec_comms(double *U, int n_a){
    MPI_Request request_send[6], request_recv[6];
    int n_req_s=0, n_req_r=0;
    #ifdef X
    if(cpu_boundary[_x_][0]){
        fill_comms_x(U, comms[_x_][0].buffer_send, n_a, 0);
        MPI_Isend (comms[_x_][0].buffer_send, BC_x_size*n_a, MPI_DOUBLE, comms[_x_][0].rank_dst, 1, Comm, request_send + n_req_s++);
        MPI_Irecv (BC_x[0], BC_x_size*n_a, MPI_DOUBLE, comms[_x_][0].rank_dst, 0, Comm, request_recv + n_req_r++);
    }
    if(cpu_boundary[_x_][1]){
        fill_comms_x(U, comms[_x_][1].buffer_send, n_a, 1);
        MPI_Isend (comms[_x_][1].buffer_send, BC_x_size*n_a, MPI_DOUBLE, comms[_x_][1].rank_dst, 0, Comm, request_send + n_req_s++);
        MPI_Irecv (BC_x[1], BC_x_size*n_a, MPI_DOUBLE, comms[_x_][1].rank_dst, 1, Comm, request_recv + n_req_r++);
    }
    #endif
    #ifdef Y
    if(cpu_boundary[_y_][0]){
        fill_comms_y(U, comms[_y_][0].buffer_send, n_a, 0);
        MPI_Isend (comms[_y_][0].buffer_send, BC_y_size*n_a, MPI_DOUBLE, comms[_y_][0].rank_dst, 1, Comm, request_send + n_req_s++);
        MPI_Irecv (BC_y[0], BC_y_size*n_a, MPI_DOUBLE, comms[_y_][0].rank_dst, 0, Comm, request_recv + n_req_r++);
    }
    if(cpu_boundary[_y_][1]){
        fill_comms_y(U, comms[_y_][1].buffer_send, n_a, 1);
        MPI_Isend (comms[_y_][1].buffer_send, BC_y_size*n_a, MPI_DOUBLE, comms[_y_][1].rank_dst, 0, Comm, request_send + n_req_s++);
        MPI_Irecv (BC_y[1], BC_y_size*n_a, MPI_DOUBLE, comms[_y_][1].rank_dst, 1, Comm, request_recv + n_req_r++);
    }
    #endif
    #ifdef Z
    if(cpu_boundary[_z_][0]){
        fill_comms_z(U, comms[_z_][0].buffer_send, n_a, 0);
        MPI_Isend (comms[_z_][0].buffer_send, BC_z_size*n_a, MPI_DOUBLE, comms[_z_][0].rank_dst, 1, Comm, request_send + n_req_s++);
        MPI_Irecv (BC_z[0], BC_z_size*n_a, MPI_DOUBLE, comms[_z_][0].rank_dst, 0, Comm, request_recv + n_req_r++);
    }
    if(cpu_boundary[_z_][1]){
        fill_comms_z(U, comms[_z_][1].buffer_send, n_a, 1);
        MPI_Isend (comms[_z_][1].buffer_send, BC_z_size*n_a, MPI_DOUBLE, comms[_z_][1].rank_dst, 0, Comm, request_send + n_req_s++);
        MPI_Irecv (BC_z[1], BC_z_size*n_a, MPI_DOUBLE, comms[_z_][1].rank_dst, 1, Comm, request_recv + n_req_r++);
    }
    #endif
    for(int i=0; i<n_req_r; i++)
        MPI_Wait(request_recv+i, MPI_STATUS_IGNORE);
    
    for(int i=0; i<n_req_s; i++)
        MPI_Wait (request_send+i, MPI_STATUS_IGNORE);
        
    MPI_Barrier(Comm);
}
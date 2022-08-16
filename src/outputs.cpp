#include "sd3d.hpp"

void Write_fields(int step){
    transform_sp_to_cv(U_sp,U_cv);
    cons_to_prim(U_cv,W_cv,size_cv);
    ofstream output;
    output.open ("W_cv_"+to_string(cpu_rank)+"_"+to_string(step)+".dat");
    output.write((char*)U_sp, total_size_cv*sizeof(double));
    output.close();
    std::cout<<endl<<"SNAPSHOT "<<step<<" at t="<<t<<std::endl;
}

void Write_array(double *array, int size, char *name){
    ofstream output;
    char outname[100];
    sprintf(outname,"/scratch/snx3000/dvelasco/sd3d/%s_%d.dat", name, cpu_rank);
    output.open(outname);
    output.write((char*)array, size*sizeof(double));
    output.close();
}

void Write_field(double *Field, char const *name, int n_output) {
    FILE *fo;
    char outname[100];
    int relay;
    int i, j, k;
    long offset;
    sprintf(outname,"/scratch/snx3000/dvelasco/sd3d/%s_%d.dat", name, n_output);
    #ifdef MPI
    if (cpu_rank > 0) // Force sequential write
        MPI_Recv (&relay, 1, MPI_INT, cpu_rank-1, 42, Comm, MPI_STATUS_IGNORE);
    #endif

    if(Master) fo = fopen(outname, "w");
    else       fo = fopen(outname, "r+");
    offset = x_i*nx_cv + y_i*ny_cv*(NX*nx_cv) + z_i*nz_cv*(NX*nx_cv)*(NY*ny_cv);
    for(k=0; k<Nz*nz_cv; k++){
        for(j=0; j<Ny*ny_cv; j++){
            offset = (x_i*nx_cv) + (j+y_i*ny_cv)*(NX*nx_cv) + (k+z_i*nz_cv)*(NX*nx_cv)*(NY*ny_cv);
            fseek(fo, offset*sizeof(double), SEEK_SET);
            fwrite(Field + NGHx*n_cv + (j+NGHy*ny_cv)*cv_x + (k+NGHz*nz_cv)*cv_x*cv_y, sizeof(double)*(Nx*nx_cv), 1, fo);
        }
    }
    fclose(fo); 
    #ifdef MPI
    if (cpu_rank < cpu_size-1)  // Force sequential write
    MPI_Send (&relay, 1, MPI_INT, cpu_rank+1, 42, Comm);
    #endif
}

void Write(int output){
    #ifdef SD
    transform_sp_to_cv(U_sp,U_cv);
    #endif
    cons_to_prim(U_cv,W_cv,size_cv);
    Write_field(W_cv, "density", output);
    #ifdef X
    Write_field(W_cv + _vx_*size_cv, "vx", output);
    #endif
    #ifdef Y
    Write_field(W_cv + _vy_*size_cv, "vy", output);
    #endif
    #ifdef Z
    Write_field(W_cv + _vz_*size_cv, "vz", output);
    #endif
    Write_field(W_cv +  _p_*size_cv, "pressure", output);
    #ifndef SD
    Write_field(troubles, "troubles", output);
    #endif
    if(Master)
        cout<<endl<<"SNAPSHOT "<<output<<" at t="<<t<<" dt="<<dt<<endl;
}

void Write_edges(){
    if(!Master)
        return;
    double edge;
    ofstream output;
    #ifdef X
    output.open("/scratch/snx3000/dvelasco/sd3d/X.dat");
    for(int j=0;j<NX;j++){
        for(int i=0;i<n_cv;i++){
            edge = ((j + x_i) + x_fp[i])*dx;
            output.write((char*)&edge, sizeof(double));
        }
    }
    output.write((char*)&boxlen_x, sizeof(double));
    //output.write((char*)&X_faces[NGH*n_cv], (Nx*n_cv+1)*sizeof(double));
    output.close();
    #endif
    #ifdef Y
    output.open("/scratch/snx3000/dvelasco/sd3d/Y.dat");
    for(int j=0;j<NY;j++){
        for(int i=0;i<n_cv;i++){
            edge = ((j + y_i) + x_fp[i])*dy;
            output.write((char*)&edge, sizeof(double));
        }
    }
    output.write((char*)&boxlen_y, sizeof(double));
    //output.write((char*)&Y_faces[NGH*n_cv], (Ny*n_cv+1)*sizeof(double));
    output.close();
    #endif
    #ifdef Z
    output.open("/scratch/snx3000/dvelasco/sd3d/Z.dat");
    for(int j=0;j<NZ;j++){
        for(int i=0;i<n_cv;i++){
            edge = ((j + z_i) + x_fp[i])*dz;
            output.write((char*)&edge, sizeof(double));
        }
    }
    output.write((char*)&boxlen_z, sizeof(double));
    //output.write((char*)Z_faces+NGH*n_cv, (Nz*n_cv+1)*sizeof(double));
    output.close();
    #endif
}
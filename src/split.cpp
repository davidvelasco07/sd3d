#include "sd3d.hpp"

void Split_Domain(){
    #ifdef MPI
    MPI_Barrier(Comm);
    //MPI Processes grid
    #ifndef _3D_
    cpu_x = 1;//int(sqrt(double(cpu_size)));
    cpu_y = cpu_size;///cpu_x;
    int res = cpu_size - cpu_x*cpu_y;
    if(Master)
        cout<<cpu_size<<" -> "<<cpu_x<<","<<cpu_y<<","<<cpu_z<<" -> "<<res<<endl;
    while(res!=0 && cpu_x>1){
        cpu_x--;
        cpu_y = cpu_size/cpu_x;
        res = cpu_size - cpu_x*cpu_y;
        if(Master)
            cout<<cpu_size<<" -> "<<cpu_x<<","<<cpu_y<<","<<cpu_z<<" -> "<<res<<endl;
    }
    #endif
    #ifdef _3D_
    cpu_z = cpu_size;
    cpu_y = 1;
    cpu_x = 1;
    int res=0;
    if(Master)
        cout<<cpu_size<<" -> "<<cpu_x<<","<<cpu_y<<","<<cpu_z<<endl;
    #endif
    //MPI Split domain 
    if (cpu_x > NX || cpu_y > NY || cpu_z > NZ){
        if(Master)
            cout<<endl<<"Error: The number of processes exceeds the number of cells"<<endl;
        finish();
    }
    #ifdef Z
    Nz = NZ/cpu_z;
    res = NZ-Nz*cpu_z;
    rank_z = cpu_rank/(cpu_x*cpu_y);
    if( rank_z < res )
       Nz++;
    z_i = Nz*rank_z;
    if(rank_z >= res)
        z_i+=res;
    #endif
    #ifdef Y
    Ny = NY/cpu_y;  
    res = NY-Ny*cpu_y;
    rank_y = (cpu_rank - rank_z*cpu_x*cpu_y)/cpu_x; 
    if( rank_y < res )
       Ny++;
    y_i = Ny*rank_y;
    if(rank_y >= res)
        y_i+=res;
    #endif
    #ifdef X
    Nx = NX/cpu_x;
    res = NX-Nx*cpu_x;
    rank_x = cpu_rank - rank_y*cpu_x - rank_z*cpu_x*cpu_y;
    if( rank_x < res )
       Nx++;
    x_i = Nx*rank_x;
    if(rank_x >= res)
        x_i+=res;
    #endif
    cout<<cpu_rank<<" -> "<<rank_x<<","<<rank_y<<","<<rank_z<<"  "<<Nx<<","<<Ny<<","<<Nz<<"  "<<x_i<<","<<y_i<<","<<z_i<<endl;
    MPI_Barrier(Comm);
#endif
}
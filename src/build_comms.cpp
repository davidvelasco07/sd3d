#include "sd3d.hpp"

comm_cpu comms[3][2];

void Build_comms(){
#ifdef MPI
    int periodic[3];
    #ifdef Periodic_X
    periodic[_x_]=1;
    #else
    periodic[_x_]=0;
    #endif
    #ifdef Periodic_Y
    periodic[_y_]=1;
    #else
    periodic[_y_]=0;
    #endif
    #ifdef Periodic_Z
    periodic[_z_]=1;
    #else
    periodic[_z_]=0;
    #endif

    for(int j=0; j<3; j++)
        for(int i=0; i<2; i++)
            cpu_boundary[j][i]=0;
   
    if(cpu_size > 1){
#ifdef X
        if(rank_x > 0)
            cpu_neighbour[_x_][0] = cpu_rank - 1;
        else
            cpu_neighbour[_x_][0] = cpu_rank + (cpu_x-1)*periodic[_x_];
        if(rank_x < cpu_x-1)
            cpu_neighbour[_x_][1] = cpu_rank + 1;
        else
            cpu_neighbour[_x_][1] = cpu_rank - (cpu_x-1)*periodic[_x_];
        cpu_boundary[_x_][0] = (cpu_neighbour[_x_][0] == cpu_rank ? 0: 1);
        cpu_boundary[_x_][1] = (cpu_neighbour[_x_][1] == cpu_rank ? 0: 1);
        if(cpu_boundary[_x_][0])
            N_comms+=comms[_x_][0].build(cpu_rank,cpu_neighbour[_x_][0],BC_x_size*n_cv);
        if(cpu_boundary[_x_][1])
            N_comms+=comms[_x_][1].build(cpu_rank,cpu_neighbour[_x_][1],BC_x_size*n_cv);
#endif
#ifdef Y
        if(rank_y > 0)
            cpu_neighbour[_y_][0] = cpu_rank - cpu_x;
        else
            cpu_neighbour[_y_][0] = cpu_rank + cpu_x*(cpu_y-1)*periodic[_y_];
        if(rank_y < cpu_y-1)
            cpu_neighbour[_y_][1] = cpu_rank + cpu_x;
        else
            cpu_neighbour[_y_][1] = cpu_rank - cpu_x*(cpu_y-1)*periodic[_y_];
        cpu_boundary[_y_][0] = (cpu_neighbour[_y_][0] == cpu_rank ? 0: 1);
        cpu_boundary[_y_][1] = (cpu_neighbour[_y_][1] == cpu_rank ? 0: 1);
        if(cpu_boundary[_y_][0])
            N_comms+=comms[_y_][0].build(cpu_rank,cpu_neighbour[_y_][0],BC_y_size*n_cv);
        if(cpu_boundary[_y_][1])
            N_comms+=comms[_y_][1].build(cpu_rank,cpu_neighbour[_y_][1],BC_y_size*n_cv);
#endif
#ifdef Z
        if(rank_z > 0)
            cpu_neighbour[_z_][0] = cpu_rank - cpu_x*cpu_y;
        else
            cpu_neighbour[_z_][0] = cpu_rank + cpu_x*cpu_y*(cpu_z-1)*periodic[_z_];
        if(rank_z < cpu_z-1)
            cpu_neighbour[_z_][1] = cpu_rank + cpu_x*cpu_y;
        else
            cpu_neighbour[_z_][1] = cpu_rank - cpu_x*cpu_y*(cpu_z-1)*periodic[_z_];
        cpu_boundary[_z_][0] = (cpu_neighbour[_z_][0] == cpu_rank ? 0: 1);
        cpu_boundary[_z_][1] = (cpu_neighbour[_z_][1] == cpu_rank ? 0: 1);
        if(cpu_boundary[_z_][0])
            N_comms+=comms[_z_][0].build(cpu_rank,cpu_neighbour[_z_][0],BC_z_size*n_cv);
        if(cpu_boundary[_z_][1])
            N_comms+=comms[_z_][1].build(cpu_rank,cpu_neighbour[_z_][1],BC_z_size*n_cv);
#endif
        //cout<<cpu_rank<<" has "<<N_comms<<" comms"<<endl;
    }
#endif
}

#include <iostream>
#include <mpi.h>
#include "common.h"
#include <math.h>
#include <vector>
#include<assert.h>


#include<stdio.h>

int GetMPIrank(){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    return rank;
}
unsigned int getGridIndex(const unsigned int i, const unsigned int j, const unsigned int NumGrid){
    // assert((i < NumGrid) and (i >= 0));
    // assert((j < NumGrid) and (j >= 0));
    return(i*NumGrid + j);
}
//

void bin_particles(std::vector<particle_t> & local_particles, particle_t * particles, grid_t * grid, double  start_proc_boundary,double end_proc_boundary, int n, int NumGrid, double dx){
    local_particles.clear();
    unsigned int idx,idy,id;
    int counter(0);

    for(int i = 0; i < n;i++){
        if((particles[i].x >= start_proc_boundary) and (particles[i].x < end_proc_boundary)){
            double x = particles[i].x - start_proc_boundary;
            idx = static_cast<unsigned  int>(x * dx);
            idy = static_cast<unsigned int>(particles[i].y * dx);
            id  = (idx*NumGrid + idy);
            grid[id].id_[grid[id].count_] = counter;
            particles[i].gid =id;

            counter++;
            grid[id].count_++;
            local_particles.push_back(particles[i]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

}
void clear_bins(grid_t* grid,int size){
    for(int i = 0; i < size; i++){
        grid[i].count_ = 0;
    }
}

void erase_buffer_particles(grid_t * grid,std::vector<particle_t> & localParticles, int rank, int nrows, int size, int NumGrid){
    int removeBuffer[2];
    removeBuffer[0] = NumGrid;
    removeBuffer[1] = (nrows - 1)*NumGrid;
    if(rank == 0){
        removeBuffer[0] = -1;
    }
    if(rank == size - 1){
        removeBuffer[1] = -1;
    }
//    int locSize = localParticles.size();

    for(int i = 0; i < localParticles.size();i++){
        if(((int)localParticles[i].gid < (removeBuffer[0]) or (localParticles[i].gid >= removeBuffer[1]))){
            std::iter_swap(&localParticles[i],&localParticles[localParticles.size() - 1]);
            localParticles.pop_back();
            i = i - 1;
        }
    }
//    if(removeBuffer[0] != -1){
//        for(int i = 0; i < NumGrid; i++){
//            for(int j = 0; j < grid[i].count_; j++) {
//                std::iter_swap(&localParticles[grid[i].id_[j]],&localParticles[localParticles.size() - 1]);
//                unsigned int partGid = localParticles[localParticles.size() - 1].gid;
//                for(int gid = 0; gid < grid[partGid].)
//                localParticles.pop_back();
//            }
//        }
//    }
//    if(removeBuffer[1] != -1){
//        for(int i = removeBuffer[1]*NumGrid; i < (removeBuffer[1] + 1)*NumGrid; i++){
//            for(int j = 0; j < grid[i].count_; j++) {
//                std::iter_swap(&localParticles[grid[i].id_[j]],&localParticles[localParticles.size() - 1]);
//                localParticles.pop_back();
//            }
//        }
//    }

}

void checkWithinBB(std::vector<particle_t>& particles, double startBB, double endBB){
    for(int i = 0; i < particles.size(); i++){
        assert((particles[i].x >= startBB) and (particles[i].x <= endBB));
    }
}
int main(int argc, char *argv[])
{

    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
 
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    int rank, n_proc;
    MPI_Comm comm;
    int dim[2];
    int period[2], reorder;
    int coord[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);


    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;


    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );

    if( rank == 0 )
    init_particles( n, particles );
    // MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(particles,n,PARTICLE,0,MPI_COMM_WORLD);


    MPI_Barrier(MPI_COMM_WORLD);
    unsigned int NumGrid = 16;//ceil(get_size()/0.05);



    int num_row_per_proc = ceil(NumGrid/(n_proc*1.0));
    NumGrid = n_proc*num_row_per_proc;
    const double dx = NumGrid/get_size();
    if(rank == 0) {
        std::cout << "Total Number of Grid = " << NumGrid << "\n";
        std::cout << "num particle = " << n << "\n";
        std::cout << "Dx = " << get_size()/NumGrid << "\n";
    }
    assert(dx > cutoff);

    int start_row = max((rank)*num_row_per_proc - 1,0);
    int end_row = min((rank+1)*num_row_per_proc+1,NumGrid);
  
    int num_row = num_row_per_proc+2;

    if((rank == 0) or (rank == n_proc - 1)){
        num_row = num_row_per_proc + 1;
    }
    double start_proc_boundary = start_row / dx;
    double end_proc_boundary = end_row / dx ;

    int startGrid = 1*NumGrid;
    double realProcBoundaryStart = (start_row+1)/dx;
    int endGrid = (num_row_per_proc + 1)*NumGrid;
    double realProcBoundaryEnd   = (end_row - 1)/dx;
    if(rank == 0){
        startGrid = 0;
        endGrid = (num_row_per_proc)*NumGrid;
        realProcBoundaryStart = 0;
    }
    if(rank == n_proc - 1){
        realProcBoundaryEnd = end_row/dx;
    }
    
     

    int exchange_up = ((rank+1) < n_proc )?(rank + 1) : MPI_PROC_NULL;
    int exchange_down = ((rank-1) >= 0 )?(rank - 1) : MPI_PROC_NULL;
    const int ParticlePerGrid = ceil(n/(NumGrid));
    MPI_Barrier(MPI_COMM_WORLD);

    for(int i = 0; i < n_proc; i++){
        if(rank == i){
            std::cout << "Rank = " << rank << " Start Row = " << start_row <<  " End row = " << end_row -1<< 
            " " << start_proc_boundary << " " << end_proc_boundary << " " << " " << realProcBoundaryStart <<
            " " << realProcBoundaryEnd << " " <<exchange_up << " " << exchange_down << "\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    grid_t * grid = (grid_t *) malloc(NumGrid* num_row*sizeof(grid_t));


    // particle_t *local_buffer_up =  (particle_t*) malloc( ParticlePerGrid*NumGrid * sizeof(particle_t) );
    // particle_t *local_buffer_down =  (particle_t*) malloc( ParticlePerGrid*NumGrid * sizeof(particle_t) );
    for(int i = 0;i < NumGrid*num_row; i++){
        grid[i].id_ = (int *) malloc(ParticlePerGrid*2* sizeof(int));
        /**Computing neigbourhood index of grid **/
        unsigned int idx = i/NumGrid;
        unsigned int idy = i%NumGrid;
        // if(rank == 0 and i == 7){
        //     std::cout << idx << " " << idy << " "<< NumGrid << " " << num_row << "\n";
        //     assert(false);
        // }
        for(int j = 0; j < 8; j++){
            grid[i].neighbours_[j] = -1;
        }
        if((idx == 0) and (idy == 0)){
            grid[i].neighbours_[N]  = getGridIndex(idx  ,idy+1,NumGrid) ;
            grid[i].neighbours_[NE] = getGridIndex(idx+1,idy+1,NumGrid) ;
            grid[i].neighbours_[E]  = getGridIndex(idx+1,idy  ,NumGrid) ;
        }
        else if((idx == num_row - 1) and (idy == 0)){
            grid[i].neighbours_[W]  = getGridIndex(idx-1,idy  ,NumGrid) ;
            grid[i].neighbours_[NW] = getGridIndex(idx-1,idy+1,NumGrid) ;
            grid[i].neighbours_[N]  = getGridIndex(idx  ,idy+1,NumGrid) ;
        }
        else if((idx == 0) and (idy == NumGrid - 1)){
            grid[i].neighbours_[S ] = getGridIndex(idx  ,idy-1,NumGrid);
            grid[i].neighbours_[SE] = getGridIndex(idx+1,idy-1,NumGrid);
            grid[i].neighbours_[E]  = getGridIndex(idx+1,idy  ,NumGrid);
        }
        else if((idx == num_row - 1) and (idy == NumGrid - 1)){
            grid[i].neighbours_[S ] = getGridIndex(idx  ,idy-1,NumGrid);
            grid[i].neighbours_[SW] = getGridIndex(idx-1,idy-1,NumGrid);
            grid[i].neighbours_[W]  = getGridIndex(idx-1,idy  ,NumGrid);
        }
        else if(idx == 0){
            grid[i].neighbours_[N ] = getGridIndex(idx  ,idy+1,NumGrid);
            grid[i].neighbours_[NE] = getGridIndex(idx+1,idy+1,NumGrid);
            grid[i].neighbours_[E]  = getGridIndex(idx+1,idy  ,NumGrid);
            grid[i].neighbours_[SE] = getGridIndex(idx+1,idy-1,NumGrid);
            grid[i].neighbours_[S]  = getGridIndex(idx  ,idy-1,NumGrid);
        }
        else if(idx == (num_row - 1)){
            grid[i].neighbours_[N ] = getGridIndex(idx  ,idy+1,NumGrid);
            grid[i].neighbours_[NW] = getGridIndex(idx-1,idy+1,NumGrid);
            grid[i].neighbours_[W]  = getGridIndex(idx-1,idy  ,NumGrid);
            grid[i].neighbours_[SW] = getGridIndex(idx-1,idy-1,NumGrid);
            grid[i].neighbours_[S]  = getGridIndex(idx  ,idy-1,NumGrid);
        }
        else if(idy == 0){
            grid[i].neighbours_[N ] = getGridIndex(idx  ,idy+1,NumGrid);
            grid[i].neighbours_[NW] = getGridIndex(idx-1,idy+1,NumGrid);
            grid[i].neighbours_[W]  = getGridIndex(idx-1,idy  ,NumGrid);
            grid[i].neighbours_[NE] = getGridIndex(idx+1,idy+1,NumGrid);
            grid[i].neighbours_[E]  = getGridIndex(idx+1,idy  ,NumGrid);
        }
        else if(idy == (NumGrid - 1)){
            grid[i].neighbours_[E]  = getGridIndex(idx+1,idy  ,NumGrid);
            grid[i].neighbours_[SE] = getGridIndex(idx+1,idy-1,NumGrid);
            grid[i].neighbours_[W]  = getGridIndex(idx-1,idy  ,NumGrid);
            grid[i].neighbours_[SW] = getGridIndex(idx-1,idy-1,NumGrid);
            grid[i].neighbours_[S]  = getGridIndex(idx  ,idy-1,NumGrid);
        }
        else{
            grid[i].neighbours_[E]  = getGridIndex(idx+1,idy  ,NumGrid);
            grid[i].neighbours_[SE] = getGridIndex(idx+1,idy-1,NumGrid);
            grid[i].neighbours_[W]  = getGridIndex(idx-1,idy  ,NumGrid);
            grid[i].neighbours_[SW] = getGridIndex(idx-1,idy-1,NumGrid);
            grid[i].neighbours_[S]  = getGridIndex(idx  ,idy-1,NumGrid);
            grid[i].neighbours_[N ] = getGridIndex(idx  ,idy+1,NumGrid);
            grid[i].neighbours_[NW] = getGridIndex(idx-1,idy+1,NumGrid);
            grid[i].neighbours_[NE] = getGridIndex(idx+1,idy+1,NumGrid);
        }

    }
    if (rank == n_proc - 1){
        // std::cout << num_row << "\n";
        
        for(int i = 0;i < NumGrid; i++){
            grid[(num_row - 2)*NumGrid + i].neighbours_[N ] = -1;
            grid[(num_row - 2)*NumGrid + i].neighbours_[NW] = -1;
            grid[(num_row - 2)*NumGrid + i].neighbours_[NE] = -1;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
   
    // std::cout << "In n = " << n << "\n";
   
    clear_bins(grid,num_row*NumGrid);
    std::vector<particle_t> local_particles;
    bin_particles(local_particles,particles,grid,start_proc_boundary,end_proc_boundary,n,NumGrid, dx);
    MPI_Barrier(MPI_COMM_WORLD);
    free(particles);
    int sum(0);
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0; i < n_proc; i++){
        if(rank == i){
            std::cout << "Rank = " << rank << " Start Grid = " << startGrid <<  " End Grid = " << endGrid<<"\n";

        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i = startGrid; i < endGrid ; i++){
            sum += grid[i].count_;
            for(int j = 0; j < grid[i].count_ ; j++){
                double x = local_particles[grid[i].id_[j]].x;
                if(not((x >= realProcBoundaryStart ) and (x <= realProcBoundaryEnd))){
                    std::cout << "X = " << x << " " << rank << " " << realProcBoundaryStart << " " << realProcBoundaryEnd << "\n";
                }
                assert((x >= realProcBoundaryStart ) and (x <= realProcBoundaryEnd));
            }
            
            
    }

    int sum_all;
    MPI_Reduce(&sum,&sum_all,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

   
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0){
        assert(sum_all == n);
    }



    /**Checking particle exchange **/
/*    std::vector<particle_t> sendParticles;
    particle_t *receive_particle = (particle_t*) malloc( 100 * sizeof(particle_t) );
    particle_t sp1,sp2;
    sp1.x = rank; sp1.y = rank; sp1.vx = -rank; sp1.vy = -rank;  
    sp2.x = 100*rank; sp2.y = 100*rank; sp2.vx = -100*rank; sp2.vy = -100*rank;  
    sendParticles.push_back(sp1);
    sendParticles.push_back(sp2);
    MPI_Request request[2];
    MPI_Status status[2];
    MPI_Isend(sendParticles.data(),sendParticles.size(),PARTICLE,exchange_up,sendParticles.size(),MPI_COMM_WORLD,&request[0]);
    std::cout << rank << " " << "Sending to  " << exchange_up << "\n";
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << rank << " " << "Receiving from  " << exchange_down << "\n";
    MPI_Recv(&receive_particle[0],100,PARTICLE,exchange_down,MPI_ANY_TAG,MPI_COMM_WORLD,&status[0]);
    if(rank == 2)
        std::cout << "Received particle = " << rank << " " << receive_particle[0].x << " " << receive_particle[0].y << " " 
        << receive_particle[1].x << " " << receive_particle[1].y << " "  << status[0].MPI_TAG << " " << status[0].MPI_SOURCE << "\n";*/
    

    for (int step = 0; step < 1; step++) {
        for(int i = startGrid; i < endGrid ; i++){
        for(int part = 0; part < grid[i].count_; part++){
                local_particles[grid[i].id_[part]].ax = 0;
                local_particles[grid[i].id_[part]].ay = 0;
                for(int partGrid = 0; partGrid < grid[i].count_; partGrid++){
                    apply_force(local_particles[grid[i].id_[part]],local_particles[grid[i].id_[partGrid]], &dmin, &davg, &navg);
                }
                for(int neigh = 0; neigh < 8; neigh++){
                    unsigned int neighbour_id =grid[i].neighbours_[neigh];
                    if( neighbour_id != -1){
                        for(int part_neigh = 0; part_neigh < grid[neighbour_id].count_; part_neigh++){
                            apply_force(local_particles[grid[i].id_[part]], local_particles[grid[neighbour_id].id_[part_neigh]], &dmin, &davg, &navg);
                        }
                    }
                }
            }
        }
        erase_buffer_particles(grid,local_particles,rank,num_row,n_proc,NumGrid);
        int num_part;
        int local_part = static_cast<int>(local_particles.size());
        MPI_Reduce(&local_part,&num_part,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
        if(rank == 0)
        assert(n == num_part);
        checkWithinBB(local_particles,realProcBoundaryStart,realProcBoundaryEnd);

        for(int part = 0; part < local_particles.size(); part++){
            move(local_particles[part]);
        }



    }
   
    //     std::cout << counter << "\n";
    //     counter = 0;
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     int count_down = 0;
    //     for(int i = startGrid; i < (num_row -1)*NumGrid ; i++){
    //         for(int part = 0; part < grid[i].count_; part++){
    //             counter++;
    //             move(particles[grid[i].id_[part]]);
    //             if(particles[grid[i].id_[part]].x < start_proc_boundary){
    //                 local_buffer_send_down.push_back(&particles[grid[i].id_[part]]);
    //             }
    //             else if(particles[grid[i].id_[part]].x > end_proc_boundary){
    //                 local_buffer_send_up.push_back(&particles[grid[i].id_[part]]);
    //             }


    //         }
    //     }
    //     std::cout << counter << "\n";

    // }

    
    MPI_Finalize();
    return 0;
}
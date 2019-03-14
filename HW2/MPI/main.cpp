#include <iostream>
#include <mpi.h>
#include "common.h"
#include <math.h>
#include <vector>
#include<assert.h>
#include <fstream>

#include<stdio.h>

int GetMPIrank() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

unsigned int getGridIndex(const unsigned int i, const unsigned int j, const unsigned int NumGrid) {
    // assert((i < NumGrid) and (i >= 0));
    // assert((j < NumGrid) and (j >= 0));
    return (i * NumGrid + j);
}
//

void bin_particles_start(std::vector<particle_t> &local_particles, particle_t *particles, grid_t *grid,
                         double start_proc_boundary, double end_proc_boundary, int n, int NumGrid, double dx) {
    local_particles.clear();
    unsigned int idx, idy, id;
    int counter(0);

    for (int i = 0; i < n; i++) {
        if ((particles[i].x >= start_proc_boundary) and (particles[i].x < end_proc_boundary)) {
            double x = particles[i].x - start_proc_boundary;
            idx = static_cast<unsigned int>(x * dx);
            idy = static_cast<unsigned int>(particles[i].y * dx);
            id = (idx * NumGrid + idy);
            grid[id].id_[grid[id].count_] = counter;
            particles[i].gid = id;

            counter++;
            grid[id].count_++;
            local_particles.push_back(particles[i]);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

}

void
bin_particles(std::vector<particle_t> &local_particles, grid_t *grid, double start_proc_boundary,
              double end_proc_boundary, int NumGrid, double OneByDx) {

    unsigned int idx, idy, id;
    int counter(0);

    for (int i = 0; i < local_particles.size(); i++) {
        assert((local_particles[i].x >= start_proc_boundary) and (local_particles[i].x <= end_proc_boundary));
        double x = local_particles[i].x - start_proc_boundary;
        idx = static_cast<unsigned int>(x * OneByDx);
        idy = static_cast<unsigned int>(local_particles[i].y * OneByDx);
        id = (idx * NumGrid + idy);
        grid[id].id_[grid[id].count_] = counter;
        local_particles[i].gid = id;

        counter++;
        grid[id].count_++;
    }
    MPI_Barrier(MPI_COMM_WORLD);

}

void
binBufferparticles(std::vector<particle_t> &local_particles, grid_t *grid, double start_proc_boundary,
                   double end_proc_boundary, int NumGrid, double OneByDx, int startN) {

    unsigned int idx, idy, id;


    for (int i = startN; i < local_particles.size(); i++) {
        assert((local_particles[i].x >= start_proc_boundary) and (local_particles[i].x <= end_proc_boundary));
        double x = local_particles[i].x - start_proc_boundary;
        idx = static_cast<unsigned int>(x * OneByDx);
        idy = static_cast<unsigned int>(local_particles[i].y * OneByDx);
        id = (idx * NumGrid + idy);
        grid[id].id_[grid[id].count_] = i;
        local_particles[i].gid = id;


        grid[id].count_++;
    }
    MPI_Barrier(MPI_COMM_WORLD);

}

void clear_bins(grid_t *grid, int size) {
    for (int i = 0; i < size; i++) {
        grid[i].count_ = 0;
    }
}

void erase_buffer_particles(grid_t *grid, std::vector<particle_t> &localParticles, int rank, int nrows, int size,
                            int NumGrid, double RealBoundaryStart, double RealBoundaryEnd, int startN) {

    for (int i = 0; i < localParticles.size(); i++) {
        if ((localParticles[i].x < RealBoundaryStart) or (localParticles[i].x >= RealBoundaryEnd)) {

            std::iter_swap(&localParticles[i], &localParticles[localParticles.size() - 1]);
            localParticles.pop_back();
            i = i - 1;
        }
    }



}

void packBuffers(std::vector<particle_t> &local_particles, std::vector<particle_t> &top_buffer,
                 std::vector<particle_t> &down_buffer,
                 double bufferStartBoundary, double bufferEndBoundary) {
    top_buffer.clear();
    down_buffer.clear();
    for (int i = 0; i < local_particles.size(); i++) {
        if (local_particles[i].x < bufferStartBoundary) {
            down_buffer.push_back(local_particles[i]);
        } else if (local_particles[i].x > bufferEndBoundary) {
            top_buffer.push_back(local_particles[i]);
        }
        ///TODO: put binning here
    }


}

void packBuffersandBin(grid_t *grid, std::vector<particle_t> &local_particles, std::vector<particle_t> &top_buffer,
                       std::vector<particle_t> &down_buffer,
                       double bufferStartBoundary, double bufferEndBoundary, double start_proc_boundary,
                       double end_proc_boundary,
                       double OneByDx, unsigned int NumGrid) {
    top_buffer.clear();
    down_buffer.clear();
    unsigned int idx, idy, id;
    int counter(0);

    for (int i = 0; i < local_particles.size(); i++) {
        if (local_particles[i].x <= bufferStartBoundary) {
            down_buffer.push_back(local_particles[i]);
        }
        if (local_particles[i].x >= bufferEndBoundary) {
            top_buffer.push_back(local_particles[i]);
        }
        /*** Binning step *******************/
        if ((local_particles[i].x >= start_proc_boundary) and (local_particles[i].x <= end_proc_boundary)) {
            double x = local_particles[i].x - start_proc_boundary;
            idx = static_cast<unsigned int>(x * OneByDx);
            idy = static_cast<unsigned int>(local_particles[i].y * OneByDx);
            id = (idx * NumGrid + idy);
            grid[id].id_[grid[id].count_] = counter;
            local_particles[i].gid = id;
            counter++;
            grid[id].count_++;
        }

    }


}


void checkWithinBB(std::vector<particle_t> &particles, double startBB, double endBB) {
    for (int i = 0; i < particles.size(); i++) {
        if (not(((particles[i].x >= startBB) and (particles[i].x <= endBB)))) {
            std::cout << GetMPIrank() << particles[i].x << " " << particles[i].y << " " << particles[i].vx << " "
                      << particles[i].vy <<
                      " " << startBB << " " << endBB << "\n";
        }
        assert((particles[i].x >= startBB) and (particles[i].x <= endBB));
    }
}

int main(int argc, char *argv[]) {

    int navg, nabsavg = 0;
    double dmin, absmin = 1.0, davg, absavg = 0.0;
    double rdavg, rdmin;
    int rnavg;

    //
    //  process command line parameters
    //
    if (find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        printf("-s <filename> to specify a summary file name\n");
        printf("-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", NULL);
    char *sumname = read_string(argc, argv, "-s", NULL);

    int rank, n_proc;
    MPI_Comm comm;
    int dim[2];
    int period[2], reorder;
    int coord[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);


    FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen(sumname, "a") : NULL;


    particle_t *particles = (particle_t *) malloc(n * sizeof(particle_t));

    MPI_Datatype PARTICLE;
    MPI_Type_contiguous(7, MPI_DOUBLE, &PARTICLE);
    MPI_Type_commit(&PARTICLE);

    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int *) malloc((n_proc + 1) * sizeof(int));
    for (int i = 0; i < n_proc + 1; i++)
        partition_offsets[i] = min(i * particle_per_proc, n);

    int *partition_sizes = (int *) malloc(n_proc * sizeof(int));
    for (int i = 0; i < n_proc; i++)
        partition_sizes[i] = partition_offsets[i + 1] - partition_offsets[i];

    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size(n);

    if (rank == 0)
        init_particles(n, particles);
    // MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);


    MPI_Barrier(MPI_COMM_WORLD);
    unsigned int NumGrid = ceil(get_size() / 0.05);


    int num_row_per_proc = ceil(NumGrid / (n_proc * 1.0));
    NumGrid = n_proc * num_row_per_proc;
    const double OneByDx = NumGrid / get_size();
    const double dx = 1.0 / OneByDx;
//    if (rank == 0) {
//        std::cout << "Total Number of Grid = " << NumGrid << "\n";
//        std::cout << "num particle = " << n << "\n";
//        std::cout << "Dx = " << get_size() / NumGrid << "\n";
//    }
    assert(OneByDx > cutoff);

    int start_row = max((rank) * num_row_per_proc - 1, 0);
    int end_row = min((rank + 1) * num_row_per_proc + 1, NumGrid);

    int num_row = num_row_per_proc + 2;

    if ((rank == 0) or (rank == n_proc - 1)) {
        num_row = num_row_per_proc + 1;
    }
    double start_proc_boundary = start_row / OneByDx;
    double end_proc_boundary = end_row / OneByDx;

    int startGrid = 1 * NumGrid;
    double realProcBoundaryStart = (start_row + 1) / OneByDx;
    int endGrid = (num_row_per_proc + 1) * NumGrid;
    double realProcBoundaryEnd = (end_row - 1) / OneByDx;
    if (rank == 0) {
        startGrid = 0;
        endGrid = (num_row_per_proc) * NumGrid;
        realProcBoundaryStart = 0;
    }
    if (rank == n_proc - 1) {
        realProcBoundaryEnd = end_row / OneByDx;
    }

    double exchangeBoundaryDown = realProcBoundaryStart + dx;
    double exchangeBoundaryUp = realProcBoundaryEnd - dx;
    if (rank == 0) {
        exchangeBoundaryDown = 0;
    } else if (rank == (n_proc - 1)) {
        exchangeBoundaryUp = get_size();
    }


    int exchange_up = ((rank + 1) < n_proc) ? (rank + 1) : MPI_PROC_NULL;
    int exchange_down = ((rank - 1) >= 0) ? (rank - 1) : MPI_PROC_NULL;
    const int ParticlePerGrid = ceil(n / (NumGrid));
    MPI_Barrier(MPI_COMM_WORLD);

//    for (int i = 0; i < n_proc; i++) {
//        if (rank == i) {
//            std::cout << "Rank = " << rank << " Start Row = " << start_row << " End row = " << end_row - 1 <<
//                      " " << start_proc_boundary << " " << end_proc_boundary << " " << " " << realProcBoundaryStart <<
//                      " " << realProcBoundaryEnd << " " << exchange_up << " " << exchange_down << " "
//                      << exchangeBoundaryUp << " "
//                      << exchangeBoundaryDown << "\n";
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);

    grid_t *grid = (grid_t *) malloc(NumGrid * num_row * sizeof(grid_t));


    for (int i = 0; i < NumGrid * num_row; i++) {
        grid[i].id_ = (int *) malloc(ParticlePerGrid * 2 * sizeof(int));
        /**Computing neigbourhood index of grid **/
        unsigned int idx = i / NumGrid;
        unsigned int idy = i % NumGrid;
        // if(rank == 0 and i == 7){
        //     std::cout << idx << " " << idy << " "<< NumGrid << " " << num_row << "\n";
        //     assert(false);
        // }
        for (int j = 0; j < 8; j++) {
            grid[i].neighbours_[j] = -1;
        }
        if ((idx == 0) and (idy == 0)) {
            grid[i].neighbours_[N] = getGridIndex(idx, idy + 1, NumGrid);
            grid[i].neighbours_[NE] = getGridIndex(idx + 1, idy + 1, NumGrid);
            grid[i].neighbours_[E] = getGridIndex(idx + 1, idy, NumGrid);
        } else if ((idx == num_row - 1) and (idy == 0)) {
            grid[i].neighbours_[W] = getGridIndex(idx - 1, idy, NumGrid);
            grid[i].neighbours_[NW] = getGridIndex(idx - 1, idy + 1, NumGrid);
            grid[i].neighbours_[N] = getGridIndex(idx, idy + 1, NumGrid);
        } else if ((idx == 0) and (idy == NumGrid - 1)) {
            grid[i].neighbours_[S] = getGridIndex(idx, idy - 1, NumGrid);
            grid[i].neighbours_[SE] = getGridIndex(idx + 1, idy - 1, NumGrid);
            grid[i].neighbours_[E] = getGridIndex(idx + 1, idy, NumGrid);
        } else if ((idx == num_row - 1) and (idy == NumGrid - 1)) {
            grid[i].neighbours_[S] = getGridIndex(idx, idy - 1, NumGrid);
            grid[i].neighbours_[SW] = getGridIndex(idx - 1, idy - 1, NumGrid);
            grid[i].neighbours_[W] = getGridIndex(idx - 1, idy, NumGrid);
        } else if (idx == 0) {
            grid[i].neighbours_[N] = getGridIndex(idx, idy + 1, NumGrid);
            grid[i].neighbours_[NE] = getGridIndex(idx + 1, idy + 1, NumGrid);
            grid[i].neighbours_[E] = getGridIndex(idx + 1, idy, NumGrid);
            grid[i].neighbours_[SE] = getGridIndex(idx + 1, idy - 1, NumGrid);
            grid[i].neighbours_[S] = getGridIndex(idx, idy - 1, NumGrid);
        } else if (idx == (num_row - 1)) {
            grid[i].neighbours_[N] = getGridIndex(idx, idy + 1, NumGrid);
            grid[i].neighbours_[NW] = getGridIndex(idx - 1, idy + 1, NumGrid);
            grid[i].neighbours_[W] = getGridIndex(idx - 1, idy, NumGrid);
            grid[i].neighbours_[SW] = getGridIndex(idx - 1, idy - 1, NumGrid);
            grid[i].neighbours_[S] = getGridIndex(idx, idy - 1, NumGrid);
        } else if (idy == 0) {
            grid[i].neighbours_[N] = getGridIndex(idx, idy + 1, NumGrid);
            grid[i].neighbours_[NW] = getGridIndex(idx - 1, idy + 1, NumGrid);
            grid[i].neighbours_[W] = getGridIndex(idx - 1, idy, NumGrid);
            grid[i].neighbours_[NE] = getGridIndex(idx + 1, idy + 1, NumGrid);
            grid[i].neighbours_[E] = getGridIndex(idx + 1, idy, NumGrid);
        } else if (idy == (NumGrid - 1)) {
            grid[i].neighbours_[E] = getGridIndex(idx + 1, idy, NumGrid);
            grid[i].neighbours_[SE] = getGridIndex(idx + 1, idy - 1, NumGrid);
            grid[i].neighbours_[W] = getGridIndex(idx - 1, idy, NumGrid);
            grid[i].neighbours_[SW] = getGridIndex(idx - 1, idy - 1, NumGrid);
            grid[i].neighbours_[S] = getGridIndex(idx, idy - 1, NumGrid);
        } else {
            grid[i].neighbours_[E] = getGridIndex(idx + 1, idy, NumGrid);
            grid[i].neighbours_[SE] = getGridIndex(idx + 1, idy - 1, NumGrid);
            grid[i].neighbours_[W] = getGridIndex(idx - 1, idy, NumGrid);
            grid[i].neighbours_[SW] = getGridIndex(idx - 1, idy - 1, NumGrid);
            grid[i].neighbours_[S] = getGridIndex(idx, idy - 1, NumGrid);
            grid[i].neighbours_[N] = getGridIndex(idx, idy + 1, NumGrid);
            grid[i].neighbours_[NW] = getGridIndex(idx - 1, idy + 1, NumGrid);
            grid[i].neighbours_[NE] = getGridIndex(idx + 1, idy + 1, NumGrid);
        }

    }
    if (rank == n_proc - 1) {
        // std::cout << num_row << "\n";

        for (int i = 0; i < NumGrid; i++) {
//            grid[(num_row - 2) * NumGrid + i].neighbours_[N] = -1;
//            grid[(num_row - 2) * NumGrid + i].neighbours_[NW] = -1;
//            grid[(num_row - 2) * NumGrid + i].neighbours_[NE] = -1;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // std::cout << "In n = " << n << "\n";

    clear_bins(grid, num_row * NumGrid);
    std::vector<particle_t> local_particles;
    bin_particles_start(local_particles, particles, grid, start_proc_boundary, end_proc_boundary, n, NumGrid, OneByDx);
    MPI_Barrier(MPI_COMM_WORLD);
    free(particles);
    int sum(0);
    MPI_Barrier(MPI_COMM_WORLD);
 /*   for (int i = 0; i < n_proc; i++) {
        if (rank == i) {
            std::cout << "Rank = " << rank << " Start Grid = " << startGrid << " End Grid = " << endGrid << "\n";

        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = startGrid; i < endGrid; i++) {
        sum += grid[i].count_;
        for (int j = 0; j < grid[i].count_; j++) {
            double x = local_particles[grid[i].id_[j]].x;
            if (not((x >= realProcBoundaryStart) and (x <= realProcBoundaryEnd))) {
                std::cout << "X = " << x << " " << rank << " " << realProcBoundaryStart << " " << realProcBoundaryEnd
                          << "\n";
            }
            assert((x >= realProcBoundaryStart) and (x <= realProcBoundaryEnd));
        }


    }*/


   /* int sum_all;
    MPI_Reduce(&sum, &sum_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        assert(sum_all == n);
    }*/

    std::vector<particle_t> send_up, send_down;
    unsigned int exchangeParticleSize = ParticlePerGrid * 2 * NumGrid;
    particle_t *receive_up = (particle_t *) malloc(exchangeParticleSize * sizeof(particle_t));
    particle_t *receive_down = (particle_t *) malloc(exchangeParticleSize * sizeof(particle_t));
    /**Checking particle exchange **/


    MPI_Request request[2];
    MPI_Status status[2];
    int startN = 0;
    for (int step = 0; step < NSTEPS; step++) {
        MPI_Irecv(&receive_up[0], exchangeParticleSize, PARTICLE, exchange_up, MPI_ANY_TAG, MPI_COMM_WORLD,&request[0]);
        MPI_Irecv(&receive_down[0], exchangeParticleSize, PARTICLE, exchange_down, MPI_ANY_TAG, MPI_COMM_WORLD,&request[1]);
        for (int i = startGrid; i < endGrid; i++) {

            for (int part = 0; part < grid[i].count_; part++) {
                local_particles[grid[i].id_[part]].ax = 0;
                local_particles[grid[i].id_[part]].ay = 0;

                for (int partGrid = 0; partGrid < grid[i].count_; partGrid++) {

                    apply_force(local_particles[grid[i].id_[part]], local_particles[grid[i].id_[partGrid]], &dmin,
                                &davg, &navg);
                }
                for (int neigh = 0; neigh < 8; neigh++) {
                    unsigned int neighbour_id = grid[i].neighbours_[neigh];
                    if (neighbour_id != -1) {
                        for (int part_neigh = 0; part_neigh < grid[neighbour_id].count_; part_neigh++) {
                            apply_force(local_particles[grid[i].id_[part]],
                                        local_particles[grid[neighbour_id].id_[part_neigh]], &dmin, &davg, &navg);

                        }
                    }
                }
            }
        }

        erase_buffer_particles(grid, local_particles, rank, num_row, n_proc, NumGrid, realProcBoundaryStart,
                               realProcBoundaryEnd,startN);


        for (int part = 0; part < local_particles.size(); part++) {
            move(local_particles[part]);
        }
        clear_bins(grid, num_row * NumGrid);
        packBuffersandBin(grid, local_particles, send_up, send_down, exchangeBoundaryDown, exchangeBoundaryUp,
                          start_proc_boundary, end_proc_boundary, OneByDx, NumGrid);


        startN = local_particles.size();

        MPI_Send(send_up.data(), send_up.size(), PARTICLE, exchange_up, send_up.size(), MPI_COMM_WORLD);
        MPI_Send(send_down.data(), send_down.size(), PARTICLE, exchange_down, send_down.size(), MPI_COMM_WORLD);
        MPI_Wait(&request[0],&status[0]);
        MPI_Wait(&request[1],&status[1]);

        for (int i = 0; i < status[0].MPI_TAG; i++) {
            assert((receive_up[i].x >= start_proc_boundary) and (receive_up[i].x <= end_proc_boundary));
            local_particles.push_back(receive_up[i]);
        }
        for (int i = 0; i < status[1].MPI_TAG; i++) {
            assert((receive_down[i].x >= start_proc_boundary) and (receive_down[i].x <= end_proc_boundary));
            local_particles.push_back(receive_down[i]);
        }


        binBufferparticles(local_particles, grid, start_proc_boundary, end_proc_boundary, NumGrid, OneByDx, startN);

    }



    MPI_Finalize();
    return 0;
}

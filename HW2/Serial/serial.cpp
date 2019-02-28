#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"


unsigned int getGridIndex(const unsigned int i, const unsigned int j, const unsigned int NumGrid){
    // assert((i < NumGrid) and (i >= 0));
    // assert((j < NumGrid) and (j >= 0));
    return(i*NumGrid + j);
}
//
//  benchmarking program
//
int main(int argc, char **argv) {
    int navg, nabsavg = 0;
    double davg, dmin, absmin = 1.0, absavg = 0.0;

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

    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname ? fopen(sumname, "a") : NULL;

    particle_t *particles = (particle_t *) malloc(n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);

    //
    //  simulate a number of time steps
    //
    /** Setting up Grid parameters**/

    const unsigned int NumGrid = ceil(get_size()/0.02);
    const double dx = NumGrid/get_size();
    const int ParticlePerGrid = ceil(n/(NumGrid));

    grid_t * grid = (grid_t *) malloc(NumGrid* NumGrid*sizeof(grid_t));

    for(int i = 0;i < NumGrid*NumGrid; i++){
        grid[i].id_ = (int *) malloc(ParticlePerGrid*2* sizeof(int));
        /**Computing neigbourhood index of grid **/
        unsigned int idx = i/NumGrid;
        unsigned int idy = i%NumGrid;
        for(int j = 0; j < 8; j++){
            grid[i].neighbours_[j] = -1;
        }
        if((idx == 0) and (idy == 0)){
            grid[i].neighbours_[N]  = getGridIndex(idx  ,idy+1,NumGrid);
            grid[i].neighbours_[NE] = getGridIndex(idx+1,idy+1,NumGrid);
            grid[i].neighbours_[E]  = getGridIndex(idx+1,idy  ,NumGrid);
        }
        else if((idx == NumGrid - 1) and (idy == 0)){
            grid[i].neighbours_[W]  = getGridIndex(idx-1,idy  ,NumGrid);
            grid[i].neighbours_[NW] = getGridIndex(idx-1,idy+1,NumGrid);
            grid[i].neighbours_[N]  = getGridIndex(idx  ,idy+1,NumGrid);
        }
        else if((idx == 0) and (idy == NumGrid - 1)){
            grid[i].neighbours_[S ] = getGridIndex(idx  ,idy-1,NumGrid);
            grid[i].neighbours_[SE] = getGridIndex(idx+1,idy-1,NumGrid);
            grid[i].neighbours_[E]  = getGridIndex(idx+1,idy  ,NumGrid);
        }
        else if((idx == NumGrid - 1) and (idy == NumGrid - 1)){
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
        else if(idx == (NumGrid - 1)){
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



    double simulation_time = read_timer();

    for (int step = 0; step < NSTEPS; step++) {

        /**Reseting the counter of number of particles at every time steps**/
        for(int i = 0;i < NumGrid*NumGrid; i++){
            grid[i].count_ = 0;
        }

        /** Assigning particles onto grid **/
        int idx,idy;
        int id;
        for(int i = 0; i < n; i++){
            idx = static_cast<int>(particles[i].x * dx);
            idy = static_cast<int>(particles[i].y * dx);
            id  = (idx*NumGrid + idy);

            grid[id].id_[grid[id].count_] = i;
            grid[id].count_++;
            // assert(grid[id].count_ < ParticlePerGrid*2);
        }

        // /**Checking that no particle is lost**/
        // int sum = 0;
        // for(int i = 0; i < NumGrid*NumGrid;i++){
        //     sum += grid[i].count_;
        // }

//        for(int i = 0; i < NumGrid*NumGrid;i++){
//            printf("------------Location of the particle in Grid %d -----------------\n",i);
//            for(int j = 0; j < grid[i].count_; j++){
//                printf("%f %f\n",particles[grid[i].id_[j]].x,particles[grid[i].id_[j]].y);
//            }
//        }

        // if(sum != n){
        //     printf("Particles is lost during the simulation\n.Aborting\n");
        //     exit(0);
        // }
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        //
        //  compute forces
        //

       for(int i = 0; i < NumGrid*NumGrid; i++){
            for(int part = 0; part < grid[i].count_; part++){
                particles[grid[i].id_[part]].ax = 0;
                particles[grid[i].id_[part]].ay = 0;
                for(int partGrid = 0; partGrid < grid[i].count_; partGrid++){
                    apply_force(particles[grid[i].id_[part]],particles[grid[i].id_[partGrid]], &dmin, &davg, &navg);
                }
                for(int neigh = 0; neigh < 8; neigh++){
                    unsigned int neighbour_id =grid[i].neighbours_[neigh];
                    if( neighbour_id != -1){
                        for(int part_neigh = 0; part_neigh < grid[neighbour_id].count_; part_neigh++){
                            apply_force(particles[grid[i].id_[part]], particles[grid[neighbour_id].id_[part_neigh]], &dmin, &davg, &navg);
                        }
                    }
                }
            //    for (int j = 0; j < n; j++) {
            //        apply_force(particles[grid[i].id_[part]], particles[j], &dmin, &davg, &navg);
            //    }
            }
        }
    //    for (int i = 0; i < n; i++) {
    //        particles[i].ax = particles[i].ay = 0;
    //        for (int j = 0; j < n; j++) {
    //            apply_force(particles[i], particles[j], &dmin, &davg, &navg);
    //        }
    //    }

        //
        //  move particles
        //
        for (int i = 0; i < n; i++)
            move(particles[i]);

        if (find_option(argc, argv, "-no") == -1) {
            //
            // Computing statistical data
            //
            if (navg) {
                absavg += davg / navg;
                nabsavg++;
            }
            if (dmin < absmin) absmin = dmin;

            //
            //  save if necessary
            //
            if (fsave && (step % SAVEFREQ) == 0)
                save(fsave, n, particles);
        }
    }
    simulation_time = read_timer() - simulation_time;

    printf("n = %d, simulation time = %g seconds", n, simulation_time);

    if (find_option(argc, argv, "-no") == -1) {
        if (nabsavg) absavg /= nabsavg;
        //
        //  -the minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf(", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4) printf("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (absavg < 0.8) printf("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if (fsum)
        fprintf(fsum, "%d %g\n", n, simulation_time);

    //
    // Clearing space
    //
    if (fsum)
        fclose(fsum);
    free(particles);
    if (fsave)
        fclose(fsave);

    return 0;
}

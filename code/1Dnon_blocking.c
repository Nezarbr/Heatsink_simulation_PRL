#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>  
#include <mpi.h>
#include <sys/time.h>

/* AUTHOR Simulation : Charles Bouillaguet <charles.bouillaguet@lip6.fr>
   USAGE  : compile with -lm (and why not -O3)
            redirect the standard output to a text file
            gcc heatsink.c -O3 -lm -o heatsink
            ./heatsink > steady_state.txt
            then run the indicated python script for graphical rendering

   DISCLAIMER : this code does not claim to an absolute realism.
                this code could be obviously improved, but has been written so as
				to make as clear as possible the physics principle of the simulation.
*/

/* one can change the matter of the heatsink, its size, the power of the CPU, etc. */
#define ALUMINIUM
#define NORMAL       /* MEDIUM is faster, and FAST is even faster (for debugging) */
#define DUMP_STEADY_STATE

const double L = 0.15;      /* length (x) of the heatsink (m) */
const double l = 0.12;      /* height (z) of the heatsink (m) */
const double E = 0.008;     /* width (y) of the heatsink (m) */
const double watercooling_T = 20;   /* temperature of the fluid for water-cooling, (°C) */
const double CPU_TDP = 280; /* power dissipated by the CPU (W) */

/* dl: "spatial step" for simulation (m) */
/* dt: "time step" for simulation (s) */
#ifdef FAST
double dl = 0.004;
double dt = 0.004;
#endif

#ifdef MEDIUM
double dl = 0.002;
double dt = 0.002;
#endif

#ifdef NORMAL
double dl = 0.001;
double dt = 0.001;
#endif

#ifdef CHALLENGE
double dl = 0.0001;
double dt = 0.00001;
#endif

/* sink_heat_capacity: specific heat capacity of the heatsink (J / kg / K) */
/* sink_density: density of the heatsink (kg / m^3) */
/* sink_conductivity: thermal conductivity of the heatsink (W / m / K) */
/* euros_per_kg: price of the matter by kilogram */
#ifdef ALUMINIUM
double sink_heat_capacity = 897;
double sink_density = 2710;
double sink_conductivity = 237;
double euros_per_kg = 1.594;
#endif

#ifdef COPPER
double sink_heat_capacity = 385;
double sink_density = 8960;
double sink_conductivity = 390;
double euros_per_kg = 5.469;
#endif

#ifdef GOLD
double sink_heat_capacity = 128;
double sink_density = 19300;
double sink_conductivity = 317;
double euros_per_kg = 47000;
#endif

#ifdef IRON
double sink_heat_capacity = 444;
double sink_density = 7860;
double sink_conductivity = 80;
double euros_per_kg = 0.083;
#endif

const double Stefan_Boltzmann = 5.6703e-8;  /* (W / m^2 / K^4), radiation of black body */
const double heat_transfer_coefficient = 10;    /* coefficient of thermal convection (W / m^2 / K) */
double CPU_surface;


double wallclock_time()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6);
}

/* 
 * Return True if the CPU is in contact with the heatsink at the point (x,y).
 * This describes an AMD EPYC "Rome".
 */
static inline bool CPU_shape(double x, double y)
{
    x -= (L - 0.0754) / 2;
    y -= (l - 0.0585) / 2;
    bool small_y_ok = (y > 0.015 && y < 0.025) || (y > 0.0337 && y < 0.0437);
    bool small_x_ok = (x > 0.0113 && x < 0.0186) || (x > 0.0193 && x < 0.0266)
        || (x > 0.0485 && x < 0.0558) || (x > 0.0566 && x < 0.0639);
    bool big_ok = (x > 0.03 && x < 0.045 && y > 0.0155 && y < 0.0435);
    return big_ok || (small_x_ok && small_y_ok);
}

/* returns the total area of the surface of contact between CPU and heatsink (in m^2) */
double CPU_contact_surface()
{
    double S = 0;
    for (double x = dl / 2; x < L; x += dl)
        for (double y = dl / 2; y < l; y += dl)
            if (CPU_shape(x, y))
                S += dl * dl;
    return S;
}

/* Returns the new temperature of the cell (i, j, k). For this, there is an access to neighbor
 * cells (left, right, top, bottom, front, back), except if (i, j, k) is on the external surface. */
static inline double update_temperature(double *T, int u, int n, int m, int o, int i, int j, int k)
{
		/* quantity of thermal energy that must be brought to a cell to make it heat up by 1°C */
    const double cell_heat_capacity = sink_heat_capacity * sink_density * dl * dl * dl; /* J.K */
    const double dl2 = dl * dl;
    double thermal_flux = 0;

    if (i > 0)
        thermal_flux += (T[u - 1] - T[u]) * sink_conductivity * dl; /* neighbor x-1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (i < n - 1)
        thermal_flux += (T[u + 1] - T[u]) * sink_conductivity * dl; /* neighbor x+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (j > 0)
        thermal_flux += (T[u - n] - T[u]) * sink_conductivity * dl; /* neighbor y-1 */
    else {
        /* Bottom cell: does it receive it from the CPU ? */
        if (CPU_shape(i * dl, k * dl))
            thermal_flux += CPU_TDP / CPU_surface * dl2;
        else {
            thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
            thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
        }
    }

    if (j < m - 1)
        thermal_flux += (T[u + n] - T[u]) * sink_conductivity * dl; /* neighbor y+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (k > 0)
        thermal_flux += (T[u - n * m] - T[u]) * sink_conductivity * dl; /* neighbor z-1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (k < o - 1)
        thermal_flux += (T[u + n * m] - T[u]) * sink_conductivity * dl; /* neighbor z+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    /* adjust temperature depending on the heat flux */
    return T[u] + thermal_flux * dt / cell_heat_capacity;
}

/* Run the simulation on the k-th xy plane.
 * v is the index of the start of the k-th xy plane in the arrays T and R. */



int main()
{
    MPI_Init(NULL, NULL);

    int rank, size;
    	
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);



    CPU_surface = CPU_contact_surface();
    double V = L * l * E;
    int n = ceil(L / dl);
    int m = ceil(E / dl);
    int o = ceil(l / dl);

    if(rank ==0) {
    fprintf(stderr, "HEATSINK\n");
    fprintf(stderr, "\tDimension (cm) [x,y,z] = %.1f x %.1f x %.1f\n", 100 * L, 100 * E, 100 * l);
    fprintf(stderr, "\tVolume = %.1f cm^3\n", V * 1e6);
    fprintf(stderr, "\tWeight = %.2f kg\n", V * sink_density);
    fprintf(stderr, "\tPrice = %.2f €\n", V * sink_density * euros_per_kg);
    fprintf(stderr, "SIMULATION\n");
    fprintf(stderr, "\tGrid (x,y,z) = %d x %d x %d (%.1fMo)\n", n, m, o, 7.6293e-06 * n * m * o);
    fprintf(stderr, "\tdt = %gs\n", dt);
    fprintf(stderr, "CPU\n");
    fprintf(stderr, "\tPower = %.0fW\n", CPU_TDP);
    fprintf(stderr, "\tArea = %.1f cm^2\n", CPU_surface * 10000);

    }

    /* starting timer */
	double start = wallclock_time();

    double *T_final = malloc(n * m * o * sizeof(*T_final)); 
    
    double o_double = (double)o; // Conversion de o et m en double
 
    /*------------- Calcul du nombre de processus nécessaires -------------*/ 
    int kpl_per_process = ceil(o_double / size);
    int nbr_processk = ceil(o_double / kpl_per_process);

    /*------------- ------------- ------------- ------------- -------------*/

    // creating a subset communicator of ranks ranging from 0 to nbr_processk -1 
    // It is used in MPI_Gather as it only takes process data with same number of data


    // Check if the process should be part of the new subset communicator.
    int include_in_subset = (rank < nbr_processk ) ? 1 : MPI_UNDEFINED;

    // Create a subset of MPI_COMM_WORLD
    MPI_Comm subset_comm;
    MPI_Comm_split(MPI_COMM_WORLD, include_in_subset, rank, &subset_comm);

    /*------------- ------------- ------------- ------------- -------------*/

    //MPI_Status status;
  
    ////////////////////////////////////////////////////////


    if(rank< nbr_processk){
        double t = 0;
        int n_steps = 0;
        int convergence = 0;

    
        int kstart_pl = rank * kpl_per_process;
        int kend_pl = (rank == size - 1) ? o : kstart_pl + kpl_per_process;

        fprintf(stderr, "\t o= %d | kpl_per_process = %d| kstart_pl= %d | kend_pl= %d \n", o,kpl_per_process,kstart_pl,kend_pl);

     

        /* temperature of each cell, in degree Kelvin depends on decomposition  */
        double *T = malloc(n * m * (kend_pl-kstart_pl + 2)  * sizeof(*T)); // +2 to take into account the plan needed for edges (T for moment t)
        double *T_del = malloc(n * m * (kend_pl-kstart_pl + 2)  * sizeof(*T)); // (T_del for moment t+1)

        if (T == NULL || T_del == NULL) {
            perror("T or T_del could not be allocated");
            exit(1);
        }

        /* initially the heatsink is at the temperature of the water-cooling fluid */
        for (int u = 0; u < n * m * (kend_pl-kstart_pl+2 ) ; u++){
            T_del[u] = T[u] = watercooling_T + 273.15;
        }

        /* simulating time steps */
        while (convergence == 0) {
            /* fix the edges at this moment by doing send recv ops*/
            if(rank==0){


                MPI_Request request1, request2;
                MPI_Status status[2];

                MPI_Isend(T+(kend_pl-kstart_pl -1 + 1)*n*m , n*m, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &request1);
                MPI_Irecv( T+(1+kend_pl-kstart_pl)*n*m , n*m, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &request2);

                // Check the completion of non-blocking operations
                int flag1, flag2;
                MPI_Test(&request1, &flag1, &status[0]);
                MPI_Test(&request2, &flag2, &status[1]);

                // Continue with computation until both send and receive are completed
                while (!flag1 || !flag2) {
                    // Perform computation while waiting for communication to complete

                    // Check the completion of non-blocking operations in each iteration
                    MPI_Test(&request1, &flag1, &status[0]);
                    MPI_Test(&request2, &flag2, &status[1]);
                }
            
            
            }
            else if((rank>0)&(rank<size-1)){

                MPI_Request request, request2, request3, request4;
                MPI_Status status[2];

                // Non-blocking send and receive from the left neighbor
                MPI_Isend(T + n * m, n * m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request);
                MPI_Irecv(T, n * m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request2);

                // Non-blocking send and receive from the right neighbor
                MPI_Isend(T + (kend_pl - kstart_pl) * n * m, n * m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &request3);
                MPI_Irecv(T + (kend_pl - kstart_pl + 1) * n * m, n * m, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &request4);

                // Check the completion of non-blocking operations
                int flag1, flag2;
                MPI_Test(&request, &flag1, &status[0]);
                MPI_Test(&request2, &flag2, &status[1]);

                // Continue with computation while waiting for communication with the left neighbor to complete
                while (!flag1 || !flag2) {
                    // Perform computation while waiting for communication to complete

                    // Check the completion of non-blocking operations in each iteration
                    MPI_Test(&request, &flag1, &status[0]);
                    MPI_Test(&request2, &flag2, &status[1]);
                }

                // Reset flags for the next set of non-blocking operations
                flag1 = flag2 = 0;

                // Check the completion of non-blocking operations with the right neighbor
                MPI_Test(&request3, &flag1, &status[0]);
                MPI_Test(&request4, &flag2, &status[1]);

                // Continue with computation while waiting for communication with the right neighbor to complete
                while (!flag1 || !flag2) {
                    // Perform computation while waiting for communication to complete

                    // Check the completion of non-blocking operations in each iteration
                    MPI_Test(&request3, &flag1, &status[0]);
                    MPI_Test(&request4, &flag2, &status[1]);
                }

            }
            else{
                MPI_Request request3, request4;
                MPI_Status status[2];

                // Non-blocking send and receive from the left neighbor
                MPI_Isend(T + n * m, n * m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request3);
                MPI_Irecv(T, n * m, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request4);

                // Check the completion of non-blocking operations
                int flag1, flag2;
                MPI_Test(&request3, &flag1, &status[0]);
                MPI_Test(&request4, &flag2, &status[1]);

                // Continue with computation while waiting for communication with the left neighbor to complete
                while (!flag1 || !flag2) {
                    // Perform computation while waiting for communication to complete
                    // Check the completion of non-blocking operations in each iteration
                    MPI_Test(&request3, &flag1, &status[0]);
                    MPI_Test(&request4, &flag2, &status[1]);
                }

            }
                /* Update all cells. xy planes are processed, for increasing values of z. */
            for (int k = kstart_pl; k < kend_pl ; k++) {   // z 
                if(k==0){
                    k=k+1; // if k == 0 that plan is not processed
                }
                for (int j = 0; j < m; j++) {   // y
                    for (int i = 0; i < n; i++) {   // x
                        int u = (k-kstart_pl+1) * n * m + j * n + i;
                        T_del[u] = update_temperature(T, u, n, m, o, i, j, k);
                    }
                }
            }
            double delta_T = 0;
            double sum;

                /* each second, we test the convergence, and print a short progress report */
            if (n_steps % ((int)(1 / dt)) == 0) {
                double max = -INFINITY;
                for (int u = n * m  ; u < n * m * (kend_pl-kstart_pl); u++) {
                    delta_T += (T_del[u] - T[u]) * (T_del[u] - T[u]);
                    if (T_del[u] > max)
                        max = T_del[u];
                }

                MPI_Allreduce(&delta_T, &sum , 1, MPI_DOUBLE, MPI_SUM, subset_comm); // sum of delta_T on all processes and affects it to sum in each process
                
                sum = sqrt(sum) / dt;
                    
                fprintf(stderr, "t = %.1fs ; T_max = %.1f°C ; convergence = %g\n", t, max - 273.15, sum);
                fprintf(stderr, "delta = %g\n", sum);

                if (sum < 0.1)
                    convergence = 1;
            }

                /* the new temperatures are in T_del */
            double *tmp = T_del;
            T_del = T;
            T = tmp;
            t += dt;
            n_steps += 1;

                // Wait for the completion of all non-blocking operations.
            MPI_Barrier(subset_comm);
            }

        if(rank==0){
        // double *T_final = malloc(n * m * o * sizeof(*T_final));
        for (int k = 0; k <kend_pl ; k++) {   // z
                for (int j = 0; j < m; j++) {   // y
                    for (int i = 0; i < n; i++) {   // x
                        T_final[k*n*m+j*n+i] = T[k*n*m+j*n+i] ;
                    }
                }
            }
        }
        if (include_in_subset != MPI_UNDEFINED) {
            MPI_Gather(T+n*m, n*m*(kend_pl-kstart_pl), MPI_DOUBLE, T_final, n*m*(kend_pl-kstart_pl), MPI_DOUBLE, 0, subset_comm);      
        }

        if(rank== size-1){
            MPI_Send(T+n*m, n*m*(kend_pl-kstart_pl), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        double end = wallclock_time();
        fprintf(stderr, "Total computing time: %g sec\n", end - start);
        
        if(rank ==0){
            MPI_Recv(T_final+n*m*(size-1)*(kpl_per_process), n*m*(o-(size-1)*(kpl_per_process)), MPI_DOUBLE, size-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            #ifdef DUMP_STEADY_STATE
                printf("###### STEADY STATE; t = %.1f\n", t);
                for (int k = 0; k < o; k++) {   // z
                    printf("# z = %g\n", k * dl);
                    for (int j = 0; j < m; j++) {   // y
                        for (int i = 0; i < n; i++) {   // x
                            printf("%.1f ", T_final[k * n * m + j * n + i] - 273.15);
                        }
                        printf("\n");
                    }
                }
                printf("\n");
                fprintf(stderr, "For graphical rendering: python3 rendu_picture_steady.py [filename.txt] %d %d %d\n", n, m, o);
            #endif
                exit(EXIT_SUCCESS);
        }

  
    }

    free(T_final);

    MPI_Finalize();
}


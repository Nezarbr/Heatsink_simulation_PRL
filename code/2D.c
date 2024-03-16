#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>  
#include <mpi.h>
#include <sys/time.h>

/* AUTHOR simulation : Charles Bouillaguet <charles.bouillaguet@lip6.fr>
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
#define FAST     /* MEDIUM is faster, and FAST is even faster (for debugging) */
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
static inline double update_temperature(double *T, int u, int n, int m, int o, int i, int j, int k, int yst, int yend)
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
        thermal_flux += (T[u - n * (yend-yst+2)] - T[u]) * sink_conductivity * dl; /* neighbor z-1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (k < o - 1)
        thermal_flux += (T[u + n * (yend-yst+2)] - T[u]) * sink_conductivity * dl; /* neighbor z+1 */
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    /* adjust temperature depending on the heat flux */
    return T[u] + thermal_flux * dt / cell_heat_capacity;
}

/* Run the simulation on the k-th xy plane.
 * v is the index of the start of the k-th xy plane in the arrays T and R. */


//fonction qui retourne les diviseurs du milieu (deux fois le milieu si n est carré parfait || les deux du milieu sinon) et si n premier, effectue le traitement sur n-1
void divmid(int n, int* p, int* g) {
    int* divisors = malloc(n*sizeof(int)); 

    int count = 0; 

    for (int i = 1; i <= n; i++) {
        if (n % i == 0) {
            divisors[count] = i;
            count++;
        }
    }

    if(count == 2) {
        divmid(n-1, p, g);
        return;
    }
    if (count % 2 == 0) {
        *p = divisors[count / 2 - 1];
        *g = divisors[count / 2];
    } else {
        *p = *g = divisors[count / 2];
    }
    free(divisors);
}


int findproc(int proce_statek, int proce_statey, int sizesubset, int processy, int processk) {
    for(int i = 0; i < sizesubset; i++) {
        int row = i / processy;
        int col = i % processy;

        if ((row == proce_statek) && (col == proce_statey)) {
            return i;
        }
    }
    return -1;
}





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
    double *final = malloc(n * m * o * sizeof(*final)); 

    //////////////////////////////   case where sqrt(size) is an integer //////  //////  //////  //////  //////  ////// 

    /* let's go! we switch the CPU on and launch the simulation until it reaches a stationary state. */
    double t = 0;

    int n_steps = 0;
    int convergence = 0;

    ////////////////////////////////////////////////////////
    int grand_diviseur;
    int petit_diviseur;

    divmid(size, &petit_diviseur, &grand_diviseur);

    
    ////////////////////////////////////////////////////////
    /* Determine the number of processes needed; there may be processes that are not used */

    // Conversion de o et m en double
    double o_double = (double)o;
    double m_double = (double)m;

    // Calcul du nombre de processus nécessaires
    int kpl_per_process = ceil(o_double / grand_diviseur);
    int ypl_per_process = ceil(m_double / petit_diviseur);

    int nbr_processy = ceil(m_double / ypl_per_process);
    int nbr_processk = ceil(o_double / kpl_per_process);




    /*------------- ------------- ------------- ------------- -------------*/
    // creating a subset communicator of ranks ranging from 0 to size -1 
    // It is used in MPI_Gather as it only takes process data with same number of data


    // Check if the process should be part of the new subset communicator.
    int include_in_subset = (rank < nbr_processk*nbr_processy ) ? 1 : MPI_UNDEFINED;

    // Create a subset of MPI_COMM_WORLD
    MPI_Comm subset_comm;
    MPI_Comm_split(MPI_COMM_WORLD, include_in_subset, rank, &subset_comm);

    /*------------- ------------- ------------- ------------- -------------*/

  
    ////////////////////////////////////////////////////////

    // maj
    grand_diviseur = nbr_processk;
    petit_diviseur = nbr_processy;

    if(rank< nbr_processk*nbr_processy){

        // Allocate memory for gathered data
        int* recvcounts = (int*)malloc(nbr_processk*nbr_processy * sizeof(int));
        int* displs = (int*)malloc(nbr_processk*nbr_processy * sizeof(int));

        int* yend_lock = (int*)malloc(nbr_processk*nbr_processy * sizeof(int));
        int* ystart_lock  = (int*)malloc(nbr_processk*nbr_processy * sizeof(int));

        int* kend_lock = (int*)malloc(nbr_processk*nbr_processy * sizeof(int));
        int* kstart_lock  = (int*)malloc(nbr_processk*nbr_processy * sizeof(int));

        if (rank == 0) {
            

            // Calculate receive counts and displacements
            for (int i = 0; i < nbr_processk*nbr_processy; i++) {

                int kstart_pl = (int)(i / petit_diviseur) * kpl_per_process;
                int kend_pl = (i >= (nbr_processk*nbr_processy - petit_diviseur)) ? o : kstart_pl + kpl_per_process;

                int ystart_pl = (i % petit_diviseur) * ypl_per_process;
                int yend_pl =  ((i+1)%(petit_diviseur) == 0 ) ? m : ystart_pl + ypl_per_process;


                yend_lock[i] = yend_pl  ;
                ystart_lock[i] =  ystart_pl;
                kend_lock[i] =  kend_pl;
                kstart_lock[i] = kstart_pl;
                recvcounts[i] = (yend_pl-ystart_pl)*(kend_pl-kstart_pl)*n ;
                displs[i] = (i > 0) ? (displs[i - 1] + recvcounts[i - 1]) : 0;
                

            }
        }
        

        //////////////////////////////////////////////////////// partition sur k

        int kstart_pl = (int)(rank / petit_diviseur) * kpl_per_process;
        int kend_pl = (rank >= (nbr_processk*nbr_processy - petit_diviseur)) ? o : kstart_pl + kpl_per_process;

        //////////////////////////////////////////////////////// partition sur y 

        int ystart_pl = (rank % petit_diviseur) * ypl_per_process;
        int yend_pl =  ((rank+1)%(petit_diviseur) == 0 ) ? m : ystart_pl + ypl_per_process;

        fprintf(stderr,"\n rank is %d | kstart_pl is %d | kend_pl is %d  | kpl_per_process  is %d | ystart_pl is %d | yend_pl is %d  | ypl_per_process  is %d \n",rank ,kstart_pl ,kend_pl ,kpl_per_process ,ystart_pl ,yend_pl ,ypl_per_process );

        /* temperature of each cell, in degree Kelvin depends on decomposition  */
        double *T = malloc(n * (yend_pl-ystart_pl+2) * (kend_pl-kstart_pl + 2)  * sizeof(*T)); // +2 to take into account the plan needed for edges (T for moment t)
        double *T_del = malloc(n * (yend_pl-ystart_pl+2) * (kend_pl-kstart_pl + 2)  * sizeof(*T_del)); // (T_del for moment t+1)

        double *T_fin = malloc(n * (yend_pl-ystart_pl) * (kend_pl-kstart_pl )  * sizeof(*T_fin)); 

        double *T_cote1 = malloc(n * (kend_pl-kstart_pl)  * sizeof(*T_cote1));  // T_cote is the one that will be sent and received in following the y ax
        double *T_cote2 = malloc(n * (kend_pl-kstart_pl)  * sizeof(*T_cote2));  // T_cote is the one that will be sent and received in following the y ax

        if (T == NULL || T_del == NULL) {
            perror("T or T_del could not be allocated");
            exit(1);
        }

        /* initially the heatsink is at the temperature of the water-cooling fluid */
        for (int u = 0; u < n * (yend_pl-ystart_pl+2) * (kend_pl-kstart_pl+2 ) ; u++){
            T_del[u] = T[u] = watercooling_T + 273.15;
        }

       


        /*------------- ------------- ------------- ------------- -------------*/
        // creating a subset communicator of ranks ranging from 0 to size -1 
        // It is used in MPI_Gather as it only takes process data with same number of data

        /*------------- ------------- ------------- ------------- -------------*/

   
        /* simulating time steps */
        while (convergence == 0) {
            // process rank k communicates with k - 1 , k + 1, k + sqrt(size), k - sqrt(size) if they exist
            /* fix the edges at this moment by doing send recv ops*/
            if (rank % (petit_diviseur)==0){
                // T_cote takes elements from T_del to send 
                for(int lp = kstart_pl;lp<kend_pl;lp++){
                    for(int p=0;p<n;p++){
                        T_cote2[p+n*(lp-kstart_pl)] = T[ n*(yend_pl-ystart_pl) + n*(yend_pl-ystart_pl+2)*(lp-kstart_pl+1)+p ] ;
                    }
                }
                MPI_Sendrecv(T_cote2 , n*(kend_pl-kstart_pl), MPI_DOUBLE, rank+1, 0, T_cote1, n*(kend_pl-kstart_pl), MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
                // T_del takes elements from T_cotr to receive
                for(int lp = kstart_pl;lp<kend_pl;lp++){
                    for(int p=0;p<n;p++){
                        T[ n*(yend_pl-ystart_pl+1) + n*(yend_pl-ystart_pl+2)*(lp-kstart_pl+1)+p  ]   =  T_cote1[p+n*(lp-kstart_pl)]    ;        
                    }
                }
      
            }

            else if( ( (rank+1) % (petit_diviseur) )!=0){
     
                // T_cote takes elements from T_del to send 
                for(int lp = kstart_pl;lp<kend_pl;lp++){
                    for(int p=0;p<n;p++){
                        T_cote1[p+n*(lp-kstart_pl)] = T[ n + n*(yend_pl-ystart_pl+2)*(lp-kstart_pl+1)  + p ] ; 
                    }
                }
                MPI_Sendrecv(T_cote1 , n*(kend_pl-kstart_pl), MPI_DOUBLE, rank-1, 0, T_cote2, n*(kend_pl-kstart_pl), MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                // T_del takes elements from T_cotr to receive
                for(int lp = kstart_pl;lp<kend_pl;lp++){
                    for(int p=0;p<n;p++){

                        T[ p + n*(yend_pl-ystart_pl+2)*(lp-kstart_pl+1) ]=  T_cote2[ p+n*(lp-kstart_pl) ] ; 


                        T_cote1[p+n*(lp-kstart_pl)] = T[ n*(yend_pl-ystart_pl) + n*(yend_pl-ystart_pl+2)*(lp-kstart_pl+1)+p  ]; // what he will send
                    }
                }

                MPI_Sendrecv(T_cote1 , n*(kend_pl-kstart_pl), MPI_DOUBLE, rank+1, 0, T_cote2, n*(kend_pl-kstart_pl), MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for(int lp = kstart_pl;lp<kend_pl;lp++){
                    for(int p=0;p<n;p++){
                        T[ n*(yend_pl-ystart_pl+1) + n*(yend_pl-ystart_pl+2)*(lp-kstart_pl+1)+p   ]   =  T_cote2[p+n*(lp-kstart_pl)]   ;         // T_del + n*(yend_pl-ystart_pl + 1)*(l-kstart_pl+1)+p = T_cote2[p] ;
                    }
                }
            }

            else{
                // T_cote takes elements from T_del to send 
                for(int lp = kstart_pl;lp<kend_pl;lp++){
                    for(int p=0;p<n;p++){
                        T_cote1[p+n*(lp-kstart_pl)] = T[ p + n + n*(yend_pl-ystart_pl+2)*(lp-kstart_pl+1) ] ; 
                    }
                }
        
                MPI_Sendrecv(T_cote1 , n*(kend_pl-kstart_pl), MPI_DOUBLE, rank-1, 0, T_cote2, n*(kend_pl-kstart_pl), MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // T_del takes elements from T_cotr to receive
                for(int lp = kstart_pl;lp<kend_pl;lp++){
                    for(int p=0;p<n;p++){
                        T[ p + n*(yend_pl-ystart_pl+2)*(lp-kstart_pl+1) ]=  T_cote2[p+n*(lp-kstart_pl)] ;

                    }
                }
            }

            ///////////////////// send following k ax

            /* fix the edges at this moment by doing send recv ops*/
            if(rank<petit_diviseur){
                MPI_Sendrecv(T  + (kend_pl-kstart_pl -1 + 1)*n*(yend_pl-ystart_pl+2) , n*(yend_pl-ystart_pl+2), MPI_DOUBLE, rank+(petit_diviseur),
                0, T + (kend_pl-kstart_pl + 1)*n*(yend_pl-ystart_pl+2), n*(yend_pl-ystart_pl+2), MPI_DOUBLE, rank+(petit_diviseur), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else if(rank< (petit_diviseur*grand_diviseur - petit_diviseur)){
                MPI_Sendrecv(T  + n*(yend_pl-ystart_pl+2), n*(yend_pl-ystart_pl+2), MPI_DOUBLE, rank-(petit_diviseur), 0, T , 
                n*(yend_pl-ystart_pl+2), MPI_DOUBLE, rank-(petit_diviseur), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                MPI_Sendrecv(T  +(kend_pl-kstart_pl)*n*(yend_pl-ystart_pl+2) , n*(yend_pl-ystart_pl+2), MPI_DOUBLE, rank+(petit_diviseur), 0, T  +(kend_pl-kstart_pl+1)*n*(yend_pl-ystart_pl+2), 
                n*(yend_pl-ystart_pl+2), MPI_DOUBLE, rank+(petit_diviseur), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else{
      
                MPI_Sendrecv(T + n*(yend_pl-ystart_pl+2)  , n*(yend_pl-ystart_pl+2), MPI_DOUBLE, rank-(petit_diviseur), 0, T , n*(yend_pl-ystart_pl+2), 
                MPI_DOUBLE, rank-(petit_diviseur), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            }
            /* Update all cells. xy planes are processed, for increasing values of z. */
            for (int k = kstart_pl; k < kend_pl ; k++) {   // z 
                if(k==0){
                    k=k+1; // if k == 0 that plan is not processed
                }
                for (int j = ystart_pl; j < yend_pl; j++) {   // y
                    for (int i = 0; i < n; i++) {   // x

                        int u =  (k-kstart_pl+1) * n * (yend_pl-ystart_pl+2) + (j-ystart_pl+1) * n + i;    
                        
                        T_del[u] = update_temperature(T, u, n, m, o, i, j, k,ystart_pl,yend_pl);
                   
                        
                    }
                }
            }

            double delta_T = 0;
            double sum;

            /* each second, we test the convergence, and print a short progress report */
            if (n_steps % ((int)(1 / dt)) == 0) {
                double max = -INFINITY;

                /* Update all cells. xy planes are processed, for increasing values of z. */
                for (int k = kstart_pl; k < kend_pl ; k++) {   // z 

                    if(k==0){
                        k=k+1; // if k == 0 that plan is not processed
                    }
                    for (int j = ystart_pl; j < yend_pl; j++) {   // y
                        for (int i = 0; i < n; i++) {   // x

                            int u =  (k-kstart_pl+1) * n * (yend_pl-ystart_pl+2) + (j-ystart_pl+1) * n + i;    

                            delta_T += (T_del[u] - T[u]) * (T_del[u] - T[u]);
                            
                            
                            if (T_del[u] > max)
                                max = T_del[u];
                            
                        }
                    }

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


            // Attendre la fin de toutes les opérations non bloquantes
            MPI_Barrier(subset_comm);
        }

        /* Update all cells. xy planes are processed, for increasing values of z. */
            for (int k = kstart_pl+1; k < kend_pl+1 ; k++) {   // z 
                for (int j = ystart_pl+1; j < yend_pl+1; j++) {   // y
                    for (int i = 0; i < n; i++) {   // x  
                        T_fin[i+(j-ystart_pl-1)*n + (k-kstart_pl-1)*n*(yend_pl-ystart_pl)] = T_del[i + (j-ystart_pl)*n+ (k-kstart_pl)*n*(yend_pl-ystart_pl+2)];
                    }
                }
            }

        double end = wallclock_time();
        fprintf(stderr, "Total computing time: %g sec\n", end - start);

        // Use MPI_Gatherv to gather data with varying sizes
        MPI_Gatherv(T_fin, n * (yend_pl-ystart_pl) * (kend_pl-kstart_pl ),  MPI_DOUBLE, T_final, recvcounts, displs, MPI_DOUBLE, 0, subset_comm);

        if(rank ==0){

            for (int k = 0; k < o; k++) {   // z

            /// connaitre processes 
            int proce_state1 = (int)k/kpl_per_process; // 
                for (int j = 0; j < m; j++) {   // y
                    int proce_state2 = (int)j/ypl_per_process;
                    int result_rank = findproc(proce_state1, proce_state2, nbr_processk*nbr_processy, nbr_processy, nbr_processk);
                    for (int i = 0; i < n; i++) {   // x
                        final[i + n*j + n*m*k] =   T_final[displs[result_rank]+n*(j-ystart_lock[result_rank])+(k-kstart_lock[result_rank])*(yend_lock[result_rank]-ystart_lock[result_rank])*n + i];
                            
                    }
                        
                }
                    
            }

            #ifdef DUMP_STEADY_STATE
                printf("###### STEADY STATE; t = %.1f\n", t);
                for (int k = 0; k < o; k++) {   // z
                    printf("# z = %g\n", k * dl);
                    for (int j = 0; j < m; j++) {   // y
                        for (int i = 0; i < n; i++) {   // x
                            printf("%.1f ", final[k * n * m + j * n + i] - 273.15);
                           
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

    MPI_Finalize();
}

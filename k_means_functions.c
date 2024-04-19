#include "k_means_functions.h"



#define MAX 10000          // max value of your points set
#define MIN 0               // min value of your points set
#define DEV 1000            // deviation I will use to create the database starting from centroids

#define MAX_DIST sqrt(2)    // I will normalize the dataset, so the max distance between 2 different point
                            // will be 2^(1/2)

#define PRECISION 0.0000001  // It will be used to evaluate the precision of the founded centroids




Point* centroids_initialize(int k){
    /*
     * k = number of centroids you want to initialize
     *
     * DESCRIPTION:
     * It initializes a struct Point array of k dimension with normalized random point (in [0;1] interval)
     */

    Point *centroids = (Point*) malloc(k*sizeof (Point));

    for(int i = 0; i < k; i++){
        centroids[i].x = (double)((rand()%(MAX-MIN+1))+MIN)/(MAX-MIN);
        centroids[i].y = (double)((rand()%(MAX-MIN+1))+MIN)/(MAX-MIN);
    }

    return centroids;
}

Point* new_2D_database( int n, int k, Point *centroids) {
    /*
     * n = number of point in the database
     * k = number of centroids you want to use
     * centroids = struct Point array with centroids you want to use
     *
     * DESCRIPTION:
     * It takes the true centroids you want to use for your database and for each of them
     * it produces n/k point with deviation DEV (x_1 +/- DEV; y_1 +/- DEV).
     * To use the function you need OpenMP
     */

    Point *points = (Point*) malloc(n*sizeof (Point));
    if (points == NULL) return NULL;

    // if you have more centroids than point return NULL
    if (k > n){
        printf("ATTENTION! YOU HAVE MORE CENTROIDS THAN POINTS");
        return NULL;
    }

#pragma omp parallel default(none) shared(points, n, k) firstprivate(centroids)
    {

        int id = omp_get_thread_num();
        int tot_thr = omp_get_num_threads();
        int dim_c = (int) (n / k), dim_m = (int) (n / tot_thr);     // Divide points array in 'tot_thr' parts
        int start = id * dim_m, stop = (id + 1)*dim_m;              // Preparing start and stop for the next loop
        if(id == tot_thr - 1) stop = n;                             // For the last thread stop == n

        int n_chunk = (int) (start / dim_c);                        // First centroid every thread will use


        for (int i = start; i < stop; ++i ){
            points[i].x = centroids[n_chunk].x + (double)((rand()%(2*DEV+1))-DEV)/(MAX-MIN);
            points[i].y = centroids[n_chunk].y + (double)((rand()%(2*DEV+1))-DEV)/(MAX-MIN);

            // When n/k points are created for a single centroid, function passes to the next centroid
            if ((i+1)%dim_c == 0) {
                if (n_chunk + 1 < k) n_chunk++;
            }
        }
    }
    return points;
}

void find_2D_centroids(Point* db, Point* centroids, int n, int k, int num_thr) {
    /*
     * db = struct Point array of dimension n. It is the database the function will use
     * centroids = struct Point array of dimension k. These are the starting random centroids
     * n = number of points in database
     * k = number of centroids
     * num_thr = number of threads you want to activate
     *
     * DESCRIPTION:
     * It takes the database db and the starting random centroids, and it returns the calculated centroids.
     * To use the function you need OpenMP.
     */


#pragma omp parallel num_threads(num_thr) default(none) shared(n, num_thr, centroids, db, k)
    {
#pragma omp single
        {

            int counter, iterations = 0, my_exit = 0;       // Control variables.
            int single_id = omp_get_thread_num();           // id of the thread that will create tasks.

            int single_private_counters[num_thr][k];        // Matrix in which we will save the number of points associated
                                                            // with every centroid (for each thread).

            double single_private_square_dist[num_thr];     // Array in which we will save the sum of all the  square
                                                            // distances between a point and its centroid (for each thread).

            Point single_private_centroids[num_thr][k];     // Matrix in which we will calculate the new centroids
                                                            //  (for each thread).

            double max_dist = 2*pow(MAX_DIST,2);
            double precision = PRECISION;
            double old_mean_square_dist = max_dist;


            // First while: it runs until the stop condition is reached.
            while(my_exit !=1){

                iterations++;
                counter = 0;

                // Second while: the single_id thread starts to create tasks, one for each active thread.
                while (counter < num_thr) {

#pragma omp task default(none) shared(n, num_thr, k, single_private_counters, single_private_centroids,\
                                        single_private_square_dist, iterations) firstprivate(db, centroids, max_dist)

                    {
                        //I've divided the points array in n/num_thr parts, one for each thread.
                        int id = omp_get_thread_num();
                        int dim_m = (int) (n / num_thr);

                        //For each thread I find the start and the stop for the next loop.
                        int start = id * dim_m;
                        int stop = (id + 1) * dim_m;

                        //If the thread is the last one, its stop == n.
                        if (id == num_thr - 1) stop = n;

                        //Temporary variables
                        double d_min, d_temp_x, d_temp_y;
                        int num_centroid;

                        //Local thread variables. They will be used in the next loop.
                        Point thr_private_centr[k];
                        int thr_private_count[k];
                        double thr_private_square_dist = 0.0;


                        // Zero initialization for all the variables we will increase.
                        for (int i = 0; i < k; ++i) {
                            thr_private_count[i] = 0;
                            thr_private_centr[i].x = 0.0;
                            thr_private_centr[i].y = 0.0;
                        }

                        // Here I start to evaluate the cartesian distances.
                        for (int i = start; i < stop; ++i) {
                            d_min = max_dist;
                            //At the end of these second loop I will find the nearest centroid.
                            for (int j = 0; j < k; ++j) {
                                d_temp_x = (db[i].x - centroids[j].x) * (db[i].x - centroids[j].x);
                                d_temp_y = (db[i].y - centroids[j].y) * (db[i].y - centroids[j].y);

                                if (d_temp_x + d_temp_y < d_min) {
                                    d_min = d_temp_x + d_temp_y;
                                    num_centroid = j;                       // num_centroid is the position of the nearest centroid
                                                                            // in centroids array.
                                }
                            }

                            thr_private_square_dist += d_min;               // Increase the sum of the distances.
                            thr_private_count[num_centroid]++;              // Increase the number of point of the cluster.
                            thr_private_centr[num_centroid].x += db[i].x;   // Increase the sum of the x coordinates.
                            thr_private_centr[num_centroid].y += db[i].y;   // Increase the sum of the y coordinates.
                        }


                        // Here I save the private threads results in single_id thread variables.
                        single_private_square_dist[id] = thr_private_square_dist;
                        for (int i = 0; i < k; ++i) {
                            single_private_counters[id][i] = thr_private_count[i];
                            single_private_centroids[id][i].x = thr_private_centr[i].x;
                            single_private_centroids[id][i].y = thr_private_centr[i].y;
                        }
                    }

                    counter++;       //Ready for a new task.
                }

#pragma omp taskwait

                //I will use the counters of single_id thread to gather all the results.
                for (int i = 0; i < k; ++i) {

                    for (int j = 0; j < num_thr; ++j) {
                        if (j == single_id) continue;

                        // For all the other threads I increase the counters.
                        single_private_counters[single_id][i] += single_private_counters[j][i];
                        single_private_centroids[single_id][i].x += single_private_centroids[j][i].x;
                        single_private_centroids[single_id][i].y += single_private_centroids[j][i].y;

                        //For the first iteration, I sum all the square distances.
                        if (i == 0) single_private_square_dist[single_id] += single_private_square_dist[j];
                    }

                    //If a cluster is not empty I calculate the new centroids with the mass centre formula.
                    if (single_private_counters[single_id][i] != 0){
                        centroids[i].x = single_private_centroids[single_id][i].x/single_private_counters[single_id][i];
                        centroids[i].y = single_private_centroids[single_id][i].y/single_private_counters[single_id][i];
                    }

                    //If a cluster is empty, I randomly chose a new centroid.
                    else{
                        centroids[i].x = db[rand()%n+1].x;
                        centroids[i].y = db[rand()%n+1].y;
                    }
                }

                single_private_square_dist[single_id] = single_private_square_dist[single_id]/n;

                // - If you want to evaluate the performances, it is better to use a condition that forces the function
                //   to do the same number of cycles for every call. I've used:
                if (iterations == 40) my_exit = 1;
                // - If you don't want to evaluate the performances, you could use the precision condition:
                //if (fabs(old_mean_square_dist - single_private_square_dist[single_id]) < precision) my_exit = 1;

                old_mean_square_dist = single_private_square_dist[single_id];
            }
        }
    }
}

CentroidsVec *centroids_init_vec(){
    /*
     * DESCRIPTION:
     * It initializes a struct CentroidsVec array of N_CENTROIDS dimension with normalized random point (in [0;1] interval)
     */

    CentroidsVec *centroids;
    centroids = (CentroidsVec*) malloc(sizeof (CentroidsVec));

    for(int i = 0; i < N_CENTROIDS; i++){
        centroids->x[i] = (double)((rand()%(MAX-MIN+1))+MIN)/(MAX-MIN);
        centroids->y[i] = (double)((rand()%(MAX-MIN+1))+MIN)/(MAX-MIN);
    }

    return centroids;
}

CentroidsVec *random_centroids_selection(PointVec *db){
    /*
     * DESCRIPTION:
     * It initializes a struct CentroidsVec array of N_CENTROIDS dimension with normalized random point (in [0;1] interval)
     */

    CentroidsVec *centroids;
    centroids = (CentroidsVec*) malloc(sizeof (CentroidsVec));
    int temp;

    for(int i = 0; i < N_CENTROIDS; i++){
        temp = rand() % N;
        centroids->x[i] = db->x[temp];
        centroids->y[i] = db->y[temp];
        centroids->n_points[i] = 0;
    }

    return centroids;
}

PointVec *new_2D_database_vec(CentroidsVec *centroids) {
    /*
     * n = number of point in the database
     * centroids = struct CentroidsVec array with centroids you want to use
     *
     * DESCRIPTION:
     * It takes the true centroids you want to use for your database and for each of them
     * it produces n/k point with deviation DEV (x_1 +/- DEV; y_1 +/- DEV).
     * To use the function you need OpenMP
     */

    PointVec *points;
    points = (PointVec*) malloc(sizeof (PointVec));


#pragma omp parallel default(none) shared(points) firstprivate(centroids)
    {

        int id = omp_get_thread_num();
        int tot_thr = omp_get_num_threads();
        int dim_c = (int) (N / N_CENTROIDS), dim_m = (int) (N / tot_thr);     // Divide points array in 'tot_thr' parts
        int start = id * dim_m, stop = (id + 1)*dim_m;                        // Preparing start and stop for the next loop
        if(id == tot_thr - 1) stop = N;                                       // For the last thread stop == n
        int n_chunk = (int) (start / dim_c);                                  // First centroid every thread will use


        for (int i = start; i < stop; ++i ){
            points->x[i] = centroids->x[n_chunk] + (double)((rand()%(2*DEV+1))-DEV)/(MAX-MIN);
            points->y[i] = centroids->y[n_chunk] + (double)((rand()%(2*DEV+1))-DEV)/(MAX-MIN);

            // When n/k points are created for a single centroid, function passes to the next centroid
            if ((i+1)%dim_c == 0) {
                if (n_chunk + 1 < N_CENTROIDS) n_chunk++;
            }
        }
    }
    return points;
}







void k_meansFV(PointVec *db, CentroidsVec *centroids, int num_thr) {
    /*
     * db = struct PointVec array of dimension N. It is the database the function will use
     * centroids = struct CentroidsVec array of dimension N_CENTROIDS. These are the starting random centroids
     * num_thr = number of threads you want to activate
     *
     * DESCRIPTION:
     * It takes the database db and the starting random centroids, and it returns the calculated centroids.
     * To use the function you need OpenMP.
     */


#pragma omp parallel num_threads(num_thr) default(none) shared(num_thr, centroids, db)
    {
#pragma omp master
        {
            //int my_temp;
            int counter, iterations = 0, my_exit = 0;                   // Control variables.
            int master_id = omp_get_thread_num();                       // id of the thread that will create tasks.
            int dim_m = (int) (N / num_thr);                            // dimension of each block

            double master_private_square_dist[num_thr];                 // Array in which we will save the sum of all the  square
                                                                        // distances between a point and its centroid (for each thread).

            CentroidsVec master_private_centroids[num_thr];             // Matrix in which we will calculate the new centroids
                                                                        //  (for each thread).
            // Zero initialization
            for (int i = 0; i < num_thr; ++i) {
                for (int j = 0; j < N_CENTROIDS; ++j) {
                    master_private_centroids[i].x[j] = 0.0;
                    master_private_centroids[i].y[j] = 0.0;
                    master_private_centroids[i].n_points[j] = 0;
                }
            }

            double max_dist = 2 * pow(MAX_DIST,2);
            double precision = PRECISION;
            double old_mean_square_dist = max_dist;

            // First while: it runs until the stop condition is reached.
            while(my_exit !=1){

                iterations++;
                counter = 0;

                // Second while: the single_id thread starts to create tasks, one for each active thread.
                while (counter < num_thr) {

#pragma omp task default(none) shared(num_thr, master_private_centroids,\
                                        master_private_square_dist, iterations, dim_m) firstprivate(db, centroids, max_dist, counter)

                    {
                        //For each thread I find the start and the stop for the next loop.
                        int start = counter * dim_m;
                        int stop = (counter + 1) * dim_m;

                        //If the thread is the last one, its stop == N.
                        if (counter == num_thr - 1) stop = N;

                        //Temporary variables
                        double d_min, d_temp_x, d_temp_y;
                        int num_centroid;

                        //Local thread variables. They will be used in the next loop.
                        double thr_private_square_dist = 0.0;
                        CentroidsVec thr_private_centr;

                        // Zero initialization for all the variables we will increase.
                        for (int i = 0; i < N_CENTROIDS; ++i) {
                            thr_private_centr.n_points[i] = 0;
                            thr_private_centr.x[i] = 0.0;
                            thr_private_centr.y[i] = 0.0;
                        }

                        // Here I start to evaluate the cartesian distances.
                        for (int i = start; i < stop; ++i) {
                            d_min = max_dist;
                            //At the end of these second loop I will find the nearest centroid.
                            for (int j = 0; j < N_CENTROIDS; ++j) {
                                d_temp_x = (db->x[i] - centroids->x[j]) * (db->x[i] - centroids->x[j]);
                                d_temp_y = (db->y[i] - centroids->y[j]) * (db->y[i] - centroids->y[j]);

                                if (d_temp_x + d_temp_y < d_min) {
                                    d_min = d_temp_x + d_temp_y;
                                    num_centroid = j;                       // num_centroid is the position of the nearest centroid
                                                                            // in centroids array.
                                }
                            }

                            thr_private_square_dist += d_min;                // Increase the sum of the distances.
                            thr_private_centr.n_points[num_centroid]++;      // Increase the number of point of the cluster.
                            thr_private_centr.x[num_centroid] += db->x[i];   // Increase the sum of the x coordinates.
                            thr_private_centr.y[num_centroid] += db->y[i];   // Increase the sum of the y coordinates.
                        }


                        // Here I save the private threads results in master_id thread variables.
                        master_private_square_dist[counter] = thr_private_square_dist;

                        for (int i = 0; i < N_CENTROIDS; ++i) {
                            master_private_centroids[counter].n_points[i] = thr_private_centr.n_points[i];
                            master_private_centroids[counter].x[i] = thr_private_centr.x[i];
                            master_private_centroids[counter].y[i] = thr_private_centr.y[i];
                        }
                    }
                    counter++;
                    //Ready for a new task.
                }

#pragma omp taskwait

                //I will use the counters of master_id thread to gather all the results.
                for (int i = 0; i < num_thr; ++i) {
                    //printf("THR %d\t", i);

                    for (int j = 0; j < N_CENTROIDS; ++j) {

                        // For all the other threads I increase the counters.
                        if (i != master_id){
                            master_private_centroids[master_id].n_points[j] += master_private_centroids[i].n_points[j];
                            master_private_centroids[master_id].x[j] += master_private_centroids[i].x[j];
                            master_private_centroids[master_id].y[j] += master_private_centroids[i].y[j];
                        }

                        //For the first iteration, I sum all the square distances.
                        if (i == 0) master_private_square_dist[master_id] += master_private_square_dist[i];
                    }
                }

                for (int i = 0; i < N_CENTROIDS; ++i) {
                    //I calculate new centroids with the mass centre formula.
                    centroids->x[i] = master_private_centroids[master_id].x[i]/master_private_centroids[master_id].n_points[i];
                    centroids->y[i] = master_private_centroids[master_id].y[i]/master_private_centroids[master_id].n_points[i];

                    master_private_centroids[master_id].n_points[i] = 0;
                    master_private_centroids[master_id].x[i] = 0.0;
                    master_private_centroids[master_id].y[i] = 0.0;
                }

                master_private_square_dist[master_id] = master_private_square_dist[master_id]/N;

                // - If you want to evaluate the performances, it is better to use a condition that forces the function
                //   to do the same number of cycles for every call. I've used:
                if (iterations == 100) my_exit = 1;
                // - If you don't want to evaluate the performances, you could use the precision condition:
                //if (fabs(old_mean_square_dist - master_private_square_dist[master_id]) < precision) my_exit = 1;

                old_mean_square_dist = master_private_square_dist[master_id];
            }
        }
    }
}

void k_meansSV(PointVec *db, CentroidsVec *centroids, int num_thr){

    /*
     * db = struct PointVec array of dimension N. It is the database the function will use
     * centroids = struct CentroidsVec array of dimension N_CENTROIDS. These are the starting random centroids
     * num_thr = number of threads you want to activate
     *
     * DESCRIPTION:
     * It takes the database db and the starting random centroids, and it returns the calculated centroids.
     * To use the function you need OpenMP.
     */

    // Reduction variables
    double x[N_CENTROIDS], y[N_CENTROIDS];
    int n_points[N_CENTROIDS];


    double tempx, tempy;
    int n_closer_centroid, n_steps = 0;
    double max_dist;

    // Zero initialization
    for (int i = 0; i < N_CENTROIDS; ++i) {
        x[i] = 0.0; y[i] = 0.0; n_points[i] = 0;
    }

#pragma omp parallel num_threads(num_thr) default(none) shared(db, centroids, x, y, n_points)\
                                                            private(tempx, tempy, n_closer_centroid, max_dist, n_steps)
    {
        while (n_steps < 100){
            #pragma omp for reduction(+: y[:N_CENTROIDS], x[:N_CENTROIDS], n_points[:N_CENTROIDS])
            for (int i = 0; i < N; ++i) {
                max_dist = 2 * pow(MAX_DIST,2);
                for (int j = 0; j < N_CENTROIDS; ++j) {
                    tempx = (db->x[i] - centroids->x[j]) * (db->x[i] - centroids->x[j]);
                    tempy = (db->y[i] - centroids->y[j]) * (db->y[i] - centroids->y[j]);
                    if (tempx + tempy < max_dist){
                        n_closer_centroid = j;
                        max_dist = tempx + tempy;
                    }
                }
                x[n_closer_centroid] += db->x[i];
                y[n_closer_centroid] += db->y[i];
                n_points[n_closer_centroid] += 1;
            }
            #pragma omp barrier

            #pragma omp single
            {
                for (int i = 0; i < N_CENTROIDS; ++i) {
                    centroids->x[i] = x[i] / n_points[i];
                    centroids->y[i] = y[i] / n_points[i];
                    x[i] = 0.0; y[i] = 0.0; n_points[i] = 0;
                }
            }
            n_steps++;
        }
    }
}


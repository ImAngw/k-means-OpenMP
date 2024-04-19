#include "k_means_functions.h"



int main()
{
    printf("MAX THR %d\n", MAX_THR);
    printf("N PARTICLES %d\n", N);
    printf("N CENTROIDS %d\n\n", N_CENTROIDS);


    //#################################################################################################################################
    // For a single execution you can run this:

    srand(time(NULL));
    double start, stop, time;
    CentroidsVec *true_centroids = centroids_init_vec();                            // True centroids used to create the database
    PointVec *database;
    database = new_2D_database_vec(true_centroids);
    CentroidsVec *centroids= random_centroids_selection(database);              // First random initialization of the centroids



    // Sequential version SV, only one active thread
    start = omp_get_wtime();
    k_meansSV(database, centroids, 1);
    stop = omp_get_wtime();
    time = stop - start;
    printf("SEQUENTIAL TIME SV : %lf s\n", time);


    // Parallel version SV
    centroids = random_centroids_selection(database);
    start = omp_get_wtime();
    k_meansSV(database, centroids, MAX_THR);
    stop = omp_get_wtime();
    printf("PARALLEL TIME SV : %lf s\n", stop - start);
    printf("SPEEDUP %.2f\n", time / (stop - start));

    //graph_vec(database, centroids, true_centroids);


    printf("\n");


    // Sequential version FV, only one active thread
    centroids = random_centroids_selection(database);
    start = omp_get_wtime();
    k_meansFV(database, centroids, 1);
    stop = omp_get_wtime();
    printf("SEQUENTIAL TIME FV: %lf s\n", stop - start);
    time = stop - start;


    // Parallel version FV
    centroids = random_centroids_selection(database);
    start = omp_get_wtime();
    k_meansFV(database, centroids, MAX_THR);
    stop = omp_get_wtime();
    printf("PARALLEL TIME FV: %lf s\n", stop - start);
    printf("SPEEDUP %.2f\n", time / (stop - start));

    //graph_vec(database, centroids, true_centroids);


/*
    printf("FIRST VERSION\n");
    for (int j = 0; j < MAX_THR + 1; j++) {
        if ( j % 4 == 0){
            centroids = random_centroids_selection(database);
            start = omp_get_wtime();

            if (j == 0) k_meansFV(database, centroids, 1);
            else k_meansFV(database, centroids, j);

            stop = omp_get_wtime();
            printf("%f, ", stop - start);
        }
    }

    printf("\n");
    printf("SECOND VERSION\n");
    for (int j = 0; j < MAX_THR + 1; j++) {
        if ( j % 4 == 0){
            centroids = random_centroids_selection(database);
            start = omp_get_wtime();

            if (j == 0) k_meansSV(database, centroids, 1);
            else k_meansSV(database, centroids, j);

            stop = omp_get_wtime();
            printf("%f, ", stop - start);
        }
    }
*/


    free(true_centroids);
    free(database);
    free(centroids);

    return 0;
}

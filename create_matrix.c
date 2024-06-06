#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 4

int get_random_number() {
    return (rand() % 20) - 10;
}

void create_matrix(char* filename) {
    MPI_File fh;
    if (MPI_SUCCESS != MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh)) {
        printf(" ERROR : The file could not be opened. \n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int x;
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            int x = get_random_number();
            //x = 0;
            MPI_File_write(fh, &x, 1, MPI_INT, MPI_STATUS_IGNORE);
            printf(" %d", x);
        }
        printf("\n");
    }

    MPI_File_close(&fh);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    srand(time(NULL));

    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        char matrixA[] = "/home/varvara/Документы/skipod/matrixA.txt";
        printf("Filling in the matrix A\n");
        create_matrix(matrixA);
        printf("Finish\n");

        char matrixB[] = "/home/varvara/Документы/skipod/matrixB.txt";
        printf("Filling in the matrix B\n");
        create_matrix(matrixB);
        printf("Finish\n");
    }

    MPI_Finalize();

    return 0;
}
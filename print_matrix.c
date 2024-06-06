#include <stdio.h>
#include <mpi.h>

#define SIZE 4

void print_matrix(char* filename) {
    MPI_File fh;
    if (MPI_SUCCESS != MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh)) {
        printf(" ERROR : The file could not be opened. \n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            int x;
            MPI_File_read(fh, &x, 1, MPI_INT, MPI_STATUS_IGNORE);
            printf(" %d", x);
        }
        printf("\n");
    }
    MPI_File_close(&fh);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        char matrixA[] = "/home/varvara/Документы/skipod/matrixA.txt";
        char matrixB[] = "/home/varvara/Документы/skipod/matrixB.txt";
        char matrixC[] = "/home/varvara/Документы/skipod/matrixC.txt";

        printf("Matrix A:\n");
        print_matrix(matrixA);

        printf("Matrix B:\n");
        print_matrix(matrixB);

        printf("Matrix C = A*B:\n");
        print_matrix(matrixC);

    }

    MPI_Finalize();

    return 0;
}
#include <stdio.h>
#include <mpi.h>
#include <signal.h>
#include <mpi-ext.h>

#define SIZE 4
#define WIDTH 2

const char matrixA[] = "/home/varvara/Документы/skipod/matrixA.txt";
const char matrixB[] = "/home/varvara/Документы/skipod/matrixB.txt";
const char matrixC[] = "/home/varvara/Документы/skipod/matrixC.txt";

MPI_Comm comm;


MPI_File open_file(const char* filename, int amode) {
    MPI_File fh;
    if (MPI_SUCCESS != MPI_File_open(comm, filename, amode, MPI_INFO_NULL, &fh)) {
        printf(" ERROR : The file could not be opened. \n");
        fflush(stdout);
        MPI_Abort(comm, 1);
    }

    return fh;
}

MPI_Datatype get_lane_type(){
    int rank;
    MPI_Comm_rank(comm, &rank);

    MPI_Datatype lane_type;
    int c_sizes[2] = { SIZE, SIZE };
    int c_subsizes[2] = { WIDTH, SIZE };
    int c_starts[2] = { rank * WIDTH, 0 };
    MPI_Type_create_subarray(2, c_sizes, c_subsizes, c_starts, MPI_ORDER_C, MPI_INT, &lane_type);
    MPI_Type_commit(&lane_type);

    return lane_type;
}

MPI_Datatype get_column_type(){
    int rank;
    MPI_Comm_rank(comm, &rank);

    MPI_Datatype column_type;
    int l_sizes[2] = { SIZE, SIZE };
    int l_subsizes[2] = { SIZE, WIDTH };
    int l_starts[2] = { 0, 0 };
    MPI_Type_create_subarray(2, l_sizes, l_subsizes, l_starts, MPI_ORDER_C, MPI_INT, &column_type);
    MPI_Type_commit(&column_type);

    return column_type;
}

MPI_Datatype get_elem_type(){
    int rank;
    MPI_Comm_rank(comm, &rank);

    MPI_Datatype elem_type;
    int e_sizes[2] = { SIZE, SIZE };
    int e_subsizes[2] = { WIDTH, WIDTH };
    int e_starts[2] = { 0, 0 };
    MPI_Type_create_subarray(2, e_sizes, e_subsizes, e_starts, MPI_ORDER_C, MPI_INT, &elem_type);
    MPI_Type_commit(&elem_type);

    return elem_type;
}

void multiply_matrix() {
    int comm_size;
    int rank;
    int first = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Request request;
    int count = SIZE / WIDTH; // количество полос
    int lane[WIDTH][SIZE];
    int column[SIZE][WIDTH];
    MPI_Datatype elem_type = get_elem_type();
    MPI_Datatype lane_type = get_lane_type();
    MPI_Datatype column_type = get_column_type();

    MPI_File fhA = open_file(matrixA, MPI_MODE_RDONLY);
    MPI_File fhB = open_file(matrixB, MPI_MODE_RDONLY);
    MPI_File fhC = open_file(matrixC, MPI_MODE_WRONLY);

    MPI_Offset offset = (rank % count) * WIDTH * sizeof(int);
    MPI_File_set_view(fhA, 0, MPI_INT, lane_type, "native", MPI_INFO_NULL);
    MPI_File_set_view(fhB, offset, MPI_INT, column_type, "native", MPI_INFO_NULL);


    //считываем начальные полосу и столбец
    if (rank < count) {
        MPI_File_read(fhA, lane, SIZE * WIDTH, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_read(fhB, column, SIZE * WIDTH, MPI_INT, MPI_STATUS_IGNORE);
    }

    MPI_File_close(&fhB);
    MPI_File_close(&fhA);

    for (int n = 0; n < count; n++) {
        MPI_Offset offsetС = rank * WIDTH * SIZE * sizeof(int) + ((n + rank) % count) * WIDTH * sizeof(int);
        MPI_File_set_view(fhC, offsetС, MPI_INT, elem_type, "native", MPI_INFO_NULL);

        if (rank < count) {
            int dest = (rank + count - 1) % count;
            MPI_Isend(column, WIDTH * SIZE, MPI_INT, dest, 0, comm, &request);

            int elem[WIDTH][WIDTH];
            for (int i = 0; i < WIDTH; i++) {
                for (int j = 0; j < WIDTH; j++) {
                    elem[i][j] = 0;
                    for (int k = 0; k < SIZE; k++) {
                        elem[i][j] += lane[i][k] * column[k][j];
                    }
                }
            }
            if (n == 0 && rank == 0 && first == 1) {
                printf("Process %d failed\n", rank);
                fflush(stdout);
                raise(SIGKILL);
            }

            MPI_File_write(fhC, elem, WIDTH * WIDTH, MPI_INT, MPI_STATUS_IGNORE);

        }

        //точка сбора и основных пройессов, и запасных
        int err = MPI_Barrier(comm);
        if (err != MPI_SUCCESS) {
            first = 0;//чтоб умер только один процесс
            MPIX_Comm_shrink(comm, &comm);
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &comm_size);

            if(comm_size < count) {
                printf("There is nothing to replace the fallen percentage\n");
                fflush(stdout);
                MPI_Abort(comm, 1);
            }

            //считываем соответствующую строку из A и столбец из B
            MPI_Datatype lane_type = get_lane_type();
            MPI_Datatype column_type = get_column_type();

            fhA = open_file(matrixA, MPI_MODE_RDONLY);
            fhB = open_file(matrixB, MPI_MODE_RDONLY);

            MPI_File_set_view(fhA, 0, MPI_INT, lane_type, "native", MPI_INFO_NULL);
            MPI_Offset offsetB = ((n + rank) % count) * WIDTH * sizeof(int);
            MPI_File_set_view(fhB, offsetB, MPI_INT, column_type, "native", MPI_INFO_NULL);

            if (rank < count) {
                MPI_File_read(fhA, lane, SIZE * WIDTH, MPI_INT, MPI_STATUS_IGNORE);
                MPI_File_read(fhB, column, SIZE * WIDTH, MPI_INT, MPI_STATUS_IGNORE);
            }

            MPI_File_close(&fhB);
            MPI_File_close(&fhA);
            n--; //не переходим на следующую итерацию

        } else if (rank < count) {
            int source = (rank + 1) % count;
            MPI_Recv(column, WIDTH * SIZE, MPI_INT, source, 0, comm, MPI_STATUS_IGNORE);
        }
    }

    MPI_File_close(&fhC);
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int comm_size;
    if (comm_size < SIZE / WIDTH) {
        printf(" ERROR  \n");
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int rank;
    comm = MPI_COMM_WORLD;
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
    MPI_Comm_size(comm, &comm_size);
    MPI_Comm_rank(comm, &rank);

    printf("My rank is %d, %d processes in work\n", rank, comm_size);
    fflush(stdout);

    multiply_matrix();

    printf("Finish: my rank is %d \n", rank);
    fflush(stdout);


    MPI_Finalize();

    return 0;
}
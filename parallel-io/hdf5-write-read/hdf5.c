#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <hdf5.h>
#include <mpi.h>
#include <assert.h>

#define DATASIZE   64

void h5_writer(int, int *, int);
void h5_reader(int, int *, int);

int main(int argc, char *argv[])
{
    int my_id, ntasks, i, localsize;
    int *writevector, *readvector;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    if (ntasks > 64) {
        fprintf(stderr, "Datasize (64) should be divisible by number "
                "of tasks.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (DATASIZE % ntasks != 0) {
        fprintf(stderr, "Datasize (64) should be divisible by number "
                "of tasks.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    localsize = DATASIZE / ntasks;
    writevector = (int *) malloc(localsize * sizeof(int));
    readvector = (int *) malloc(localsize * sizeof(int));
    for (i = 0; i < localsize; i++)
        writevector[i] = i + 1 + localsize * my_id;

    h5_writer(my_id, writevector, localsize);

    h5_reader(my_id, readvector, localsize);

    for (i = 0; i < localsize; i++) 
	 assert(writevector[i]==readvector[i]);

    free(writevector);
    free(readvector);
    MPI_Finalize();
    return 0;
}

void h5_writer(int my_id, int *localvector, int localsize)
{
    herr_t status;
    hid_t plist_id, dset_id, filespace, memspace, file_id;
    hsize_t dims, counts, offsets;

    /* Create the handle for parallel file access property list
       and create a new file */

    // Create a new property list for file access
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    // Store MPI IO communicator info to the file access property list
    H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);

    // Create a new HDF5 file named "vector_<my_id>.h5"
    char* filename = malloc(sizeof(char)*10 + 10); // TODO: add better malloc logic, this should be good enough though.
    printf(filename, "vector_%i.h5",my_id);
    hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);

    /* Create the dataset */

    // Create a new simple dataspace for the file and open for access
    hid_t dataspace = H5Screate_simple(1, (const hsize_t[]){(hsize_t)localsize}, NULL);
    // Create a new dataset named "Vectors" for 'file'
    hid_t dataset = H5Dcreate(file, "Vectors", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Select a hyperslab of the file dataspace */

    // Number of blocks to be included in the hyperslab region
    hsize_t count[] = {1};
    // Select a hyperslab region of the file dataspace
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, (const hsize_t[]){(hsize_t)localsize}, NULL, count, NULL);

    /* Now we can write our local data to the correct position in the
       dataset. Here we use collective write, but independent writes are
       also possible. */

    // Create a new simple dataspace for the memory buffer and open for access
    hid_t memspace = H5Screate_simple(1, count, NULL);
    plist_id = H5Pcreate(H5P_DATASET_XFER);

    // Set data access to collective
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, plist_id, &localvector);

    /* Close all opened HDF5 handles */
    H5Pclose(plist_id);

}


void h5_reader(int my_id, int *localvector, int localsize)
{
    herr_t status;
    hid_t plist_id, dset_id, filespace, memspace, file_id;
    hsize_t dims, counts, offsets;

    /* Create the handle for parallel file access property list
       and open a file for reading */

    // Create a new property list for file access
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    // Store MPI IO communicator info to the file access property list
    H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);

    char* filename = malloc(sizeof(char)*10 + 10); // TODO: add better malloc logic, this should be good enough though.
    printf(filename, "vector_%i.h5",my_id);

    hid_t file = H5Dopen(filename, "Vectors", plist);

    /* Open the dataset and get the filespace id */

    /* Select a hyperslab of the file dataspace */

    /* Now we can read our local data from the correct position in the
       dataset. Here we use collective read but independent reads are
       also possible. */

    /* Close all opened HDF5 handles */
}

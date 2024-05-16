#pragma once

#include "mpi.h"

namespace mpi_util
{
    class MPIUtil
    {
    public:
        MPIUtil();
        ~MPIUtil();

        MPIUtil(int argc, char **argv);

        //void init(int argc, char **argv);
        void finalize();

        int get_rank();
        int get_size();
        void barrier();
        void abort(int error_code);
        void send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag);
        void recv(void *buf, int count, MPI_Datatype datatype, int source, int tag);
        void allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op);
        void gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root);
        void scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root);
        void allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype);


    private:
        MPI_Comm mpi_communicator;
        int mpi_rank;
        int mpi_size;
        int mpi_error_code;
        MPI_Status mpi_status;
        MPI_Request mpi_request;
    };

    const int n_mpi_processes(MPI_Comm mpi_communicator);
    const int this_mpi_process(MPI_Comm mpi_communicator);
}

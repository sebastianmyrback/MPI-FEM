#include "../include/utilities/mpi_util.hpp"

namespace mpi_util
{
    MPIUtil::MPIUtil()
    {
        mpi_error_code = MPI_Init(NULL, NULL);
        if (mpi_error_code != MPI_SUCCESS)
        {
            throw std::runtime_error("Error initializing MPI");
        }

        mpi_communicator = MPI_COMM_WORLD;
        MPI_Comm_rank(mpi_communicator, &mpi_rank);
        MPI_Comm_size(mpi_communicator, &mpi_size);
    }

    MPIUtil::~MPIUtil()
    {
        MPI_Finalize();
    }

    MPIUtil::MPIUtil(int argc, char **argv)
    {
        mpi_error_code = MPI_Init(&argc, &argv);
        if (mpi_error_code != MPI_SUCCESS)
        {
            throw std::runtime_error("Error initializing MPI");
        }

        // mpi_communicator = MPI_COMM_WORLD;
        // MPI_Comm_rank(mpi_communicator, &mpi_rank);
        // MPI_Comm_size(mpi_communicator, &mpi_size);
    }

    void MPIUtil::finalize()
    {
        MPI_Finalize();
    }

    int MPIUtil::get_rank()
    {
        return mpi_rank;
    }

    int MPIUtil::get_size()
    {
        return mpi_size;
    }

    void MPIUtil::barrier()
    {
        MPI_Barrier(mpi_communicator);
    }

    void MPIUtil::abort(int error_code)
    {
        MPI_Abort(mpi_communicator, error_code);
    }

    void MPIUtil::send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag)
    {
        MPI_Send(buf, count, datatype, dest, tag, mpi_communicator);
    }

    void MPIUtil::recv(void *buf, int count, MPI_Datatype datatype, int source, int tag)
    {
        MPI_Recv(buf, count, datatype, source, tag, mpi_communicator, &mpi_status);
    }

    void MPIUtil::allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op)
    {
        MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, mpi_communicator);
    }

    void MPIUtil::gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root)
    {
        MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, mpi_communicator);
    }

    void MPIUtil::scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root)
    {
        MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, mpi_communicator);
    }

    void MPIUtil::allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype)
    {
        MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, mpi_communicator);
    }


    const int n_mpi_processes(MPI_Comm mpi_communicator)
    {
        int n_mpi_processes;
        MPI_Comm_size(mpi_communicator, &n_mpi_processes);
        return n_mpi_processes;
    }

    const int this_mpi_process(MPI_Comm mpi_communicator)
    {
        int this_mpi_process;
        MPI_Comm_rank(mpi_communicator, &this_mpi_process);
        return this_mpi_process;
    }
}
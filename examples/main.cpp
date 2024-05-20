#include "../include/poisson.hpp"
#include "../include/utilities/mpi_util.hpp"

#include <chrono>

#define PARALLEL
//#define SERIAL

int main(int argc, char **argv)
{

    const double a = 0., b = 1.;
    //const int nintervals = 5e4;
    const int nintervals = 100;

    // Start timer
    #ifdef PARALLEL

    using namespace parallel_poisson;
    
    mpi_util::MPIUtil mpi_env(argc, argv);
    const int this_mpi_process = mpi_env.get_rank();
    const int n_mpi_processes  = mpi_env.get_size();

    double start = MPI_Wtime();

    Poisson1D poisson(a, b, nintervals);    

    double end_mesh = MPI_Wtime();

    poisson.setup_system();
    poisson.assemble_system();
    poisson.exchange_shared();
    
    double end_assemble = MPI_Wtime();
    
    poisson.solve();
    
    double end_solve = MPI_Wtime();
    
    poisson.output_solution("solution");
    
    double end = MPI_Wtime();

    double elapsed_time_total = end - start;
    double elapsed_time_mesh = end_mesh - start;
    double elapsed_time_assemble = end_assemble - end_mesh;
    double elapsed_time_solve = end_solve - end_assemble;
    double elapsed_time_output = end - end_solve;

    double max_time_mesh, max_time_assemble, max_time_solve, max_time_output, max_time;
    
    MPI_Reduce(&elapsed_time_mesh, &max_time_mesh, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time_assemble, &max_time_assemble, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time_solve, &max_time_solve, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time_output, &max_time_output, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time_total, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (this_mpi_process == 0) {
        std::cout << "Elapsed total time: " << max_time << " s" << std::endl;
        std::cout << "Mesh time: " << max_time_mesh << " s" << std::endl;
        std::cout << "Assemble time: " << max_time_assemble << " s" << std::endl;
        std::cout << "Solve time: " << max_time_solve << " s" << std::endl;
        std::cout << "Output time: " << max_time_output << " s" << std::endl;
    }
    #else
    using namespace serial_poisson;

    auto start = std::chrono::high_resolution_clock::now();
    
    Poisson1D poisson(a, b, nintervals);   

    auto end_mesh = std::chrono::high_resolution_clock::now();

    poisson.setup_system();
    poisson.assemble_system();

    auto end_assemble = std::chrono::high_resolution_clock::now();
    
    poisson.solve();

    auto end_solve = std::chrono::high_resolution_clock::now();
   
    poisson.output_solution("solution");

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Elapsed total time: " << std::chrono::duration<double>(end - start).count() << " s" << std::endl;
    std::cout << "Mesh time: " << std::chrono::duration<double>(end_mesh - start).count() << " s" << std::endl;
     std::cout << "Assemble time: " << std::chrono::duration<double>(end_assemble - end_mesh).count() << " s" << std::endl;
    std::cout << "Solve time: " << std::chrono::duration<double>(end_solve - end_assemble).count() << " s" << std::endl;
    std::cout << "Output time: " << std::chrono::duration<double>(end - end_solve).count() << " s" << std::endl;

    #endif

    return 0;
}

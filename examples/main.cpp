#include "../include/poisson.hpp"
#include "../include/utilities/mpi_util.hpp"

#include <chrono>

#define PARALLEL
//#define SERIAL

int main(int argc, char **argv)
{

    #ifdef PARALLEL
    using namespace parallel_poisson;
    mpi_util::MPIUtil mpi_env(argc, argv);
    const int this_mpi_process = mpi_env.get_rank();
    const int n_mpi_processes  = mpi_env.get_size();
    #else
    const int this_mpi_process = 0;
    using namespace serial_poisson;
    #endif

    const double a = 0., b = 1.;
    const int nintervals = 1e4;
    //const int nintervals = 10;

    // Start timer
    #ifdef PARALLEL
    double start = MPI_Wtime();

    Poisson1D poisson(a, b, nintervals);    
    
    poisson.setup_system();

    poisson.assemble_system();
        
    poisson.exchange_shared();
    
    double end_assemble = MPI_Wtime();
    
    const size_t n_iterations = poisson.solve();

    //if (this_mpi_process == 0)
      //std::cout << "Number of CG iterations: " << n_iterations << std::endl;
    
    double end_solve = MPI_Wtime();
    
    poisson.output_solution("solution");
    
    double end = MPI_Wtime();

    double elapsed_time_total = end - start;
    double elapsed_time_assemble = end_assemble - start;
    double elapsed_time_solve = end_solve - end_assemble;
    double elapsed_time_output = end - end_solve;

    double max_time_assemble, max_time_solve, max_time_output, max_time;
    
    MPI_Reduce(&elapsed_time_assemble, &max_time_assemble, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time_solve, &max_time_solve, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time_output, &max_time_output, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&elapsed_time_total, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (this_mpi_process == 0) {
        std::cout << "Elapsed total time: " << max_time << " s" << std::endl;
        std::cout << "Assemble time: " << max_time_assemble << " s" << std::endl;
        std::cout << "Solve time: " << max_time_solve << " s" << std::endl;
        std::cout << "Output time: " << max_time_output << " s" << std::endl;
    }
    #else
    auto start = std::chrono::high_resolution_clock::now();
    Poisson1D poisson(a, b, nintervals);    
   
    poisson.setup_system();
   
    poisson.assemble_system();

    auto end_assemble = std::chrono::high_resolution_clock::now();
    
    const size_t n_iterations = poisson.solve();
    //std::cout << "Number of CG iterations: " << n_iterations << std::endl;

    auto end_solve = std::chrono::high_resolution_clock::now();
   
    poisson.output_solution("solution");

    auto end = std::chrono::high_resolution_clock::now();

    // std::cout << "Elapsed total time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << std::endl;
    // std::cout << "Constructor time: " << std::chrono::duration_cast<std::chrono::seconds>(end_constructor - start).count() << " s" << std::endl;
    // std::cout << "Setup time: " << std::chrono::duration_cast<std::chrono::seconds>(end_setup - start_setup).count() << " s" << std::endl;
    // std::cout << "Assemble time: " << std::chrono::duration_cast<std::chrono::seconds>(end_assemble - start_assemble).count() << " s" << std::endl;
    // std::cout << "Solve time: " << std::chrono::duration_cast<std::chrono::seconds>(end_solve - start_solve).count() << " s" << std::endl;
    // std::cout << "Output time: " << std::chrono::duration_cast<std::chrono::seconds>(end_output - start_output).count() << " s" << std::endl;

    std::cout << "Elapsed total time: " << std::chrono::duration<double>(end - start).count() << " s" << std::endl;
     std::cout << "Assemble time: " << std::chrono::duration<double>(end_assemble - start).count() << " s" << std::endl;
    std::cout << "Solve time: " << std::chrono::duration<double>(end_solve - end_assemble).count() << " s" << std::endl;
    std::cout << "Output time: " << std::chrono::duration<double>(end - end_solve).count() << " s" << std::endl;

    #endif

    return 0;
}

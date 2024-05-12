#include <iostream>
#include "quadrature.hpp"
#include "mesh.hpp"
#include "fem.hpp"
#include "export.hpp"
#include "cg.hpp"
#include "norm.hpp"
#include <chrono>

namespace parallel_poisson 
{
    typedef std::map<std::pair<int, int>, double> SparseMatrix;
    typedef std::vector<double> Vector;

    class Poisson1D
    {
    public:
        
        Poisson1D();
        
        void run();

    private:

        void setup_system();
        void assemble_system();
        void solve();
        void compute_errors();
        void output_results() const;

        MPI_Comm mpi_communicator;

        const std::size_t n_mpi_processes;
        const std::size_t this_mpi_process;

        Mesh1D mesh;
        const quadrature::QuadratureRule<1> qr;
        const P1Lagrange1D<1> psi;

        const double f(const Point<1>);
        const double u(const Point<1>);

        SparseMatrix system_matrix;
        Vector system_rhs;
        Vector solution;

    };

    Poisson1D::Poisson1D()
        : mpi_communicator(MPI_COMM_WORLD),
          n_mpi_processes(MPI::n_mpi_processes(mpi_communicator)),
          this_mpi_process(MPI::this_mpi_process(mpi_communicator)),
          mesh(-1.5, 1.5, 10),
          qr(quadrature::midpoint)
    {}

    // define a function to be used as a source term f(x) = 8*pi^2*sin(2*pi*x)
    const double Poisson1D::f(const Point<1> & x) {
        //return 8.0 * M_PI * M_PI * sin(2.0 * M_PI * x[0]);
        return (M_PI*M_PI*(361*cos((19*M_PI*x[0])/10) - 441*cos((21*M_PI*x[0])/10)))/100;
    }

    const double Poisson1D::u(const Point<1> & x) {
        //return 2*sin(2.0 * M_PI * x[0]);
        return 2*sin(2*M_PI*x[0])*sin(M_PI*x[0]/10) + 10;
    }

    Poisson1D::setup_system()
    {
        mesh.partition(n_mpi_processes);

        // system_matrix.clear();
        // system_rhs.clear();
        // solution.clear();

        // const std::size_t n_dofs = mesh.n_cells() + 1;
        // system_matrix.resize(n_dofs, n_dofs, 0.0);
        // system_rhs.resize(n_dofs, 0.0);
        // solution.resize(n_dofs, 0.0);
    }

}   // namespace parallel_poisson




int main(int argc, char **argv)
{
    MPI::Environment env(argc, argv);

    parallel_poisson::Poisson1D poisson;
    poisson.run();

    return 0;
}

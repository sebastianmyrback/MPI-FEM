#ifndef MATVIRT_HPP
#define MATVIRT_HPP

#include "mesh.hpp"
#include "quadrature.hpp"
#include "PkLagrange.hpp"

bool dirichlet(const int );

/*
  Virtual Matrix
  We just need to define the matrix product operation
*/
template<class Mesh>
class  BaseMatrix {  
public:
    int n, m;                     // number of rows and columns
    std::shared_ptr<Mesh> Th;     // pointer to the matrix

    virtual double * mat_vec(double *x,double *ax) = 0; // the matrix vector product

    BaseMatrix(int);
    BaseMatrix(int, int, std::shared_ptr<Mesh>);
} ;

template<class Mesh>
MatVirt<Mesh>::MatVirt(int nn,int mm,Mesh*a)
{
  this->n=nn;
  this->m=mm;
  this->Th=a;  
}
  
template<class Mesh>
MatVirt<Mesh>::MatVirt(int n)
{
  this->n=n;
  this->m=n;
  this->Th=0;
}

/*
  Class for the production with stiffness matrix
*/
template<class Mesh>
class  MatStiff : public MatVirt<Mesh> {
public:

  MatStiff(int n,int m,Mesh* Th) : MatVirt<Mesh>(n,m,Th) {}
  double * Pmat(double *x,double *ax); 
};

/*
  Class for the production with mass matrix
*/
template<class Mesh>
class  MatMass : public MatVirt<Mesh> {
public:
  MatMass(int n,int m,Mesh* Th) : MatVirt<Mesh>(n,m,Th) { }
  double * Pmat(double *x,double *ax);  
};

/*
  Class for the production with identity matrix
*/
template<class Mesh>
class  MatId : public MatVirt<Mesh> {
public:
  MatId(int n) : MatVirt<Mesh>(n) {}
  double * Pmat(double *u,double *au){
    for( int i=0;i<this->n;++i) au[i] = u[i] ; return au;}      
};

double L2norm(Mesh1 *Th, double *u);
double H1norm(Mesh1 *Th, double *u);
double L2norm(Mesh2 *Th, double *u);
double H1norm(Mesh2 *Th, double *u);

template<int DIM>
void buildRHS(typename TypeMesh<DIM>::Mesh *Th, double *b, double (*f)(const double *),
	      const std::string & rule = "default") {

  typename LagrangePolynomial<DIM>::PLagrange Pk;          // the Lagrange Polynomials

  
  typedef typename TypeQuadrature<DIM>::Quadrature Quadrature;
  const Quadrature quadrature;
  const typename Quadrature::QF & qf = quadrature.setQuadrature(rule);


  const int n = Th->nv;
  for(int i=0; i< n; ++i) b[i]=0.0;                        // initialization

    
  const int nve = DIM+1;                                   // number of nodes in an element
  R phi[nve];

  for(int k=0; k< Th->nt; ++k)                             // loop over elements
    {
      const Element<DIM> & K(Th->t[k]);                    // the kth element

      int iK[nve];                                         // global indices of local nodes
      for(int i = 0; i < nve; ++i) iK[i] = K.v[i]->nu; 
      
      for(int ipq = 0; ipq < qf.n; ++ipq)  {               // loop over quadrature rule

	typename TypeQuadrature<DIM>::QP ip(qf[ipq]);      // the quadrature point
	double mip[DIM];
	K(ip, mip);
	
	
	Pk.shape(K, &ip.x, phi);                           // calculate phi(xq)
	const double Cint = K.mesure * ip.a;	           // integration constant
	
  	for( int i=0; i<nve; ++i) {                        // loop over local nodes

	  if( !dirichlet(K.v[i]->lab))                     // not derichlet
	    b[iK[i]] += Cint * f(mip) * phi[i];

	}
	  
      }      
    }  
}

#endif

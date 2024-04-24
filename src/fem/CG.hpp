#ifndef GC_HPP
#define GC_HPP
#include <assert.h>
#include <iostream>


static double* comblineaire(int n, double a, double *A,
                     double b, double *B,double *R )
{
    for(int i=0;i< n;++i)
        R[i] = a*A[i] + b*B[i];
    return R;
}
static double sdot(int n,double *A, double *B)
{
    double s=0;
    for(int i=0;i< n;++i)
        s += A[i]*B[i];
    return s;
}


template<class Mesh>
int GradienConjugue( MatVirt<Mesh> *A, // fonction et pointeur data pour A
		     MatVirt<Mesh> *C, // fonction et pointeur data pour C
		     double *b, // second membre
		     double * x, // solution qui contient une initialisation
		     int nbitermax,double eps) {

  int n = A->n;
  double * G  = new double[n];
  double * H  = new double[n];
  double * CG = new double[n];
  double * AH = CG;
  double * Ax = CG;
  double Cgg, Cggp, rho, gamma ;
  assert( n == A->m && n == C->n && n == C->m);
  comblineaire(n,1.,A->Pmat(x,Ax) , -1, b,G);
  C->Pmat(G,CG);
    
  for(int i=0;i< n;++i)  H[i] = - CG[i];
  Cgg= sdot(n,G,CG);
  if( Cgg > eps)
    {
      for(int iter=0; iter <nbitermax ; ++ iter)
        {
	  A->Pmat(H,AH);
	  rho = - sdot(n,G,H)/sdot(n,H,AH);          
	  comblineaire(n, 1.,x, rho, H,x);
	  comblineaire(n, 1.,G, rho, AH,G);
	  Cggp = Cgg;
	  C->Pmat(G,CG);
	  Cgg =sdot(n,G,CG);
	  gamma = Cgg / Cggp;
	  comblineaire(n, -1., CG, gamma, H, H);
	  if( Cgg < eps) {
	    printf(" - %d/%d  gg= %f %f rho=%f\n",iter,nbitermax,Cgg,gamma,rho);
	    break;
	  }	   
        }
    }
  
  delete [] G;
  delete [] H;
  delete [] CG;
    
  return Cgg < eps;
}

#endif

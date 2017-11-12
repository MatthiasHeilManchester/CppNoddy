/// \file EVP_real_sparse.cpp
/// \ingroup Examples
/// \ingroup EVP
/// Solves a 2x2 generalised eigenvalue problem
/// \f[ A_{2x2} \,{\underline x}_i = \lambda_i\, B_{2x2}\, {\underline x}_i \f]
/// for the 2 eigenvalues \f$ \lambda_i \f$, \f$i=1,2.\f$.
/// As a test case we use the SLEPc library.
/// In this case \f$ A_{2x2} \f$ and \f$ B_{2x2} = I \f$ are such that
/// the eigenvalues are \f$ 5,\, -1 \f$. The computation
/// requests eigenvalues that satisfy \f$ \vert\lambda\vert < 10\f$.

#include <cassert>
#include <iostream>

#include <Types.h>
#ifdef SLEPC
#include <SparseLinearEigenSystem.h>
#endif

using namespace CppNoddy;
using namespace std;

int main()
{
  cout << "\n";
  cout << "=== EVP: real eigenvalue problem                 ====\n";
  cout << "\n";

#ifndef PETSC_D
  cout << " PETSC-double support has not been included\n";
  cout << "\033[1;33;48m  * SKIPPED \033[0m\n";
#else
#ifndef SLEPC
  cout << " SCLEPC/PETSC support has not been included\n";
  cout << "\033[1;33;48m  * SKIPPED \033[0m\n";
#else

  SparseMatrix<double> a( 2, 2 );
  SparseMatrix<double> b( 2, 2 );
  a( 0, 0 ) = 1;
  a( 0, 1 ) = 2;
  a( 1, 0 ) = 4;
  a( 1, 1 ) = 3;
  b( 0, 0 ) = 1;
  b( 1, 1 ) = 1;

  // a vector for the eigenvalues
  DenseVector<D_complex> lambdas;
  // eigenvalues are: 5,-1

  SparseLinearEigenSystem<double> system( &a, &b );
  try
  {
    system.eigensolve();
  }
  catch ( std::runtime_error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    assert( false );
  }

  //// tag any eigenvalues within a distance of 2 of the point 4
  system.set_shift( 4 );
  system.tag_eigenvalues_disc( +1, 2 );
  // get those tagged eigenvalues
  lambdas = system.get_tagged_eigenvalues();
  const double tol = 1.e-13;
  //lambdas.dump();
  if ( std::abs( lambdas[ 0 ] - 5.0 ) > tol )
  {
    cout << "\033[1;31;48m  * FAILED \033[0m\n";
    cout.precision( 12 );
    cout << "    Final error = " << std::abs( lambdas[ 0 ] - 5.0 ) << "\n";
    lambdas.dump();
  }
  else
  {
    cout << "\033[1;32;48m  * PASSED \033[0m\n";
  }


#endif // SLEPC test
#endif // PETSC_D test
}
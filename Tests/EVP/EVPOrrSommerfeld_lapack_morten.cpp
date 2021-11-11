/// \file EVPOrrSommerfeld_lapack.cpp
/// \ingroup Tests
/// \ingroup EVP
/// Solves the following linear eigenvalue problem for values \f$ c \f$
/// that satisfy :
/// \f[ \phi''(y) - \alpha^2 \phi(y) - \psi(y) = 0\,, \f]
/// \f[ \psi''(y) - \alpha^2 \psi(y) - i \alpha Re \left \{ ( U(y) - c ) \psi(y) - U''(y) \phi \right \} = 0\,, \f]
/// subject to \f$ \phi(\pm 1) = \phi'(\pm 1) = 0 \f$ where
/// \f$ \alpha = 1.02 \f$, \f$ Re = 5772.2 \f$ and \f$ U(y) = 1 - y^2 \f$.
/// The matrix problem is constructed manually in this case, using second-order
/// finite differences.
/// These values approximately correspond to the first neutral temporal mode
/// in plane Poiseuille flow, therefore the test to be satisfied is that an eigenvalue
/// exists with \f$ \Im ( c ) \approx 0 \f$.

#include <EVP_bundle.h>

using namespace CppNoddy;
using namespace std;

// hacky include; stolen from oomph-lib; no include guards and stuff. 
#include "command_line_args.h"

// Other includes
#include <cfloat>    /* dbl_max */
#include <cmath>     /* erf */
  
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

namespace Parameters
{

 /// Which base flow profile are we using
 enum {Poiseuille,Cylinder_wake};

 /// Enumerated base flow
 unsigned Base_flow=Cylinder_wake; //Poiseuille;

 // Base wavenumber (to fit into 2D box for Jesper)
 double Base_wave_number_alpha=1.0;
 
 /// Number of waves in box
 unsigned N_waves_in_box=1;
 
 // Reynolds number
 double Reynolds_number=10.0; // 5772.2;
 
 /// Multiplier for exponential term in Cylinder wake profile
 // (max. velocity in free-stream (1.0) is depressed by that amount)
 double B_cylinder_wake=0.5;

 // Number of nodes in y direction
 unsigned Nnod_y=301;
 
 // Number of nodes in x direction (used for output only)
 unsigned Nnod_x=200;
 
 /// Height of lower boundary
 double Y_min=-5.0;

 /// Height of upper boundary
 double Y_max=5.0;
 
 /// compute: base flow profile u(y), its second deriv. d^2u/dy^2, the
 /// associated stream function psi and the vorticity omega.
 void base_flow(const double& y,
                double& u,
                double& u_dash_dash,
                double& psi,
                double& omega)
 {
  if (Base_flow==Poiseuille)
   {
    // Paranoia
    if ((Y_min!=-1.0)||(Y_max!=1.0))
     {
      cout << "Wrong geometry for Poiseuille flow\n";
      abort();
     }
    u = ( 1.0 - y * y );
    u_dash_dash = -2.0;
    psi=y-1.0/3.0*y*y*y+2.0/3.0; // constant chosen so that psi=0 at lower wall.
    omega=-2.0*y;
   }
  else if (Base_flow==Cylinder_wake)
   {
    // Veloc non-dim on max. velocity (outside wake) and length on "width"
    // of wake.
    u = 1.0-B_cylinder_wake*exp(-y*y);
    u_dash_dash = 2.0*B_cylinder_wake*exp(-y*y)*(1.0-2.0*y*y);
    psi=y-B_cylinder_wake/2.0*sqrt(4.0*atan(1.0))*erf(y);
    omega=-2.0*y*B_cylinder_wake*exp(-y*y);
   }
  else
   {
    cout << "Never get here. Base_flow = " << Base_flow << std::endl;
    abort();
   }
 }

}



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



int main(int argc, char* argv[])
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Number of nodes in y direction
 CommandLineArgs::specify_command_line_flag("--nnod_y",
                                            &Parameters::Nnod_y);

 // Number of nodes in x directio (for output only)
 CommandLineArgs::specify_command_line_flag("--nnod_x",
                                            &Parameters::Nnod_x);

 // Upper end of domain
 CommandLineArgs::specify_command_line_flag("--y_max",
                                            &Parameters::Y_max);
 
 // Lower end of domain
 CommandLineArgs::specify_command_line_flag("--y_min",
                                            &Parameters::Y_min);
 
 // Base wavenumber (to fit into box for Jesper)
 CommandLineArgs::specify_command_line_flag(
  "--base_wave_number_alpha",
  &Parameters::Base_wave_number_alpha);

 // Integer number of waves to fix into box
 CommandLineArgs::specify_command_line_flag("--n_waves_in_box",
                                            &Parameters::N_waves_in_box);

 // Reynolds number
 CommandLineArgs::specify_command_line_flag("--reynolds_number",
                                            &Parameters::Reynolds_number);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Set actual parameters
 double wave_number_alpha=Parameters::Base_wave_number_alpha*
  double(Parameters::N_waves_in_box);
 
 cout << "Wavenumber: " << wave_number_alpha << std::endl;
 

 //--------------------------------------
 // From here it's basically Rich's code
 //--------------------------------------

  // discretise with these many nodal points
 const std::size_t nodes( Parameters::Nnod_y);
  
  // we'll solve as TWO second order problems
  const std::size_t N( 2 * nodes );

  // domain boundaries
  const double left  = Parameters::Y_min;
  const double right = Parameters::Y_max;

  // spatial step for a uniform mesh
  const double d = ( right - left ) / ( nodes - 1 );

  // matrices for the EVP, initialised with zeroes
  DenseMatrix<D_complex> a( N, N, 0.0 );
  DenseMatrix<D_complex> b( N, N, 0.0 );


  // streamwise wavenumber and Reynolds number
  const double alpha ( wave_number_alpha *
                       double(Parameters::N_waves_in_box) );
  const double Re ( Parameters::Reynolds_number ); 
  const D_complex I( 0.0, 1.0 );

  // boundary conditions at the left boundary
  a( 0, 0 ) = 1.0;           // phi( left ) = 0
  a( 1, 0 ) = -1.5 / d;      // phi'( left ) = 0
  a( 1, 2 ) = 2.0 / d;
  a( 1, 4 ) = -0.5 / d;
  // fill the interior nodes
  for ( std::size_t i = 1; i <= nodes - 2; ++i )
  {
    // position in the channel
    const double y = left + i * d;

    // Base flow
    double U=0.0;
    double Udd=0.0;
    double psi_base=0.0;
    double omega_base=0.0;
    Parameters::base_flow(y,U,Udd,psi_base,omega_base);

    // the first quation at the i'th nodal point
    std::size_t row = 2 * i;
    a( row, row ) = -2.0 / ( d * d ) - alpha * alpha;
    a( row, row - 2 ) = 1.0 / ( d * d );
    a( row, row + 2 ) = 1.0 / ( d * d );
    a( row, row + 1 ) = -1.0;

    row += 1;
    // the second equation at the i'th nodal point
    a( row, row ) = -2.0 / ( d * d ) - alpha * alpha - I * alpha * Re * U;
    a( row, row - 2 ) = 1.0 / ( d * d );
    a( row, row + 2 ) = 1.0 / ( d * d );
    a( row, row - 1 ) = I * alpha * Re * Udd;

    b( row, row ) = - I * alpha * Re;
  }
  // boundary conditions at right boundary
  a( N - 2, N - 2 ) = 1.5 / d;
  a( N - 2, N - 4 ) = -2.0 / d;
  a( N - 2, N - 6 ) = 0.5 / d; // psi'( right ) = 0
  a( N - 1, N - 2 ) = 1.0;     // psi( right ) = 0
  // a vector for the eigenvalues
  DenseVector<D_complex> lambdas;
  DenseLinearEigenSystem<D_complex> system( &a, &b );

  try
  {
    system.eigensolve();
  }
  catch (const std::runtime_error &error )
  {
    cout << " \033[1;31;48m  * FAILED THROUGH EXCEPTION BEING RAISED \033[0m\n";
    return 1;
  }

  // This was Rich's setup for Poiseuille flow
  // // tag any eigenvalues with imaginary part > -0.1
  //system.set_shift( D_complex( 0.0, -0.1 ) );
  // system.tag_eigenvalues_upper( + 1 );

  // Get all eigenvalues with an imaginary part less than 0.1
  system.set_shift( D_complex( 0.0, 0.1 ) );
  system.tag_eigenvalues_lower( + 1 );
  lambdas = system.get_tagged_eigenvalues();

  cout << "Getting eigenvectors\n";
  DenseMatrix<D_complex> eig_vecs=system.get_tagged_eigenvectors();
  cout << "Done: Got " << eig_vecs.nrows() << " eigenvectors of length "
       << eig_vecs.ncols() << endl;

  std::string dirname("./DATA");
  mkdir( dirname.c_str(), S_IRWXU );

  TrackerFile spectrum( "./DATA/spectrum.dat" );
  spectrum.push_ptr( &lambdas, "evs" );
  spectrum.update();
 

 //--------------------------------------
 // End of Rich's code
 //--------------------------------------

  
  // Eigenvecs are rows of matrix -- how many do we have?
  unsigned n_ev=eig_vecs.nrows();

  // Do post-processing for the one with the smallest growth/decay rate
  unsigned i_ev=0; 
  double smallest_abs_imag=DBL_MAX;
  for (unsigned i=0;i<n_ev;i++)
   {
    if (abs(lambdas[i].imag())<smallest_abs_imag)
     {
      smallest_abs_imag=abs(lambdas[i].imag());
      i_ev=i;
      // cout << "New ev with smallest absolute imag part: " << i_ev << " "
      //      <<  smallest_abs_imag << std::endl;
     }
   }
  complex<double> omega_freq=lambdas[i_ev]*wave_number_alpha;
  cout << "Outputting fields for eigenvalue (c=omega/alpha): "
       << lambdas[i_ev] << " corresponding to omega = "
       << omega_freq << std::endl;

  // Normalise on peturbation to streamfct 
  double norm=0.0;
  for (unsigned i=0;i<Parameters::Nnod_y;i++)
   {
    norm+=
     pow(eig_vecs[i_ev][2*i].real(),2)+
     pow(eig_vecs[i_ev][2*i].imag(),2);
   }
  norm=sqrt(norm/double(Parameters::Nnod_y));
  for (unsigned i=0;i<2*Parameters::Nnod_y;i++)
   {
    eig_vecs[i_ev][i]/=norm;
   }


  // Output
  ofstream outfile;
  std::string filename=dirname+"/one_d_field.dat";
  outfile.open(filename.c_str());
  outfile << "# Eigenvalue (c=omega/alpha): " << lambdas[i_ev]
          << " ; omega = " << omega_freq << std::endl;
  for (unsigned i=0;i<Parameters::Nnod_y;i++)
   {
    const double y = left + double(i) * d;

    // Base flow
    double U=0.0;
    double Udd=0.0;
    double psi_base=0.0;
    double omega_base=0.0;
    Parameters::base_flow(y,U,Udd,psi_base,omega_base);


    // We're plotting quantities across width of domain:  y, phi(y) (complex),
    // vort(y) (complex), u_base, u_base'', phi_base, vort_base 
    outfile << y << " "
            << eig_vecs[i_ev][2*i].real() << " "
            << eig_vecs[i_ev][2*i].imag() << " "
            << -eig_vecs[i_ev][2*i+1].real() << " "
            << -eig_vecs[i_ev][2*i+1].imag() << " "
            << U << " "
            << Udd << " "
            << psi_base << " "
            << omega_base << " " 
            << std::endl;
   }
  outfile.close();
  
  
  // Output 2D streamfunction and vorticity at grid points (as initial
  // condition for Jesper)

  // Randomly do it at t=0
  double time=0.0;
  filename=dirname+"/two_d_field.dat";
  outfile.open(filename.c_str());
  outfile << "ZONE I=" << Parameters::Nnod_x << ", J=" << Parameters::Nnod_y 
          << std::endl; 
  for (unsigned i=0;i<Parameters::Nnod_y;i++)
   {
    const double y = left + double(i) * d;
    complex<double> phi=eig_vecs[i_ev][2*i];
    complex<double> psi=eig_vecs[i_ev][2*i+1];
    
    // Base flow
    double U=0.0;
    double Udd=0.0;
    double psi_base=0.0;
    double omega_base=0.0;
    Parameters::base_flow(y,U,Udd,psi_base,omega_base);
    Parameters::base_flow(y,U,Udd,psi_base,omega_base);
    
    for (unsigned j=0;j<Parameters::Nnod_x;j++)
     {
      const double x = double(j)/double(Parameters::Nnod_x-1)*
       (2.0*4.0*atan(1.0))/Parameters::Base_wave_number_alpha;
      
      // We're outputting: x, y, Phi_pert, Vort_pert, Phi_base, Vort_base
      outfile << x << " "
              << y << " "
              <<  (phi*exp(I*(wave_number_alpha*x-omega_freq*time))).real()
              << " " 
              << -(psi*exp(I*(wave_number_alpha*x-omega_freq*time))).real()
              << " "
              << psi_base << " "
              << omega_base << " " 
              << std::endl;
     }
   }
  outfile.close();

  
  return 0;
  
}

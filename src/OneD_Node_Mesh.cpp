/// \file OneD_Node_Mesh.cpp
/// Implementation of the one dimensional uniformly distributed mesh object.

#include <vector>
#include <string>

#include <DenseVector.h>
#include <OneD_Node_Mesh.h>
#include <Exceptions.h>

namespace CppNoddy
{

  template < typename _Type, typename _Xtype >
  void OneD_Node_Mesh<_Type, _Xtype>::set_nodes_vars( const std::size_t node, const DenseVector<_Type>& U )
  {
#ifdef PARANOID
    if ( U.size() != NV )
    {
      std::string problem;
      problem = " The OneD_Node_Mesh.set_node method is trying to add \n";
      problem += " an NVector of variables of a different size to that \n";
      problem += " stored at other nodal points. \n";
      throw ExceptionGeom( problem, NV, U.size() );
    }
#endif
    for ( std::size_t var = 0; var < U.size(); ++var )
    {
      VARS[ node * NV + var ] = U[ var ];
    }
  }

  template < typename _Type, typename _Xtype >
  DenseVector<_Type> OneD_Node_Mesh<_Type, _Xtype>::get_nodes_vars( const std::size_t &node ) const
  {
#ifdef PARANOID
    if ( ( node >= X.size() ) || ( node < 0 ) )
    {
      std::string problem;
      problem = " The OneD_Node_Mesh.get_nodes_vars method is trying to access \n";
      problem += " a nodal point outside of the range stored. \n";
      throw ExceptionRange( problem, X.size(), node );
    }
#endif
    DenseVector<_Type> nodes_vars;
    for ( std::size_t var = 0; var < NV; ++var )
    {
      nodes_vars.push_back( VARS[ node * NV + var ] );
    }
    return nodes_vars;
  }

  template < typename _Type, typename _Xtype >
  std::size_t OneD_Node_Mesh<_Type, _Xtype>::get_nnodes() const
  {
    return X.size();
  }

  template < typename _Type, typename _Xtype >
  std::size_t OneD_Node_Mesh<_Type, _Xtype>::get_nvars() const
  {
    return NV;
  }

  template < typename _Type, typename _Xtype >
  const DenseVector<_Xtype>& OneD_Node_Mesh<_Type, _Xtype>::nodes() const
  {
    return X;
  }

  template <>
  DenseVector<double> OneD_Node_Mesh<double, double>::find_roots1( const std::size_t &var, double value ) const
  {
    DenseVector<double> roots;
    for ( std::size_t node = 0; node < X.size() - 1; ++node )
    {
      std::size_t offset( node * NV + var );
      // find bracket nodes
      if ( ( VARS[ offset ] - value ) * ( VARS[ offset + NV ] - value ) < 0.0 )
      {
        double deriv = ( VARS[ offset + NV ] - VARS[ offset ] ) / ( X[ node + 1 ] - X[ node ] );
        double x = X[ node ] + ( value - VARS[ offset ] ) / deriv;
        // add the left hand node to the roots vector
        roots.push_back( x );
      }
    }
    return roots;
  }

  template < typename _Type, typename _Xtype >
  const DenseVector<_Type>& OneD_Node_Mesh<_Type, _Xtype>::vars_as_vector() const
  {
    return VARS;
  }

  template < typename _Type, typename _Xtype >
  void OneD_Node_Mesh<_Type, _Xtype>::set_vars_from_vector( const DenseVector<_Type>& vec )
  {
#ifdef PARANOID
    if  ( vec.size() != NV * X.size() )
    {
      std::string problem;
      problem = "The set_vars_from_vector method has been passed a vector\n";
      problem += "of a length that is of an incompatible size for this mesh object\n";
      throw ExceptionRuntime( problem );
    }
#endif
    VARS = vec;
  }

  template <>
  void OneD_Node_Mesh<double, double>::remesh1( const DenseVector<double>& newX )
  {
#ifdef PARANOID
    if ( std::abs( X[ 0 ] - newX[ 0 ] ) > 1.e-10 ||
         std::abs( X[ X.size() - 1 ] - newX[ newX.size() - 1 ] ) > 1.e-10 )
    {
      std::string problem;
      problem = " The OneD_Node_Mesh.remesh method has been called with \n";
      problem += " a passed coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw ExceptionRuntime( problem );
    }

    for ( std::size_t i = 0; i < newX.size() - 1; ++i )
    {
      if ( newX[ i ] >= newX[ i + 1 ] )
      {
        std::string problem;
        problem = " The OneD_Node_Mesh.remesh method has been passed \n";
        problem += " a non-monotonic coordinate vector. \n";
        throw ExceptionRuntime( problem );
      }
    }
#endif
    // copy current state of this mesh
    DenseVector<double> copy_of_vars( VARS );
    // resize the local storage
    VARS.resize( newX.size() * NV );

    // first nodal values are assumed to be untouched
    // loop thru destination mesh node at a time
    for ( std::size_t node = 1; node < newX.size() - 1; ++node )
    {
      // loop through the source mesh and find the bracket-nodes
      for ( std::size_t i = 0; i < X.size(); ++i )
      {
        if ( ( X[ i ] <= newX[ node ] ) && ( newX[ node ] < X[ i + 1 ] ) )
        {
          // linearly interpolate each variable in the mesh
          for ( std::size_t var = 0; var < NV; ++var )
          {
            double dX = newX[ node ] - X[ i ];
            double dvarsdX = ( copy_of_vars[ ( i+1 )*NV + var ] - copy_of_vars[ i*NV + var ] ) / ( X[ i + 1 ] - X[ i ] );
            VARS[ node * NV + var ] = copy_of_vars[ i * NV + var ] + dX * dvarsdX;
          }
        }
      }
    }

    // add the last nodal values to the resized vector
    for ( std::size_t var = 0; var < NV; ++var )
    {
      VARS[ ( newX.size() - 1 ) * NV + var ] = copy_of_vars[ ( X.size() - 1 ) * NV + var ];
    }
    // replace the old nodes with the new ones
    X = newX;
  }

  template <>
  void OneD_Node_Mesh<std::complex<double>, double>::remesh1( const DenseVector<double>& newX )
  {
#ifdef PARANOID
    if ( std::abs( X[ 0 ] - newX[ 0 ] ) > 1.e-10 ||
         std::abs( X[ X.size() - 1 ] - newX[ newX.size() - 1 ] ) > 1.e-10 )
    {
      std::string problem;
      problem = " The OneD_Node_Mesh.remesh method has been called with \n";
      problem += " a passed coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw ExceptionRuntime( problem );
    }

    for ( std::size_t i = 0; i < newX.size() - 1; ++i )
    {
      if ( newX[ i ] >= newX[ i + 1 ] )
      {
        std::string problem;
        problem = " The OneD_Node_Mesh.remesh method has been passed \n";
        problem += " a non-monotonic coordinate vector. \n";
        throw ExceptionRuntime( problem );
      }
    }
#endif
    // copy current state of this mesh
    DenseVector<std::complex<double> > copy_of_vars( VARS );
    // resize the local storage
    VARS.resize( newX.size() * NV );

    // first nodal values are assumed to be untouched
    // loop thru destination mesh node at a time
    for ( std::size_t node = 1; node < newX.size() - 1; ++node )
    {
      // loop through the source mesh and find the bracket-nodes
      for ( std::size_t i = 0; i < X.size(); ++i )
      {
        if ( ( X[ i ] <= newX[ node ] ) && ( newX[ node ] < X[ i + 1 ] ) )
        {
          // linearly interpolate each variable in the mesh
          for ( std::size_t var = 0; var < NV; ++var )
          {
            double dX = newX[ node ] - X[ i ];
            // if the paranoid checks above are satisfied, then the X[ i + 1 ] should still be in bounds
            std::complex<double> dvarsdX = ( copy_of_vars[ ( i+1 )*NV + var ] - copy_of_vars[ i*NV + var ] ) / ( X[ i + 1 ] - X[ i ] );
            VARS[ node * NV + var ] = copy_of_vars[ i * NV + var ] + dX * dvarsdX;
          }
        }
      }
    }

    // add the last nodal values to the resized vector
    for ( std::size_t var = 0; var < NV; ++var )
    {
      VARS[ ( newX.size() - 1 ) * NV + var ] = copy_of_vars[ ( X.size() - 1 ) * NV + var ];
    }
    // replace the old nodes with the new ones
    X = newX;

  }

  template <>
  void OneD_Node_Mesh<std::complex<double>, std::complex<double> >::remesh1( const DenseVector<std::complex<double> >& z )
  {
    std::string problem;
    problem = " The OneD_Node_Mesh.remesh method has been called with \n";
    problem += " a complex data set on a complex mesh.\n";
    throw ExceptionRuntime( problem );
  }

  template <>
  DenseVector<double> OneD_Node_Mesh<double, double>::get_interpolated_vars( const double& x_pos ) const
  {
    for ( unsigned node = 0; node < X.size() - 1; ++node )
    {
      // find bracketing nodes - incl shameless hack for evaluations at the boundary
      if ( ( X[ node ] < x_pos  || std::abs( X[ node ] - x_pos ) < 1.e-7 ) &&
           ( X[ node + 1 ] > x_pos || std::abs( X[ node + 1 ] - x_pos ) < 1.e-7 ) )
      {
        // distance from left node
        double delta_x( x_pos - X[ node ] );
        // empty data to return
        DenseVector<double> left;
        DenseVector<double> right;
        DenseVector<double> deriv;
        // interpolate data linearly
        left = get_nodes_vars( node );
        right = get_nodes_vars( node + 1 );
        deriv = ( right - left ) / ( X[ node + 1 ] - X[ node ] );
        // overwrite right
        right = left + deriv * delta_x;
        return right;
      }
    }
    std::cout << "You asked for a position of " << x_pos << " in a range " << X[ 0 ] << " to " << X[ X.size() - 1 ] << "\n";
    std::string problem;
    problem = "You have asked the OneD_Node_Mesh class to interpolate data at\n";
    problem += "a point that is outside the range covered by the mesh object.\n";
    throw ExceptionRuntime( problem );
  }

  template <>
  DenseVector<std::complex<double> > OneD_Node_Mesh<std::complex<double>, std::complex<double> >::get_interpolated_vars( const std::complex<double>& pos ) const
  {
    double x_pos( pos.real() );
#ifdef PARANOID
    std::cout << "WARNING: You are interpolating complex data on a complex mesh with 'get_interpolated_vars'.\n";
    std::cout << " This does a simple piecewise linear interpolating assuming a single valued path. \n";
#endif
    for ( unsigned node = 0; node < X.size() - 1; ++node )
    {
      // find bracketing nodes - incl shameless hack for evaluations at the boundary
      if ( ( X[ node ].real() < x_pos  || std::abs( X[ node ].real() - x_pos ) < 1.e-7 ) &&
           ( X[ node + 1 ].real() > x_pos || std::abs( X[ node + 1 ].real() - x_pos ) < 1.e-7 ) )
      {
        // distance from left node -- real coordinate is given. We also need to
        // interpolate between the two complex nodes -- hence imaginary coordinate is implict from
        // bracketing (complex) nodes
        std::complex<double> delta_z = ( X[ node + 1 ] - X[ node ] ) * ( x_pos - X[ node ].real() ) / ( X[ node + 1 ].real() - X[ node ].real() );
        // empty data to return
        DenseVector<std::complex<double> > left;
        DenseVector<std::complex<double> > right;
        DenseVector<std::complex<double> > deriv;
        // interpolate data linearly
        left = get_nodes_vars( node );
        right = get_nodes_vars( node + 1 );
        // derivative of the data
        deriv = ( right - left ) / ( X[ node + 1 ] - X[ node ] );
        // overwrite right
        right = left + deriv * delta_z;
        return right;
      }
    }
    std::cout << "You asked for a position of " << x_pos << " in a range " << X[ 0 ] << " to " << X[ X.size() - 1 ] << "\n";
    std::string problem;
    problem = "You have asked the OneD_Node_Mesh class to interpolate data at\n";
    problem += "a point that is outside the range covered by the mesh object.\n";
    problem += "Even for complex nodes we assume the path is single valued.\n";
    throw ExceptionRuntime( problem );
  }

  template <>
  DenseVector<double> OneD_Node_Mesh<std::complex<double>, double >::find_roots1( const std::size_t &var, double value ) const
  {
    std::string problem;
    problem = " The OneD_Node_Mesh.find_roots1 method has been called with \n";
    problem += " a mesh containing complex data.\n";
    throw ExceptionRuntime( problem );
  }

  template < typename _Type, typename _Xtype >
  _Type OneD_Node_Mesh<_Type, _Xtype>::integral2( std::size_t var ) const
  {
    _Type sum = 0.0;
    _Xtype dx = 0.0;
    // sum interior segments
    for ( std::size_t node = 0; node < X.size() - 1; ++node )
    {
      dx = ( X[ node + 1 ] - X[ node ] );
      sum += 0.5 * dx * ( VARS[ node * NV + var ] + VARS[ ( node+1 ) * NV + var ] );
    }
    // return the value
    return sum;
  }

  template < typename _Type, typename _Xtype >
  _Xtype OneD_Node_Mesh<_Type, _Xtype>::squared_integral2( std::size_t var ) const
  {
    _Xtype sum = 0.0;
    double dx = 0.0;
    // sum interior segments
    for ( std::size_t node = 0; node < X.size() - 1; ++node )
    {
      dx = std::abs( X[ node + 1 ] - X[ node ] );
      sum += 0.5 * dx * ( std::pow( std::abs( VARS[ node * NV + var ] ), 2 )
                          + std::pow( std::abs( VARS[ ( node + 1 ) * NV + var ] ), 2 ) );
    }
    // return the value
    return std::sqrt( sum );
  }

  template < typename _Type, typename _Xtype >
  _Type OneD_Node_Mesh<_Type, _Xtype>::integral4( std::size_t var ) const
  {
    if ( ( X.size() ) % 2 == 0 )
    {
      std::string problem;
      problem = " The OneD_Node_Mesh.Simpson_integral method is trying to run \n";
      problem += " on a mesh with an even number of points. \n";
      throw ExceptionRuntime( problem );
    }

    _Type f0, f1, f2;
    _Xtype x0, x1, x2;

    _Type sum = 0.0;

    // sum interior segments
    for ( std::size_t node = 0; node < X.size() - 2; node += 2 )
    {
      x0 = X[ node ];
      x1 = X[ node + 1 ];
      x2 = X[ node + 2 ];
      f0 = VARS[ node * NV + var ];
      f1 = VARS[ ( node+1 ) * NV + var ];
      f2 = VARS[ ( node+2 ) * NV + var ];
      sum += ( x2 - x0 )
             * (
               f1 * pow( x0 - x2, 2 ) + f0 * ( x1 - x2 ) * ( 2. * x0 - 3. * x1 + x2 )
               - f2 * ( x0 - x1 ) * ( x0 - 3. * x1 + 2. * x2 )
             )
             / ( 6. * ( x0 - x1 ) * ( x1 - x2 ) );
      // sum += (x1-x0)*( f0 + 4*f1 + f2 ) / 3.0 for equal spacing
    }
    // return the value
    return sum;
  }

  template < typename _Type, typename _Xtype >
  void OneD_Node_Mesh<_Type, _Xtype>::dump() const
  {
    std::cout << "Number of nodes = " << X.size() << "\n";
    std::cout << "Nodal positions :\n";
    X.dump();
    std::cout << "\n";
    std::cout << "Number of vars = " << NV << "\n";
    std::cout << "Interleaved mesh data : \n";
    VARS.dump();
    std::cout << "Mesh dump complete\n";
  }

  //the templated versions we require are:
  template class OneD_Node_Mesh<double>
  ;
  template class OneD_Node_Mesh<std::complex<double> >
  ;
  template class OneD_Node_Mesh<std::complex<double>, std::complex<double> >
  ;
}
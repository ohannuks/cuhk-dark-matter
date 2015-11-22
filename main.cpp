#include <gsl/gsl_errno.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <cstdio>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <cassert>
#include <fstream>
#include <algorithm>

#include <iomanip>      // std::setprecision
#include "rk.h"

#include <cstring>

typedef double Real;
using namespace std;

const Real ERROR = numeric_limits<Real>::quiet_NaN();
const bool check( const Real ERROR ) { 
  if( ERROR != ERROR) return false;
  return true;
}

//Constants
// c=1
const Real a = 20000; // pc
const Real M = 1e12; // Solar mass
const Real m = 4e6; // solar mass
const Real G = 4.79e-14; // Parcek per solar mass
const Real mp = 1;
const Real c = 1;
const Real Rs = 2.*G*m/pow(c,2); // pc, Schwarzschild radius

/*
// When c is set:
//const Real a = 20000;
//const Real M = 1e12;
//const Real m = 4e6;
//const Real G = 4.49e-15;
//const Real c = 0.3064;
//const Real mp = 1e-3;
//const Real Rs=2*G*m/pow(c,2);
*/

// Calculate radial action:
namespace Hernquist {
  ofstream output;

  //Solves limits of the equation
  void R_limits( const Real eps, const Real L, Real Rlimits[3] ) {
    // Check the input values
    {
      assert( eps > 0 );
      assert( abs(L) > 0 );
      assert( -20. + 1./eps - 8.*eps + pow(1. + 8.*eps,1.5)/eps >= 0 );
      assert( abs(L) < sqrt(-20. + 1./eps - 8.*eps + pow(1. + 8.*eps,1.5)/eps)/2. );
    }
    // eps = dimensionless newtonian energy
    // L = dimensionless angular momentum
    // Solve third order polynomial:
    double results[MAXDEGREE] = {0};
    double results_im[MAXDEGREE] = {0};
    double coeffs[MDP1];
    int degree = 3;
    
    // See hernquist.nb
    // We first solve limits for r
    // c[0]u^3+...
    //THESE COEFFICIENTS WERE CHECKED 22.11.
    coeffs[0] = 2.*eps; // 2 eps

    coeffs[1] = (-2.*a+2.*a*eps); // -2a+2a*eps

    coeffs[2] = pow(a*L,2); //a^2L^2

    coeffs[3] = a*pow(a*L,2); // a^3 L^2
    
    // Solve:
    rpoly_ak1(coeffs, &degree, results, results_im);
    sort( results, results + 3 );

    // Check results:
    for( int  i = 0; i < 3; ++i ) { 
      assert( results_im == 0 );
      if( results_im[i] != 0) { Rlimits[0]=ERROR; Rlimits[1]=ERROR; cout << "IMAGINARY VALUES" << endl; return;} 
    }

    assert( results[0] <= 0 );

    // Pick the limits; note R2<r<R3:
    Rlimits[0] = results[0];
    Rlimits[1] = results[1];
    Rlimits[2] = results[2];
    assert( Rlimits[2] >= Rlimits[1] );
    return;
  }
   
  Real integrand_function( Real r, void * params ){
    assert( r > 0 );
    // Get parameters
    const double * parameters;
    parameters = (double*) params;
    //These are inputted via the integration function (GSL)
    // These are now supposed to be dimensionless.
    const Real eps = parameters[0];
    const Real L = parameters[1];
    // Check the input values
    {
      assert( eps > 0 );
      assert( abs(L) > 0 );
      assert( -20. + 1./eps - 8.*eps + pow(1. + 8.*eps,1.5)/eps >= 0 );
      assert( abs(L) < sqrt(-20. + 1./eps - 8.*eps + pow(1. + 8.*eps,1.5)/eps)/2. );
    }

    // See Newtonian.nb
    const Real sqrtParameter = G*M*((-2*eps)/a - (a*pow(L,2))/pow(r,2) + 2/(a + r));
    // If imaginary, no contribution to the integral
    assert( sqrtParameter >= 0 );
    if( sqrtParameter < 0 ) { cerr << "Bad sqrtParameter "; cout << r << " " << eps << " " << L << " " << sqrtParameter << endl; return 0; }
    const Real fourVelocity = sqrt(sqrtParameter);
    const Real integrand = fourVelocity;

    assert( integrand >= 0 );
    
    //TODO: Check whether integrand is velocity or momentum!
    return integrand;
  }

  Real II_radial( const Real eps, const Real newtonian_dimensionless_L ) {
    // Note: eps = dimensionless hernquist energy per unit mass
    // And: L = dimensionless angular momentum per unit mass
    assert( eps>0 );
    Real Rlimits[3]={ERROR};
    R_limits(eps, newtonian_dimensionless_L, Rlimits); // Get limits

    assert( check(Rlimits[0]) );

    Real integral = 0;
    
    // Integrate with a singularity:
    {
      // <GSL>
      gsl_integration_workspace * work_space = 
                            gsl_integration_workspace_alloc(1000);
     
      double abs_error = 0;
      double rel_error = 1.e-6;
      size_t limit = 1000;
      double result, error;
     
      gsl_function integrand;
     
      double parameters[2] = {eps,newtonian_dimensionless_L};
      integrand.function = &integrand_function;
      integrand.params   = &parameters[0];
     
      // Please refer to hernquist.nb and Sadeghian article for these
      const Real Rmin = Rlimits[1];
      const Real Rmax = Rlimits[2];
      assert( Rmin>0 && Rmax>0 );
      assert( Rmin < Rmax );
      
      gsl_integration_qags( &integrand, Rmin, Rmax,
                            abs_error, rel_error, limit, work_space,
                            &result, &error );

      assert( abs(result) > 1e5*abs(error) );
     
      gsl_integration_workspace_free( work_space ); // Free memory
      // </GSL>
      integral = result;
    }

    return integral;
  }
}

// Calculate radial action:
namespace BlackHole {
  ofstream output;

  // Define dimensionless energy and angular momentum per particle mass with c=1
  // Input is always dimensionless eps!
  // COMMENT: dont change these, they've been changed way too many times already
  Real epsmax( const Real r, const Real eps ) {
    assert( r > 2.*G*m );
    Real eps_dimensional = ERROR;
    if( r >= 6.*G*m ) {
      eps_dimensional = (1. + (2.*G*m)/r)/sqrt(1. + (6.*G*m)/r);
    } else if( r >= 4.*G*m ) {
      eps_dimensional = (1. - (2.*G*m)/r)/sqrt(1. - (3.*G*m)/r);
    } else {
      assert( 1==0 );
    }
    return a/(G*M) * (1.-eps_dimensional);
  }
  
  Real epsmin( const Real r, const Real eps ) { 
    assert( r > 2.*G*m );
    // Dimensionless
    return 0;
  }
  
  Real Lmin( const Real r, const Real eps ) {
    assert( eps >= 0 && eps >= epsmin(r,eps) );
    assert( eps <= 1 );
    assert( r > 2.*G*m );
    // Note: we use dimensionless eps;
    const Real dimensional_eps = 1.-(G*M/a) * eps;
    Real L_dimensional = ERROR;
    L_dimensional = 4.*sqrt(2.)*sqrt((pow(G,2)*pow(m,2))/
     (-8. + 36.*pow(dimensional_eps,2) - 27.*pow(dimensional_eps,4) + 
       dimensional_eps*pow(-8. + 9.*pow(dimensional_eps,2),1.5)));
    return L_dimensional / sqrt(a*G*M);
  }
  
  Real Lmax( const Real r, const Real eps ) { 
    assert( eps >= 0 && eps >= epsmin(r,eps) );
    assert( eps <= 1 );
    assert( r > 2.*G*m );
    assert( -1. + eps/(1. - (2.*G*m)/r) > 0 );
    // Note: we use dimensionless eps;
    const Real dimensional_eps = 1.-(G*M/a) * eps;
    Real dimensional_L = sqrt(-1. + dimensional_eps/(1. - (2.*G*m)/r))*r;
    return dimensional_L / sqrt(a*G*M);
  }
  


  void u_limits( const Real eps, const Real L, Real ulimits[2] ) {
    //Note: eps, L are dimensionless
    const Real dimensional_eps = 1.-(G*M/a) * eps;
    const Real dimensional_L = L*sqrt(a*G*M);
    // Check ulimits
    {
    assert( dimensional_eps < 1 );
    assert( 1 < 3.*dimensional_eps/(2.*sqrt(2.)) );
    // Better be safe than sorry, see ulimits appendix
    assert( (sqrt((8. - 36.*pow(dimensional_eps,2) + 27.*pow(dimensional_eps,4) - 
          dimensional_eps*sqrt(pow(-8. + 9.*pow(dimensional_eps,2),3.)))/
        (-1. + pow(dimensional_eps,2)))*Rs)/(2.*sqrt(2.)) > abs(dimensional_L) );
    }

    // Solve third order polynomial:
    double results[MAXDEGREE] = {0};
    double results_im[MAXDEGREE] = {0};
    double coeffs[MDP1];
    int degree = 3;
    
    // See FourVelocity.nb
    // We first solve limits for u
    // c[0]u^3+...
    coeffs[0] = pow(dimensional_L,2)*Rs;

    coeffs[1] = -1.*pow(dimensional_L,2);

    coeffs[2] = Rs;

    coeffs[3] = -1.+pow(dimensional_eps,2);
    
    // Solve:
    rpoly_ak1(coeffs, &degree, results, results_im);
    // Sort
    sort(results, results+3);
    // Check results:
    for( int  i = 0; i < 3; ++i ) { 
      assert( results_im[i] == 0 );
      if( results_im[i] != 0) { ulimits[0]=ERROR; ulimits[1]=ERROR; return;} 
    }
    assert( results[0] <= 0 );
    if( results[0] <= 0 ) { ulimits[0]=ERROR; ulimits[1]=ERROR; cout << "bad data" << endl; return;}
    if( results[0] > results[1]){ ulimits[0]=ERROR; ulimits[1]=ERROR; cout << "errordata" << endl; return; }
    ulimits[0] = results[0];
    ulimits[1] = results[1];
    assert( ulimits[1] >= ulimits[0] );
    return;
  }
  
  Real integrand_function( Real u, void * params ){
    assert( u > 0 );
    // Get parameters
    const double * parameters;
    parameters = (double*) params;
    // Dimensionless!
    const Real eps = parameters[0];
    const Real L = parameters[1];
    // Dimensional! with mp=1, c=1
    const Real dimensional_eps = 1.-(G*M/a) * eps;
    const Real dimensional_L = L*sqrt(a*G*M);
    {
    // Check ulimits
    assert( dimensional_eps < 1 );
    assert( 1 < 3.*dimensional_eps/(2.*sqrt(2.)) );
    // Better be safe than sorry, see ulimits appendix
    assert( (sqrt((8. - 36.*pow(dimensional_eps,2) + 27.*pow(dimensional_eps,4) - 
          dimensional_eps*sqrt(pow(-8. + 9.*pow(dimensional_eps,2),3.)))/
        (-1. + pow(dimensional_eps,2)))*Rs)/(2.*sqrt(2.)) > abs(dimensional_L) );
    }
    
    // See Schwarzschild.nb
    const Real sqrtParameter = pow(dimensional_eps,2) - (1 - Rs*u)*(1 + pow(dimensional_L,2)*pow(u,2));
    assert( sqrtParameter >= 0 );
    const Real fourVelocity = sqrt(sqrtParameter);
    // Note: There are two minus signs which cancel each other due to: change of integration limtis u2, u1 to u1, u2 and dr/du
    const Real integrand = fourVelocity / pow(u,2);
    assert( integrand >= 0 );
    
    //TODO: Check the division!
    return integrand;
  }
  
  Real II_radial( const Real eps, const Real L ) {
    // Check input
    {
      const Real dimensional_eps = 1.-(G*M/a) * eps;
      const Real dimensional_L = L*sqrt(a*G*M);
      assert( dimensional_eps < 1 );
      assert( 1 < 3.*dimensional_eps/(2.*sqrt(2.)) );
      // Better be safe than sorry, see ulimits appendix
      assert( (sqrt((8. - 36.*pow(dimensional_eps,2) + 27.*pow(dimensional_eps,4) - 
            dimensional_eps*sqrt(pow(-8. + 9.*pow(dimensional_eps,2),3.)))/
          (-1. + pow(dimensional_eps,2)))*Rs)/(2.*sqrt(2.)) > abs(dimensional_L) );
    }

    // Note: eps, L are dimensionless variables
    
    Real ulimits[2]={ERROR};
    u_limits(eps, L, ulimits); // Get limits
    assert( !check(ulimits[0]) );
    if( !check(ulimits[0]) ) return ERROR; // Check for bad results
    
    Real integral = 0;
    
    // Integrate with a singularity:
    {
      // <GSL>
      gsl_integration_workspace * work_space = 
                            gsl_integration_workspace_alloc(1000);
     
      double abs_error = 0;
      double rel_error = 1.e-6;
      size_t limit = 1000;
      double result, error;
     
      gsl_function integrand;
     
      double parameters[2] = {eps, L};
      integrand.function = &integrand_function;
      integrand.params   = &parameters[0];
     
      gsl_integration_qags( &integrand, ulimits[0], ulimits[1],
                            abs_error, rel_error, limit, work_space,
                            &result, &error );
     
      assert( abs(result) > 1e5*abs(error) );
     
      gsl_integration_workspace_free( work_space ); // Free memory
      // </GSL>
      integral = result;
    }
    // Write output
    output << eps << " " << L << " " << ulimits[0] << " " << ulimits[1] << " " << integral << endl;
    return integral;
  }
}

inline
Real f_H( const Real dimless_eps ) {
  const Real result =   sqrt(dimless_eps)/pow(1-dimless_eps,2)*
  ((1.-2.*dimless_eps)*(8*dimless_eps*dimless_eps-8*dimless_eps-3.)+3.*asin(sqrt(dimless_eps))/sqrt(dimless_eps*(1-dimless_eps)))
;
  assert( result >= 0 );
  return result;
}

Real solve_radial_action( Real hernquist_eps, void * params ) {
  // hernquist_eps = dimensionless hernquist energy, as defined in the Sadeghian article
  assert( hernquist_eps > 0 );
  double * parameters = (double*)params;

  // Returns I_bh - I_h
  Real eps = parameters[0];
  Real L = parameters[1];

  // Check input
  {
    const Real dimensional_eps = 1.-(G*M/a) * eps;
    const Real dimensional_L = L*sqrt(a*G*M);
    assert( dimensional_eps < 1 );
    assert( 1 < 3.*dimensional_eps/(2.*sqrt(2.)) );
    // Better be safe than sorry, see ulimits appendix
    assert( (sqrt((8. - 36.*pow(dimensional_eps,2) + 27.*pow(dimensional_eps,4) - 
          dimensional_eps*sqrt(pow(-8. + 9.*pow(dimensional_eps,2),3.)))/
        (-1. + pow(dimensional_eps,2)))*Rs)/(2.*sqrt(2.)) > abs(dimensional_L) );
  }


  Real radialBH = parameters[2];
  // L is not the same here! (maybe)
  Real radialHernquist = Hernquist::II_radial(hernquist_eps, L);
  
  assert( check(radialBH) );
  assert( check(radialHernquist) );
  assert( radialBH > 0 );
  assert( radialHernquist > 0 );

  return radialBH - radialHernquist;
}


Real transform_energy( const Real eps, const Real L ) {
  // Check input
  {
    const Real dimensional_eps = 1.-(G*M/a) * eps;
    const Real dimensional_L = L*sqrt(a*G*M);
    assert( dimensional_eps < 1 );
    assert( 1 < 3.*dimensional_eps/(2.*sqrt(2.)) );
    // Better be safe than sorry, see ulimits appendix
    assert( (sqrt((8. - 36.*pow(dimensional_eps,2) + 27.*pow(dimensional_eps,4) - 
          dimensional_eps*sqrt(pow(-8. + 9.*pow(dimensional_eps,2),3.)))/
        (-1. + pow(dimensional_eps,2)))*Rs)/(2.*sqrt(2.)) > abs(dimensional_L) );
  }

  //This transforms the energy from relativistic basis to newtonian basis by equating the radial action integrals
  Real I_bh = BlackHole::II_radial(eps, L);
  assert(check( I_bh ));
  // Find the new energy
  Real E_dimensionless_newtonian = ERROR;
  
  Real root;
  // Find root:
  //<GSL>
  // Define function;
  gsl_function function;
  Real parameters[3] = {eps,L,I_bh};
  function.function = &solve_radial_action;
  function.params = &parameters[0];

  // Set up solver:
  const gsl_root_fsolver_type * solver_type;
  gsl_root_fsolver * solver;
  solver_type = gsl_root_fsolver_brent; // Bisection method
  solver = gsl_root_fsolver_alloc( solver_type ); // Allocate

  // Set up the solver with limits [0,E_max]
  //save_function_output(limits, E, L); return 0;
  //solve_limits( limits, L, r );
  Real limits[2] = {ERROR};
  // Limits are from 
  //limits[0] = -1*G*M/a*(1-1e-7);
  //limits[1] = -1*G*M/a*1e-7;
  assert( check(limits[0]) );
  gsl_root_fsolver_set(solver, &function, limits[0], limits[1]);

  for( int j = 0; j < 1000; ++j ) {
    // Iterate:
    int status = gsl_root_fsolver_iterate(solver);
    
    cout << __LINE__ << endl;
    
    // Check convergence:  
    if (status == GSL_SUCCESS)
      printf ("Converged:\n");
    
    // Get Solver's current solution:
    root = gsl_root_fsolver_root(solver);
    cout << root << endl;
    if( gsl_root_fsolver_x_upper(solver) - gsl_root_fsolver_x_lower(solver) < 1.0e-30 ) {
      cout << "Breaking" << endl;
      break;
    }
    //limits[0] = gsl_root_fsolver_x_lower(solver);
    //limits[1] = gsl_root_fsolver_x_upper(solver);
  }
  
  gsl_root_fsolver_free( solver );
  //</GSL>
  
  return root;
}

const Real pi = 3.141592653589793;

Real rho(const Real r, const int N) {
  assert( r > 0 );
  Real integral = 0;
  const Real du = 1./(double)(N);
  const Real dz = 1./(double)(N);

  // Calculate integral:
//#pragma omp parallel for schedule(dynamic) reduction(+:integral)
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0+1.0e-5,1-1.0e-5);
  for( int u_i = 0; u_i < N; ++u_i ) {
    for( int z_i = 0; z_i < N; ++z_i ) {
      const Real u = dis(gen);
      const Real z = dis(gen);
      assert( 0<u && u<1 );
      assert( 0<z && z<1 );
 

      // Note: eps, L are now dimensionless as in the Sadeghian article
      const Real epsmax = BlackHole::epsmax(r, 0);
      const Real epsmin = BlackHole::epsmin(r, 0);
      assert( epsmin == 0); // Not 100% sure if this is correct
      const Real eps = epsmin + u*epsmax;
      const Real Lmax = BlackHole::Lmax(r, eps);
      const Real Lmin = BlackHole::Lmin(r, eps);
      const Real L = sqrt( z*pow(Lmax,2)+(1.-z)*pow(Lmin,2) );
      assert( Lmin < Lmax );
      assert( epsmin < epsmax );
      assert( epsmin >= 0 );
      assert( epsmax <= 1 );

      integral += (1.-(G*M / a)*epsmax*u)
                * sqrt((pow(Lmax,2)-pow(Lmin,2))/(1.-z))
                * f_H(transform_energy(u,z));
    }
  }
  //TODO: fix the integral here
  return 1./(sqrt(2.)*pow(2.*pi,2))
       * M/pow(a,3)
       * a*BlackHole::epsmax(r, 0) / (r-2.*G*m)
       * integral;
}

namespace tests {
  namespace schwarzschild {
    template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
    }
    
    void check_u_limits( const Real eps, const Real L ) {
      // This is a function to check when the 4-velocity square root parameter is non-zero;
      // The u roots we use only apply when the following is true
      {
        const Real dimensional_eps = 1.-(G*M/a) * eps;
        const Real dimensional_L = L*sqrt(a*G*M);
        assert( dimensional_eps < 1 );
        assert( 1 < 3.*dimensional_eps/(2.*sqrt(2.)) );
        // Better be safe than sorry, see ulimits appendix
        assert( (sqrt((8. - 36.*pow(dimensional_eps,2) + 27.*pow(dimensional_eps,4) - 
              dimensional_eps*sqrt(pow(-8. + 9.*pow(dimensional_eps,2),3.)))/
            (-1. + pow(dimensional_eps,2)))*Rs)/(2.*sqrt(2.)) > abs(dimensional_L) );
      }
    }
  }

  void record_black_hole_function( const Real r ) {
    ofstream file;
    file.open("./checks/black_hole_radial.check");
    assert( r>0 && r>2*G*m );
    
    // Note: dimensionless
    const Real epsmin = BlackHole::epsmin(r, 0);
    const Real epsmax = BlackHole::epsmax(r, 0);
    const int N = 1000;
    
    for( int i = 0; i < N; ++i ) for( int j = 0; j < N; ++j ) {
      const Real eps = epsmin + (epsmax-epsmin)*i/(Real)(N+1);
      const Real Lmin = BlackHole::Lmax(r,eps);
      const Real Lmax = BlackHole::Lmax(r,eps);
      const Real L = Lmin + (Lmax-Lmin)*j/(Real)(N+1);
      const Real I_bh = BlackHole::II_radial(eps, L);
      schwarzschild::check_u_limits(eps,L);
      assert( eps > 0 && L > 0 );
      assert( Lmin <= L && L <= Lmax );
      assert( epsmin <= eps && eps <= epsmax );
      assert( check(I_bh) );
      // Record values
      file << eps << " " << L << " " << I_bh << endl;
    }
    file.close();
    return;
  }
  
  void test_limits( const Real r ) {
    for( int i = 0; i < 10000; ++i ) for( int j = 0; j < 10000; ++j ) {
      
    }
    return;
  }
  
  void tests() {
    // Function for functions
    //Plot radial action
    const Real r = 16 * G * m;
    test_limits(r);
    record_black_hole_function(r);
  }
}

int main(int argc, char **argv) {
  if( argc>1 && strcmp(argv[0], "test") ) {
    //Tests
    cout << "Test:" << endl;
    cout << "Newtonian: Should be 0.154553 (last output)" << endl;
    cout << Hernquist::II_radial(-2.1e-10,7.77632e-8) << endl;
    tests::tests();
    cout << "Tests done" << endl;
  } else {
  
    // Integration:
    const Real r = 20*G*m/pow(c,2);
    cout << rho(r,500) << endl;
    {
      BlackHole::output.open("blackhole_output.dat");
  
      BlackHole::output.close();
    }
  }
    return 0;
}

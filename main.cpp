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

typedef double Real;
using namespace std;

const Real ERROR = numeric_limits<Real>::quiet_NaN();
const bool check( const Real ERROR ) { 
  if( ERROR != ERROR) return false;
  return true;
}

//Constants
const Real a = 20000;
const Real M = 1e12;
const Real m = 4e6;
const Real G = 4.49e-15;
const Real c = 0.3064;
const Real mp = 1e-3;
const Real Rs=2*G*m/pow(c,2);

// Calculate radial action:
namespace Hernquist {
  ofstream output;

  
  void R_limits( const Real eps, const Real L, Real Rlimits[3] ) {
    // Solve third order polynomial:
    double results[MAXDEGREE] = {0};
    double results_im[MAXDEGREE] = {0};
    double coeffs[MDP1];
    int degree = 3;
    
    // See hernquist.nb
    // We first solve limits for r
    // c[0]u^3+...
    //THESE COEFFICIENTS WERE CHECKED 20.11.
    coeffs[0] = 2*eps*mp;

    coeffs[1] = (2*a*eps*mp+2*G*M*pow(mp,2));

    coeffs[2] = -1*pow(L,2);

    coeffs[3] = -1*a*pow(L,2);
    
    // Solve:
    rpoly_ak1(coeffs, &degree, results, results_im);
    sort( results, results + 3 );
    //qsort( results,3,sizeof(Real),compare);
    
    // Check results:
    for( int  i = 0; i < 3; ++i ) { 
      if( results_im[i] != 0) { Rlimits[0]=ERROR; Rlimits[1]=ERROR; return;} 
    }
    if( results[0] > 0 ) {
      cout << "VALUES:" << endl;
      for( int i = 0; i < 3; ++i ) { cout << results[i] << " "; }
      cout << endl;
    }
    assert( results[0] <= 0);
    if( results[0] > 0 ) { Rlimits[0]=ERROR; Rlimits[1]=ERROR; cout << "bad data" << endl; return;}
    if( results[0] > results[1]){ Rlimits[0]=ERROR; Rlimits[1]=ERROR; cout << "errordata" << endl; return; }
    // Pick the limits as per article:
    Rlimits[0] = results[0];
    Rlimits[1] = results[1];
    Rlimits[2] = results[2];
    assert( Rlimits[2] >= Rlimits[1] );
    return;
  }
   
  Real integrand_function( Real r, void * params ){
    // Get parameters
    const double * parameters;
    parameters = (double*) params;
    //These are inputted via the integration function (GSL)
    const Real eps = parameters[0];
    const Real L = parameters[1];
    
    // See Newtonian.nb
    const Real sqrtParameter = (2*eps)/mp - pow(L,2)/(pow(mp,2)*pow(r,2)) + 
   (2*G*M)/(a + r);
    // If imaginary, no contribution to the integral
    if( sqrtParameter < 0 ) { cerr << "Bad sqrtParameter" << endl; return 0; }
    const Real fourVelocity = sqrt(sqrtParameter);
    const Real integrand = fourVelocity;
    
    assert( integrand >= 0 );
    
    //TODO: Check whether integrand is velocity or momentum!
    return integrand;
  }
  
  Real II_radial( const Real eps, const Real LL ) {
    // Note: eps = Newtonian energy
    // And: LL = Newtonian angular momentum
    assert( 0<eps && eps<1 );
    Real Rlimits[3]={ERROR};
    R_limits(eps, LL, Rlimits); // Get limits
    if( !check(Rlimits[0]) ) return ERROR; // Check for bad results
    
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
     
      double parameters[2] = {eps,LL};
      integrand.function = &integrand_function;
      integrand.params   = &parameters[0];
     
      // Please refer to hernquist.nb and Sadeghian article for these
      const Real Rmin = Rlimits[1];
      const Real Rmax = Rlimits[2];
      
      gsl_integration_qags( &integrand, Rlimits[1], Rlimits[2],
                            abs_error, rel_error, limit, work_space,
                            &result, &error );
     
      assert( abs(result) > 1e5*abs(error) );
     
      gsl_integration_workspace_free( work_space ); // Free memory
      // </GSL>
      integral = result;
    }
    // Write output
    output << eps << " " << LL << " " << Rlimits[0] << " " << Rlimits[1] << " " << integral << endl;
    return integral;
  }
}


// Calculate radial action:
namespace BlackHole {
  ofstream output;

  // Note; eps_min and eps_max have now the energy dimension (joule)
  Real eps_min( const Real r, const Real eps ) {
    assert( r>0 && r>3*G*m );
    if( r>=6*G*m ) {
      return (pow(c,2)*mp*(1 + (2*G*m)/r))/sqrt(1 + (6*G*m)/r);
    } else if( 4*G*m <= r ) {
      return (pow(c,2)*mp*(1 - (2*G*m)/r))/sqrt(1 - (3*G*m)/r);
    } else {
      assert( 1==0 );
    }
    return ERROR;
  }
  
  Real eps_max( const Real r, const Real eps ) {
    assert( r>0 && r>3*G*m);
    return mp*pow(c,2);
  }
  
  Real L_min( const Real r, const Real eps ) {
    return 4*sqrt(2)*c*sqrt((pow(G,2)*pow(m,2))/
     (-8 - (27*pow(eps,4))/(pow(c,8)*pow(mp,4)) + 
       (36*pow(eps,2))/(pow(c,4)*pow(mp,2)) + 
       (eps*pow(-8 + (9*pow(eps,2))/
         (pow(c,4)*pow(mp,2)),1.5))/(pow(c,2)*mp))) * mp;
  }
  
  Real L_max( const Real r, const Real eps ) {
    return c*mp*r*sqrt(-1 + (pow(eps,2)*r)/
      (pow(c,4)*pow(mp,2)*(-2*G*m + r)));
  }
  
  Real dimless_eps_max( const Real r, const Real eps ) {
    return (a*(-1.*eps_min(r,eps) + pow(c,2)*mp))/(G*M*mp);
  }
  
  Real dimless_eps_min( const Real r, const Real eps ) {
    return (a*(-1.*eps_max(r,eps) + pow(c,2)*mp))/(G*M*mp);
  }
  
  Real dimless_L_max( const Real r, const Real eps ) {
    return L_min(r,eps)/(sqrt(a*G*M)*mp);
  }
  
  Real dimless_L_min( const Real r, const Real eps ) {
    return L_min(r,eps)/(sqrt(a*G*M)*mp);
  }

  
  void u_limits( const Real eps, const Real L, Real ulimits[2] ) {
    // Solve third order polynomial:
    double results[MAXDEGREE] = {0};
    double results_im[MAXDEGREE] = {0};
    double coeffs[MDP1];
    int degree = 3;
    
    // See FourVelocity.nb
    // We first solve limits for u
    // c[0]u^3+...
    //Changed 20.11.
    const Real mp= 1; // Mass of the particle
    const Real c = 1; // Speed of light
    coeffs[0] = -1.*pow(L,4)*Rs/(pow(eps,2)*pow(mp,2));

    coeffs[1] = pow(L,4)/(pow(eps*mp,2));

    coeffs[2] = -1.*pow(c,2)*pow(L,2)*Rs/pow(eps,2);

    coeffs[3] = pow(c,2)*pow(L,2)/pow(eps,2)-pow(L,2)/(pow(c,2)*pow(mp,2));
    
    // Solve:
    rpoly_ak1(coeffs, &degree, results, results_im);
    // Check results:
    for( int  i = 0; i < 3; ++i ) { 
      if( results_im[i] != 0) { ulimits[0]=ERROR; ulimits[1]=ERROR; return;} 
    }
    if( results[0] <= 0 ) { ulimits[0]=ERROR; ulimits[1]=ERROR; cout << "bad data" << endl; return;}
    if( results[0] > results[1]){ ulimits[0]=ERROR; ulimits[1]=ERROR; cout << "errordata" << endl; return; }
    ulimits[0] = results[0];
    ulimits[1] = results[1];
    assert( ulimits[1] >= ulimits[0] );
    return;
  }
  
  Real integrand_function( Real u, void * params ){
    // Get parameters
    const double * parameters;
    parameters = (double*) params;
    const Real eps = parameters[0];
    const Real L = parameters[1];
    
    // See Schwarzschild.nb
    const Real sqrtParameter =pow(c,2)*(-1 + Rs*u) + 
   (pow(eps,2) + pow(c,2)*pow(L,2)*pow(u,2)*
       (-1 + Rs*u))/(pow(c,2)*pow(mp,2));
    if( sqrtParameter < 0 ) { cerr << "Bad sqrtParameter" << endl; return 0; }
    const Real fourVelocity = -1.0*sqrt(sqrtParameter);
    const Real integrand = -1.0*fourVelocity / pow(u,2);
    
    assert( integrand >= 0 );
    
    //TODO: Check the division!
    return integrand/(1.0-Rs*u);
  }
  
  Real II_radial( const Real eps, const Real LL ) {
    // Note: eps = relativistic energy per particle 0<eps<1
    // And: LL = Relativistic angular momentum per unit mass (larger than 2 Rs)
    assert( 0<eps && eps<1 );
    Real ulimits[2]={ERROR};
    u_limits(eps, LL, ulimits); // Get limits
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
     
      double parameters[2] = {eps, LL};
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
    output << eps << " " << LL << " " << ulimits[0] << " " << ulimits[1] << " " << integral << endl;
    return integral;
  }
}

inline
Real f_H( const double dimless_eps ) {
  const double result =   sqrt(dimless_eps)/pow(1-dimless_eps,2)*
  ((1.-2.*dimless_eps)*(8*dimless_eps*dimless_eps-8*dimless_eps-3.)+3.*asin(sqrt(dimless_eps))/sqrt(dimless_eps*(1-dimless_eps)))
;
  assert( result >= 0 );
  return result;
}

Real solve_radial_action( Real x, void * params ) {
  double * parameters = (double*)params;

  // Returns I_bh - I_h
  Real eps = parameters[0];
  Real L = parameters[1];
  Real radialBH = parameters[2];

  return radialBH - Hernquist::II_radial(x, L);
}


Real transform_energy( const Real eps, const Real L ) {
  //This transforms the energy from relativistic basis to newtonian basis by equating the radial action integrals
  Real I_bh = BlackHole::II_radial(eps, L);
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
  // Limits are from -potential(0) to 0
  limits[0] = -1*G*M/a;
  limits[1] = 0;
  gsl_root_fsolver_set(solver, &function, limits[0], limits[1]);

  for( int j = 0; j < 1000; ++j ) {
  // Iterate:
  int status = gsl_root_fsolver_iterate(solver);
  
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
      assert( z != 0 && u != 0 && z != 1 && u != 1 );
      
      const Real eepsmax = BlackHole::dimless_eps_max(r, 0);
      const Real eepsmin = BlackHole::dimless_eps_min(r, 0);
      // Note: eps is now in the dimensions of Joules (energy) and it is not dimensionless
      const Real eps = BlackHole::eps_min(r,0)+u*(BlackHole::eps_max(r,0)-BlackHole::eps_min(r,0));
      const Real LLmax = BlackHole::dimless_L_max(r, eps);
      const Real LLmin = BlackHole::dimless_L_min(r, eps);
      // Note: L is now in the right dimensions too!
      const Real L = sqrt(z*pow(LLmax,2)+(1-z)*pow(LLmin,2))*sqrt(a*G*M)*mp;

      integral += (1.0-(G*M/pow(c,2) / a)*eepsmax*u)
                * sqrt((pow(LLmax,2)-pow(LLmin,2))/(1-z))
                * f_H(transform_energy(u,z));
    }
  }
  return 1./(sqrt(2)*pow(2*pi,2))*(M/pow(a,3)) * (a*BlackHole::dimless_eps_max()/(r-2*g*m/pow(c,2)))
  * integral;
  return 0;
}

int main(int argc, char **argv) {
  //Tests
  cout << "Black hole" << endl;
  cout << BlackHole::II_radial(0.999,2.*Rs*sqrt(15)) << endl;
  cout << Hernquist::II_radial(a*pow(10,-26)/(G*M),2.*Rs*sqrt(15)) << endl;
  // Integration:
  const Real r = 20*G*m/pow(c,2);
  cout << rho(r) << endl;
  {
    BlackHole::output.open("blackhole_output.dat");

    BlackHole::output.close();
  }
    return 0;
}

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
const Real m = pow(10,6); // Solar mass
const Real M = pow(10,12); // Solar mass
const Real a = 20*2.09*pow(10,16); // 20 kpc to solar mass
const Real G = 1;
const Real Rs = 2*m*G; // Schwarzschild radius
const Real mm = m/M;



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
    const Real mp= 1; // Mass of the particle
    const Real c = 1; // Speed of light
    //THESE COEFFICIENTS WERE CHECKED 19.11.
    coeffs[0] = 2.0*G*M*eps/a;

    coeffs[1] = 2.0*G*M*(eps-1);

    coeffs[2] = pow(L,2);

    coeffs[3] = a*pow(L,2);
    
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
    const Real R1 = parameters[0];
    const Real R2 = parameters[1];
    const Real R3 = parameters[2];
    const Real eps = parameters[3];
    
    // See FourVelocity.nb
    const Real c=1; const Real mp=1; // just set it all to 1 :)
    const Real sqrtParameter = -((eps*G*M*(r - R1)*(r - R2)*(r - R3))/
     (a*mp*Power(r,2)*(a + r)));
    if( sqrtParameter < 0 ) { cerr << "Bad sqrtParameter" << endl; return 0; }
    const Real fourVelocity = sqrt(2.0) * sqrt(sqrtParameter);
    const Real integrand = fourVelocity;
    
    assert( integrand >= 0 );
    
    //TODO: Check whether integrand is velocity or momentum!
    return integrand;
  }
  
  Real II_radial( const Real eps, const Real LL ) {
    // Note: eps = relativistic energy per particle 0<eps<1
    // And: LL = Relativistic angular momentum per unit mass (larger than 2 Rs)
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
     
      double parameters[4] = {Rlimits[0], Rlimits[1], Rlimits[2], eps};
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

  void eps_limits( const Real r, Real epsLimits[2] ) {
    if( r>= 6.0 * G * m ) {
      epsLimits[0] = (1.+2*G*m/r)/sqrt(1.+6*G*m/r);
    } else if( r >= 4.0*G*m ) {
      epsLimits[0] = (1.-2*G*m/r)/sqrt(1.-3*G*m/r);
    } else { assert(1==0); }
    epsLimits[1] = 1;
    return;
  }

  void L_limits( const Real eps, const Real r, LLimits[2] ) {
    LLimits[0] = sqrt(32.*pow(G*m,2) / 
        36.*pow(eps,2)-27.*pow(eps,4)-8.+eps*pow(9.*pow(eps,2)-8,3./2.));
    LLimits[1] = r*sqrt(pow(eps,2)/(1.-2.*G*m/r)-1.);
    return;
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
    const Real mp= 1; // Mass of the particle
    const Real c = 1; // Speed of light
    coeffs[0] = -1.*(pow(L,4) * Rs)
               / (pow(eps,2)*pow(mp,2));

    coeffs[1] = pow(L,4)
               / (pow(eps,2)*pow(mp,2));

    coeffs[2] = -1.*pow(c,2) * pow(L,2) * Rs
               / pow(eps,2);

    coeffs[3] = ( pow(c,2) * pow(L,2) )
               / pow(eps,2) 
	        -
	        pow(L,2)
               / ( pow(c,2) * pow(mp,2) );
    
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
    
    // See FourVelocity.nb
    const Real c=1; const Real mp=1; // just set it all to 1 :)
    const Real sqrtParameter =
    pow(c,2)*(Rs*u-1.0)
    +
    (pow(c,2)*pow(L,2)*pow(u,2)*(Rs*u-1.0)+pow(eps,2))
     / (pow(c,2)*pow(mp,2))
    ;
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
Real f_h( const double eps ) {
  const double result =   sqrt(E)/pow(1-E,2)*
  ((1.-2.*E)*(8*E*E-8*E-3.)+3.*asin(sqrt(E))/sqrt(E*(1-E)))
;

  assert( result >= 0 );
  return result;
}

Real rho(const Real r, const int N) {
  assert( r > 0 );
  Real integral = 0;
  const Real du = 1./(double)(N);
  const Real dz = 1./(double)(N);

  // Calculate integral:
//#pragma omp parallel for schedule(dynamic) reduction(+:integral)
  for( int u_i = 0; u_i < N; ++u_i ) {
    for( int z_i = 0; z_i < N; ++z_i ) {
      const Real u = rRand(0,1);
      const Real z = rRand(0,1);

//      //FIXME: 5e-4
//      const Real u = rRand(0,1-8.e-4);
//      const Real z = rRand(0,1-8.e-4);
//      assert( z != 0 && u != 0 && z != 1 && u != 1 );
//      integral += (1.-(G*M/a0)*Rel::E_max(r)*u) *
//                  sqrt((Rel::L_max2(u,r)-Rel::L_min2(u,r))/(1-z)) *
//                  hernquist_f(E_initial(u,z,r)) * du*dz;
//      assert(integral >= 0);
//      //if ((u_i+z_i)%100 == 0)
//      //cout << "Iteration:" << u_i << " " << z_i << endl;
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  //Tests
  cout << "Black hole" << endl;
  cout << BlackHole::II_radial(0.999,2.*Rs*sqrt(15)) << endl;
  cout << Hernquist::II_radial(a*pow(10,-26)/(G*M),2.*Rs*sqrt(15)) << endl;
  // Integration:
  r = 12.0 * G * m;
  rho(r);
  {
    BlackHole::output.open("blackhole_output.dat");

    BlackHole::output.close();
  }
    return 0;
}

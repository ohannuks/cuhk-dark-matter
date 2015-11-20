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
namespace BlackHole {
  ofstream output;
  
  void radial_limits( const Real eps, const Real L, Real Rlimits[2] ) {
    // Solve third order polynomial:
    double results[MAXDEGREE] = {0};
    double results_im[MAXDEGREE] = {0};
    double coeffs[MDP1];
    int degree = 3;
    
    // See limits_r.nb
    // We first solve limits for u=1/Rmax
    // c[0]u^3+...
    coeffs[0] = (pow(L,4) * Rs)
               / pow(eps,2);

    coeffs[1] = -1 * pow(L,4)
               / pow(eps,2);

    coeffs[2] = pow(L,2) * Rs
               / pow(eps,2);

    coeffs[3] = ( pow(L,2)- pow(L,2)/pow(eps,2) );
    
    // Solve:
    rpoly_ak1(coeffs, &degree, results, results_im);
    // Check results:
    for( int  i = 0; i < 3; ++i ) { if( results_im[i] != 0) { 
      Rlimits[0]=ERROR; Rlimits[1]=ERROR; return;} }
    if( results[0] <= 0 ) { Rlimits[0]=ERROR; Rlimits[1]=ERROR; cout << "bad data" << endl; return;}
    if( results[0] > results[1]){ Rlimits[0]=ERROR; Rlimits[1]=ERROR; cout << "errordata" << endl; return; }
    Rlimits[0] = 1.0/results[1];
    Rlimits[1] = 1.0/results[0];
    assert( Rlimits[1] >= Rlimits[0] );
    return;
  }
  
  Real integrand_function( Real r, void * params ){
    // Get parameters
    const double * parameters;
    parameters = (double*) params;
    const Real eps = parameters[0];
    const Real LL = parameters[1];
    
    // See main.nb
    const Real result = (pow(eps,2)-(1.0-Rs/r)*(1.0+pow(LL,2)/pow(r,2)));
    
    // Check for bad result:
    if( result < 0 ) {
      cout << "Bad result in integrand of black hole" << endl;
      return 0; // Dont break the code, just dont affect the integral with bad results
    }
    
    return sqrt( result ) / (1.0-Rs/r);
  }
  
  Real II_radial( const Real eps, const Real LL ) {
    // Note: eps = relativistic energy per particle 0<eps<1
    // And: LL = Relativistic angular momentum per unit mass (larger than 2 Rs)
    assert( 0<eps && eps<1 );
    Real Rlimits[2];
    radial_limits(eps, LL, Rlimits); // Get limits
    if( !check(Rlimits[0]) ) return ERROR; // Check for bad results
    
    Real integral = 0;
    
    // Integrate with a singularity:
    {
      // <GSL>
      gsl_integration_workspace * work_space = 
                            gsl_integration_workspace_alloc(1000);
     
      double abs_error = 0;
      double rel_error = 1.e-5;
      size_t limit = 1000;
      double result, error;
     
      gsl_function integrand;
     
      double parameters[2] = {eps, LL};
      integrand.function = &integrand_function;
      integrand.params   = &parameters[0];
     
      gsl_integration_qags( &integrand, Rlimits[0], Rlimits[1],
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
namespace BlackHole2 {
  ofstream output;
  
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


int main(int argc, char **argv) {
  //Tests
  {
    BlackHole::output.open("blackhole_output.dat");
    Real rlimits[2]={0};
    BlackHole2::u_limits(0.9999,10*Rs,rlimits);
    //cout << 1./rlimits[0] << " " << 1./rlimits[1] << endl;
    BlackHole::radial_limits(0.9999,10*Rs,rlimits );
    cout << rlimits[0] << " " << rlimits[1] << endl;
    const Real epsmin=sqrt(8./9.);
    const Real LLmin=2*Rs*sqrt(12);
    const Real LLmax=2*Rs*sqrt(29);
    Real sum=0;
    Real sum2=0;
    for( Real eps = epsmin; eps < 1; eps+= (1.-epsmin)/1000. ) {
      for( Real LL = LLmin; LL < LLmax; LL+= (LLmax-LLmin)/1000. ) {
	const Real term1=BlackHole::II_radial( eps, LL );
	const Real term2=BlackHole2::II_radial( eps, LL );
//	if( check(term1) ) { sum += term1; }
//	if( check(term2) ) { sum2+= term2; }
//        if( check(term1) || check(term2) ) {
//          cout << "R action: ";
//          if( check(term1) ) { cout << term1 << " "; }
//          cout << "u action: ";
//          if( check(term2) ) { cout << term2 << " "; }
//          cout << endl;
//        }
      }
    }
    cout << sum << endl;
    cout << sum2 << endl;
    BlackHole::output.close();
  }
    return 0;
}

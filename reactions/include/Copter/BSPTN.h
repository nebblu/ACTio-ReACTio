#ifndef BSPTN_H
#define BSPTN_H

#include <vector>
#include <gsl/gsl_sf_bessel.h>
#include "SpecialFunctions.h"
#include "BeyondLCDM.h"


class BSPTN {
public:
  // pars: base parameters. Currently used: 0: scale factor, 1: total matter fraction today
  // extpars: extended parameters. extpars[0] = Omega_rc for nDGP and fr0 for f(R)
  // epars: ODE error parameters - see specific function or examples/bs.cpp example

  // tree level  equation solvers : DGP, f(R)
  void initnb0_dgp(double pars[], double extpars[], double k[], double x[], double kargs[]);
  void initnb0_fofr(double pars[], double extpars[], double k[], double x[], double kargs[]);

  // 1-loop level  equation solvers : DGP, f(R)
  void initnb1_dgp(double pars[], double extpars[], double epars[], double k[], double x[], double kargs[]);
  void initnb1_fr(double pars[], double extpars[], double epars[], double k[], double x[], double kargs[]);

};

// NUMERICAL KERNELS /////

/// CREATION OF KERNEL ARRAYS ////

//Kernel array sizes

// 1st ORDER
extern double F1b_k1;
extern double G1b_k1;
extern double F1b_k2;
extern double G1b_k2;
extern double F1b_k3;
extern double G1b_k3;

extern double F1b_1;
extern double G1b_1;

extern double F1b_2;
extern double G1b_2;

extern double F1b_3;
extern double G1b_3;


// 2nd order
extern double F2b_k12;
extern double G2b_k12;
extern double F2b_k13;
extern double G2b_k13;
extern double F2b_k23;
extern double G2b_k23;

extern double F2b_p2mp;
extern double G2b_p2mp;

extern double F2b_p3mp;
extern double G2b_p3mp;

extern double F2b_12a;
extern double G2b_12a;

extern double F2b_13a;
extern double G2b_13a;

extern double F2b_23a;
extern double G2b_23a;

// 3rd order

extern double F3b_12mp;
extern double G3b_12mp;

extern double F3b_13mp;
extern double G3b_13mp;

extern double F3b_21mp;
extern double G3b_21mp;

extern double F3b_23mp;
extern double G3b_23mp;

extern double F3b_31mp;
extern double G3b_31mp;

extern double F3b_32mp;
extern double G3b_32mp;

extern double F3b_1pp;
extern double G3b_1pp;

extern double F3b_2pp;
extern double G3b_2pp;

extern double F3b_3pp;
extern double G3b_3pp;

// 4th order
extern double F4b_12pp;
extern double G4b_12pp;

extern double F4b_13pp;
extern double G4b_13pp;
extern double F4b_23pp;
extern double G4b_23pp;

const double myh0sqr = pow2(1./2997.92);


#endif // BSPTN_H

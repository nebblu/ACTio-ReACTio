#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cfloat>
#include <cmath>

#include "Common.h"
#include "Quadrature.h"
#include "BeyondLCDM.h"
#include "SpecialFunctions.h"

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_min.h>


#include <math.h>       /* pow */
#include <iostream>
#include <stdlib.h>
#include <functional>

using std::cref;
using std::bind;


/* Some useful functions */

//1-mu^2
static double ker1(double u1){
	return 1.-u1*u1;
}

//beta function for nDGP
// omega_rc as described in 1606.02520  for example
inline double beta(double a, double omega0, double omegarc){
	return 1. +  HA(a,omega0)/sqrt(omegarc)*(1.+HA1(a,omega0)/(3.*HA(a,omega0)*HA(a,omega0)));}


/* All modified background and Poisson equation functions */

/* GENERAL BACKGROUND EXPANSION
Used in SpecialFunctions.cpp 's
initn_lin
initn2
initn2_pseudo
initn3
initn_rsd  */

/* See 1606.02520, 1808.01120, 2005.12184  for more details */

// model selects MG or DE model (1 = GR, MG: 2 = f(R), 3 = DGP, DE models: 4 = quintessence, 5 = CPL, 6 = CPL with dark scattering)
// 7-10: EFTofDE with PPF, unscreened, superscreened and KGB non-linear implementations.

// Throughout the values of pars, extpars and model are:
// pars: base parameters. Currently used: 0: scale factor, 1: total matter fraction today, 2: total massive neutrino fraction today
// extpars: extended parameters.
// extpars[0] = Omega_rc for nDGP, fr0 for f(R), w0 for CPL
// extpars[1] = wa for CPL
// extpars[2] = xi * h for Dark scattering
// for EFTofDE and LCDM background we have:
// extpars[0-3] = alpha_{K0},alpha_{B0},alpha_{M0},M^2 resp
// for KGB with CPL background we have:
// extpars[0-4] = w0, wa, alpha_{K0},alpha_{B0},alpha_{M0}


                                        ////
                                    /////////////
                            //////////////////////////////
                      /////////////////////////////////////////
                ///////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////



/* BACKGROUND QUANTITIES * /


/*  Scale factor evolution of alpha_i for EFTofDE */
// edit as necessary - but also edit the analytic scale factor derivative function accordingly
inline double alphai_eft(double a, double omega0, double alpha0){
	return alpha0*a;
}

// scale factor derivatives of alphai_eft
inline double dalphai_eft(double a, double omega0, double alpha0){
	return alpha0;
}


double HAg(double a, double omega0, double extpars[], int model){
	double A, omegaf, omegaL;
	switch(model) {
		case 1:
		/* LCDM */
		 return HA( a, omega0);

		 case 2:
		 /* f(R) Hu-Sawicki */
 		 return HA( a, omega0);

		 case 3:
		 /* DGP */
			return HA( a, omega0);

			case 4:
			/* QUINTESSENCE */
			A = -3.*(1.+extpars[0]);
			omegaf = pow(a,A);
			omegaL= (1.-omega0)*omegaf;
			return  sqrt(omega0/pow(a,3)+omegaL);

			case 5:
			/* CPL */
			A = -3.*(1.+extpars[0]+extpars[1]);
			omegaf = pow(a,A)*exp(3.*(-1.+a)*extpars[1]);
			omegaL= (1.-omega0)*omegaf;
			return  sqrt(omega0/pow(a,3)+omegaL);

			case 6:
			/* Dark Scattering with CPL  */
			A = -3.*(1.+extpars[0]+extpars[1]);
			omegaf = pow(a,A)*exp(3.*(-1.+a)*extpars[1]);
			omegaL= (1.-omega0)*omegaf;
			return  sqrt(omega0/pow(a,3)+omegaL);

			case 7:
			/* EFTofDE: PPF - set to LCDM background for now  */
			return HA( a, omega0);

			case 8:
			/* EFTofDE: unscreened approx */
			return  HA( a, omega0);

			case 9:
			/* EFTofDE: superscreened approx */
			return  HA( a, omega0);

			case 10:
			/* KGB with CPL background  */
			A = -3.*(1.+extpars[0]+extpars[1]);
			omegaf = pow(a,A)*exp(3.*(-1.+a)*extpars[1]);
			omegaL= (1.-omega0)*omegaf;
			return  sqrt(omega0/pow(a,3)+omegaL);

			default:
					warning("SpecialFunctions: invalid model choice, model = %d \n", model);
					return 0;
}
}

/* Must specify analytic derivatives - want to avoid numerical derivatives*/

// aH dH/da / H0^2
double HA1g(double a, double omega0, double extpars[], int model){
	double A, omegaf, omegaL;
	switch(model) {
		case 1:
		/* LCDM */
		 return HA1( a, omega0);

		 case 2:
		 /* f(R) Hu-Sawicki */
 		 return HA1( a, omega0);

		 case 3:
		 /* DGP */
			return HA1( a, omega0);

			case 4:
		/*QUINTESSENCE*/
		A = -3.*(1.+extpars[0]);
		omegaf = pow(a,A);
		omegaL= (1.-omega0)*omegaf;
		return -3*omega0/(2.*pow(a,3)) + A*omegaL/2.;

		 case 5:
	  /* CPL */
		A = -3.*(1.+extpars[0]+extpars[1]);
	  omegaf = pow(a,A)*exp(3*(-1.+a)*extpars[1]);
		omegaL= (1.-omega0)*omegaf;
	  return -3.*omega0/(2.*pow(a,3)) + (A+3.*a*extpars[1])*omegaL/2.;

		case 6:
		/* Dark Scattering with CPL  */
		A = -3.*(1.+extpars[0]+extpars[1]);
	  omegaf = pow(a,A)*exp(3*(-1.+a)*extpars[1]);
		omegaL= (1.-omega0)*omegaf;
	  return -3.*omega0/(2.*pow(a,3)) + (A+3.*a*extpars[1])*omegaL/2.;

		case 7:
		/* EFTofDE: PPF - set to LCDM background for now  */
		return HA1( a, omega0);

		case 8:
		/* EFTofDE: unscreened approximation */
		return HA1( a, omega0);

		case 9:
		/* EFTofDE: superscreened approximation */
		return HA1( a, omega0);

		case 10:
		/* KGB with CPL background  */
		A = -3.*(1.+extpars[0]+extpars[1]);
	  omegaf = pow(a,A)*exp(3*(-1.+a)*extpars[1]);
		omegaL= (1.-omega0)*omegaf;
	  return -3.*omega0/(2.*pow(a,3)) + (A+3.*a*extpars[1])*omegaL/2.;

		default:
				warning("SpecialFunctions: invalid model choice, model = %d \n", model);
				return 0;
}
}

// some useful definitions for SpecialFunctions.cpp

//HA2g = -dH/dt/H^2 = -a dH/da / H
double HA2g(double a, double omega0, double extpars[], int model){
   	return -HA1g(a,omega0,extpars,model)/pow2(HAg(a,omega0,extpars,model));
}

//3 H0^2/(2H^2) * Omega_{m,0}/a^3
double HA2g2(double a, double omega0, double extpars[], int model){
   	return 3.*omega0/(2.*pow2(HAg(a,omega0,extpars,model))*pow(a,3));
}


// Friction term appearing in the Euler equation from dark matter - dark energy interactions (see: Eq. 8 of 2111.13598)
double myfricF(double a, double omega0, double extpars[], int model){
	double deexponent,omegaf, omegaL;
	switch(model) {
		case 1:
		/* LCDM */
						return  0. ;
		case 2:
		/* f(R) Hu-Sawicki */
						return  0. ;
		case 3:
		/* DGP */
						return  0. ;
		case 4:
		/* Quintessence */
						return 0. ;
		case 5:
		/* CPL */
						return 0.;
		case 6:
		/* Dark Scattering with CPL  */
		 				deexponent = -3.*(1.+extpars[0]+extpars[1]);
		 				omegaf = pow(a,deexponent)*exp(3.*(-1.+a)*extpars[1]);
		 				omegaL= (1.-omega0)*omegaf;
						return  (1.+ (extpars[0]+(1.-a)*extpars[1])) * omegaL * extpars[2] * 0.0974655;
						// = (1+w(a)) * Omega_L * p3 * 3 H / (8 \pi G) * unit conversion (p3 = xi * h ) : Omega_L = Omega_L,0 * H_0^2 / H^2 * evolution
		case 7:
		/* EFTofDE: PPF */
					return 0.;
		case 8:
		/* EFTofDE: unscreened approx */
					return 0.;
		case 9:
		/* EFTofDE: superscreened approx */
					return 0.;
		case 10:
		/* EFTofDE: KGB with CPL background */
					return 0.;
		default:
					warning("SpecialFunctions: invalid model choice, model = %d \n", model);
				  return 0;
  }
}



/*MODIFIED GRAVITY MODEL FUNCTIONS (see 1606.02520) */


/* PERTURBATIVE POISSON EQUATION FUNCTIONS UP TO 3RD ORDER */
//p1,p2,p3 are theory parameters
//k1,k2,k3 are the mode magnitudes
//u1 cosine of the angle between k2 and k3; u2 between k1 and k2; u3 between k1 and k3

// Linear modification mu
double mu(double a, double k0, double omega0, double extpars[], int model){
	double h0 = 1./2997.92458;
	double var1, var2, alphaofa[5],dalphaofa[5];
	switch(model) {
		case 1:
		/* GR */
						return  1. ;
		case 2:
		/* f(R) Hu-Sawicki */
						var1 = pow2(k0/a);
	  				return 1. + var1/(3.*(var1+pow3(omega0/pow3(a)-4.*(omega0-1.))/(2.*extpars[0]/pow2(h0)*pow2(4.-3.*omega0))));
    case 3:
		/* nDGP */
						return 1.+1./(3.*beta(a,omega0,extpars[0]));
		case 4:
		/* QUINTESSENCE */
						return  1. ;
		case 5:
		/* CPL */
						return  1. ;
		case 6:
		/* Dark Scattering with CPL  */
						return  1. ;

		/* k->infinity limit of EFTofDE. See Eq. 26 of 1606.05339*/
		// 7-10: PPF, unscreened, superscreened, KGB
		case 7:
						alphaofa[0] = alphai_eft(a,omega0,extpars[0]); // alpha_K(a) = alpha_{K0}*a
						alphaofa[1] = alphai_eft(a,omega0,extpars[1]); // alpha_B(a) = alpha_{B0}*a
						alphaofa[2] = alphai_eft(a,omega0,extpars[2]); // alpha_M(a) = alpha_{M0}*a
						alphaofa[3] = alphai_eft(a,omega0,extpars[3]); // M^2(a)

						// scale factor derivatives of alpha_i(a)
						dalphaofa[0] = dalphai_eft(a,omega0,extpars[0]);
						dalphaofa[1] = dalphai_eft(a,omega0,extpars[1]);
						dalphaofa[2] = dalphai_eft(a,omega0,extpars[2]);
						dalphaofa[3] = dalphai_eft(a,omega0,extpars[3]);

					  var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha

						var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[2] + alphaofa[1]/2. + HA2g(a,omega0,extpars,model))
																					 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)); //cs^2


					  return (1. + pow2(alphaofa[1] + 2.*alphaofa[2])/(2.*var2*var1))/alphaofa[3];

	  case 8:
						alphaofa[0] = alphai_eft(a,omega0,extpars[0]); // alpha_K(a) = alpha_{K0}*a
						alphaofa[1] = alphai_eft(a,omega0,extpars[1]); // alpha_B(a) = alpha_{B0}*a
						alphaofa[2] = alphai_eft(a,omega0,extpars[2]); // alpha_M(a) = alpha_{M0}*a
						alphaofa[3] = alphai_eft(a,omega0,extpars[3]); // M^2(a)

						// scale factor derivatives of alpha_i(a)
						dalphaofa[0] = dalphai_eft(a,omega0,extpars[0]);
						dalphaofa[1] = dalphai_eft(a,omega0,extpars[1]);
						dalphaofa[2] = dalphai_eft(a,omega0,extpars[2]);
						dalphaofa[3] = dalphai_eft(a,omega0,extpars[3]);

						var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha

						var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[2] + alphaofa[1]/2. + HA2g(a,omega0,extpars,model))
																					 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)); //cs^2


						return (1. + pow2(alphaofa[1] + 2.*alphaofa[2])/(2.*var2*var1))/alphaofa[3];

		case 9:
						alphaofa[0] = alphai_eft(a,omega0,extpars[0]); // alpha_K(a) = alpha_{K0}*a
						alphaofa[1] = alphai_eft(a,omega0,extpars[1]); // alpha_B(a) = alpha_{B0}*a
						alphaofa[2] = alphai_eft(a,omega0,extpars[2]); // alpha_M(a) = alpha_{M0}*a
						alphaofa[3] = alphai_eft(a,omega0,extpars[3]); // M^2(a)

						// scale factor derivatives of alpha_i(a)
						dalphaofa[0] = dalphai_eft(a,omega0,extpars[0]);
						dalphaofa[1] = dalphai_eft(a,omega0,extpars[1]);
						dalphaofa[2] = dalphai_eft(a,omega0,extpars[2]);
						dalphaofa[3] = dalphai_eft(a,omega0,extpars[3]);

						var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha

						var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[2] + alphaofa[1]/2. + HA2g(a,omega0,extpars,model))
																					 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)); //cs^2


						return (1. + pow2(alphaofa[1] + 2.*alphaofa[2])/(2.*var2*var1))/alphaofa[3];
		case 10:
		// CPL background so we use extpars[0]=w0, extpars[1]=wa
						alphaofa[0] = alphai_eft(a,omega0,extpars[2]); // alpha_K(a) = alpha_{K0}*a
						alphaofa[1] = alphai_eft(a,omega0,extpars[3]); // alpha_B(a) = alpha_{B0}*a
						alphaofa[2] = alphai_eft(a,omega0,extpars[4]); // alpha_M(a) = alpha_{M0}*a
						// scale factor derivatives of alpha_i(a)
						dalphaofa[0] = dalphai_eft(a,omega0,extpars[2]);
						dalphaofa[1] = dalphai_eft(a,omega0,extpars[3]);
						dalphaofa[2] = dalphai_eft(a,omega0,extpars[4]);

					  var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha
						var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[1]/2. + HA2g(a,omega0,extpars,model))
																					 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)); //cs^2
					  return 1. + pow2(alphaofa[1])/(2.*var2*var1);

		default:
				warning("SpecialFunctions: invalid model choice, model = %d \n", model);
				return 0;
}
}


// Gamma 2
double gamma2(double a, double omega0, double k0, double k1, double k2, double u1, double extpars[], int model){
	double h0 = 1./2997.92458;
	double var1,var2,var3,var4,var5,var6;
	switch(model) {
		case 1:
		/* GR */
						return  0. ;
		case 2:
		/* f(R) Hu- Sawicki */
						var1 = pow3(a);
						var2 = pow3(var1);
						var3 = pow2(h0);
						var4 = pow2(3.*omega0-4.);
						var5 = pow3(omega0-4.*var1*(omega0-1.))/(2.*var2*extpars[0]/var3*var4);
						var6 = pow2(k0/a);
							 return -(9.*var6*pow2(omega0/var1)*pow(omega0-4.*var1*(-1.+omega0),5))/
									    (48.*pow(var1,5)*pow2(extpars[0]/var3)*pow2(HA(a,omega0))*pow2(var4)
									   *(var6+var5)
									   *(pow2(k1/a)+var5)
							 		   *(pow2(k2/a)+var5));
		case 3:
		/* DGP */
						return -1./(HA(a,omega0)*HA(a,omega0)*24.*pow(beta(a,omega0,extpars[0]),3)*extpars[0])*pow(omega0/(a*a*a),2)*ker1(u1);
		case 4:
		/* QUINTESSENCE */
						return  0. ;
		case 5:
		/* CPL */
						return  0. ;
		case 6:
		/* Dark Scattering with CPL  */
						return  0. ;
		case 7:
		/* EFTofDE PPF*/
						return  0. ;
		case 8:
		/* EFTofDE - unscreened  approximation  */
						return  0. ;
		case 9:
		/* EFTofDE - superscreened approximation  */
						return  0. ;
		case 10:
		/* KGB  with CPL background*/
						return  0. ;
		default:
		warning("SpecialFunctions: invalid model choice, model = %d \n", model);
				return 0;
}
}

// Partially symmetric gamma3 (in last 2 arguments)
double gamma3(double a, double omega0, double k0, double k1, double k2, double k3, double u1,double u2, double u3, double extpars[], int model){
	double h0 = 1./2997.92458;
	double var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,k23 ;

	switch(model) {
	  case 1:
		/* GR */
						return  0. ;
	  case 2:
		/* f(R) Hu- Sawicki */
						var1 = pow3(a);
						var2 = pow3(var1);
						var3 = pow2(h0);
						var4 = pow2(3.*omega0-4.);
						var5 = pow3(omega0-4.*var1*(omega0-1.))/(2.*var2*extpars[0]/var3*var4);
						var6 = pow(var1,5);
						var7 = var6*pow2(var1);
						k23 = sqrt(k2*k2+k3*k3+2.*k2*k3*u1);
						var8 = pow2(k0/a);
						var9 = pow2(k23/a);
						var10 = pow2(extpars[0]/var3);
						var11 = pow2(var4);

						  return 1./3.*(var8*pow2(1./HA(a,omega0))*pow3(omega0/var1)/(36.
						 										*(var8+var5)
						 										*(pow2(k1/a)+var5)
						 										*(pow2(k2/a)+var5)
						 										*(pow2(k3/a)+var5)
						 				 					  *(var9+var5))
						 		 								*(-45./8.*(var9+var5)/(var7*var10*(extpars[0]/var3))*pow(omega0-4.*var1*(omega0-1.),7)/var11/var4
						 		   							+pow2(9./(4.*var6*var10)*pow(omega0-4.*var1*(omega0-1.),5)/var11)));
	  case 3:
		/* nDGP */
							return 1./(HA(a,omega0)*HA(a,omega0)*144.*pow(beta(a,omega0,extpars[0]),5)*pow(extpars[0],2))*pow(omega0/(a*a*a),3)*ker1(u1)*ker1((k2*u2+k3*u3)/sqrt(k2*k2+2.*k2*k3*u1+k3*k3));
		case 4:
		/* QUINTESSENCE */
						return  0. ;
		case 5:
		/* CPL */
						return  0. ;
		case 6:
		/* Dark Scattering with CPL  */
						return  0. ;
		case 7:
		/* EFTofDE PPF  */
						return  0. ;
		case 8:
		/* EFTofDE - unscreened  approximation  */
						return  0. ;
		case 9:
		/* EFTofDE - superscreened approximation  */
						return  0. ;
		case 10:
		/* KGB  with CPL background*/
						return  0. ;
		default:
		warning("SpecialFunctions: invalid model choice, model = %d \n", model);
				return 0;
}
}


/* Modified non-linear Poisson equation function (F) for spherical collapse (see for example A.1 of appendix of 1812.05594) */

double mymgF(double a, double yh, double yenv, double Rth, double omega0, double extpars[], double delta_initial, int model){
	double h0 = 1./2997.92;
	double dod, dod2, dRRth, fr0, var1, var2, var3, term1, term2;
	double betadgp,xm3,xterm,delta,deltaenv,Mvir;
	double alphaofa[5],dalphaofa[5],lambda2; //EFTofDE
	switch(model) {
	  case 1:
		/* LCDM */
						return  0. ;
	  case 2:
		/* f(R) Hu-Sawicki */
						fr0 = extpars[0]/h0/h0;
						var1 = pow2(3.*omega0 -4.);

						term1 = 1./pow2(omega0/pow3(yenv * a) + 4. - 4.*omega0);
						term2 = 1./pow2(omega0/pow3(yh * a) + 4. - 4.*omega0);

						dod = fr0 * yh * a *var1 * (term1 - term2) / Rth/ Rth / omega0 ;
						dod2 = pow2(dod);

						dRRth = dod - dod2 + dod2*dod/3.;

						return gsl_min(dRRth,1./3.);

	   case 3:
	 	 /*DGP normal branch*/
		      betadgp = beta(a,omega0,extpars[0]);

		      delta = (1.+delta_initial)/pow3(yh) - 1.;

		      xterm = 9.*pow2(betadgp)*extpars[0]/(2.*omega0/pow3(a)*delta);

		      xm3 = 1./xterm;

		      return 2./(3.*betadgp)*(sqrt(1.+ xm3)-1.)/xm3;

			case 4:
			/* QUINTESSENCE */
							return  0. ;
			case 5:
			/* CPL */
							return  0. ;
			case 6:
			/* Dark Scattering with CPL  */
							return  0. ;

			case 7:
			/* EFTofDE - PPF */
			// p1-p7 : extpars[3-9] , extpars[0-2] : alpha_{K0}, alpha_{B0}, alpha_{M0}
					var1 = extpars[3]/(extpars[3]-1.)*extpars[5]; // a
					Mvir = pow3(Rth/0.1)*5.*omega0; // virial mass x Gnewton - see definition in scol_init in HALO.cpp

					delta = (1.+delta_initial)/pow3(yh) - 1.; // non-linear over density of halo
					deltaenv = (1.+delta_initial)/pow3(yenv) - 1.; // non-linear over density of environment

					var2 = pow(1./delta,1./3.); // yh with 1608.00522 definitions
					var3 = pow(1./deltaenv,1./3.);  // yenv with 1608.00522 definitions

					term1 = extpars[6]*pow(a,extpars[7])*pow(2.*h0*Mvir,extpars[8])*pow(var3/var2,extpars[9]); // y0

					xm3 = pow(term1/var2,var1); // (y0/yh)^a


							return extpars[3]*extpars[4]*(pow(1.+ xm3,1./extpars[3])-1.) / xm3 ;

			case 8:
			/* EFTofDE - unscreened approximation  */
							return mu(a,0.001,omega0,extpars,model)-1.;
			case 9:
			/* EFTofDE - superscreened approximation  */
							return 0.;

			case 10:
			/* KGB with CPL background */
							return 0.
		default:
					warning("SpecialFunctions: invalid model choice, model = %d \n", model);
				  return 0;
  }
}

// Dark energy contribution to virial theorem - see Eq.A9 of 1812.05594
double  WEFF(double a, double omega0, double extpars[], int model){
	double h2;
	switch(model) {
	 	 case 1:
		/* LCDM */
	 	 return 2.*(1.-omega0);

		 case 2:
		/* f(R) */
	 	 return 2.*(1.-omega0);

		 case 3:
		 /* DGP */
		 return 2.*(1.-omega0);

		 case 4:
		/* QUINTESSENCE*/
		return -(1.-omega0)*(1.+3.*extpars[0]);

		case 5:
		/*CPL*/
		 h2 = pow2(HAg(a,omega0,extpars,model));
		 return -(1.+3.*(extpars[0]+(1.-a)*extpars[1]))*(h2-omega0/pow3(a));

		 case 6:
		 /* Dark Scattering with CPL  */
		 h2 = pow2(HAg(a,omega0,extpars,model));
		 return -(1.+3.*(extpars[0]+(1.-a)*extpars[1]))*(h2-omega0/pow3(a));

		 case 7:
		 /* EFTofDE: PPF */
		 return 2.*(1.-omega0);

		 case 8:
		 /* EFTofDE: unscreened approximation */
		 return 2.*(1.-omega0);

		 case 9:
		 /* EFTofDE: superscreened approximation */
		 return 2.*(1.-omega0);

		 case 10:
		 /* KGB with CPL background */
		 h2 = pow2(HAg(a,omega0,extpars,model));
		 return -(1.+3.*(extpars[0]+(1.-a)*extpars[1]))*(h2-omega0/pow3(a));


		 default:
				 warning("SpecialFunctions: invalid model choice, model = %d \n", model);
				 return 0;
}
}




////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
                  /////////////////////////////////////////////////
                         //////////////////////////////////
                                  //////////////////
                                        /////////

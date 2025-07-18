#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cfloat>
#include <cmath>

#include "Common.h"
#include "Quadrature.h"
#include "BeyondLCDM.h"
#include "SpecialFunctions.h"
#include "Spline.h"
#include "SCOL.h"

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

/* See 1606.02520, 1608.00522, 1808.01120, 2005.12184  for more details */

// "model" selects MG or DE model:
// 1: GR
// 2: Hu-Sawicki f(R) [LCDM background]
// 3: normal branch DGP [LCDM background]
// 4: quintessence
// 5: CPL
// 6: CPL with dark scattering
// 7: EFTofDE with mu in k->infinity limit & nPPF (1608.00522) G_eff for non-linear scales [assumes gamma2 = gamma3 = 0]
// 8: EFTofDE with mu in k->infinity limit & linear G_eff for non-linear scales [assumes gamma2 = gamma3 = 0]
// 9: EFTofDE with mu in k->infinity limit & G_Newton for non-linear scales [assumes gamma2 = gamma3 = 0]
// 10: EFTofDE with full mu & nPPF (1608.00522) G_eff for non-linear scales [assumes gamma2 = gamma3 = 0]
// 11: EFTofDE with full mu & linear G_eff for non-linear scales [assumes gamma2 = gamma3 = 0]
// 12: EFTofDE with full mu & Phenomenological G_eff for non-linear scales [assumes gamma2 = gamma3 = 0]
// 13: Model independent parametrisation: CPL {w0,wa} for background, growth index gamma for linear perturbations, Phenomenological G_eff for nonlinear scales
// 14: Cubic Galileon
// 15: QCDM
// 16: K-mouflage as in 1403.5424
// 17: K-mouflage with nPPF [1608.00522]
// 18: Model independent parametrisation: CPL {w0,wa} for background, mu(k,z) for linear perturbations (see https://arxiv.org/pdf/2411.12026 for default form), Phenomenological G_eff for nonlinear scales (q1,q2,q3,q4) 


// Throughout the values of pars, extpars and model are:

// "pars[]" : base parameters. Currently used:
// 0: scale factor
// 1: total matter fraction today
// 2: total massive neutrino fraction today

// extpars[maxpars]: extended model parameters  (maxpars = 20 default, see BeyondLCDM.h)

// 0: Omega_rc for nDGP, fr0 for f(R), w0 for CPL
// 1: wa for CPL
// 2: xi * h for Dark scattering

// for EFTofDE (models 7-11) we have:
// 0: alpha_K(a)
// 1: alpha_B(a)
// 2: alpha_M(a)
// 3: alpha_T(a)
// 4: M^2/M_planck^2

// for models 7 & 10:
// 5+ : nPPF parameters (see eq 5.3 of 1608.00522)

// for model 12:
// 5: par1 : size of screening scale
// 6: mass dependence of screening scale
// 7: y_environment dependence of screening scale
// 8: Yukawa suppression scale mass dependence

// for model 13:
// Background CPL:
// 0: w0
// 1: wa

// linear perturbations, growth index :
// 2: gamma

// nonlinear, error function parametrisation:
// 3: size of screening scale
// 4: mass dependence of screening scale
// 5: y_environment dependence of screening scale
// 6: Yukawa suppression scale mass dependence


// for model 14-15 CG and QCDM:
// 0 : s
// 1 : q

// for model 16: K-mouflage
// 0: n
// 1: lambda
// 2: K0
// 3: beta0


// for model 17: K-mouflage with nPPF
// 0: n
// 1: lambda
// 2: K0
// 3: beta0
// 5+ - nPPF params

// for model 18: Model independent parametrisation with CPL background
// See Eq. 3.6 of https://arxiv.org/pdf/2411.12026 for linear Geff(k,z) parametrisation and Eq.41 of https://arxiv.org/pdf/2210.01094 nonlinear Geff(k,z) 
// 0: w0
// 1: wa
// 2: mu0 
// 3: c1 
// 4: lambda 
// 5: size of screening scale
// 6: mass dependence of screening scale
// 7: y_environment dependence of screening scale
// 8: Yukawa suppression scale mass dependence
 

/* A) Functions to edit for new models in Reaction computation: */

// HAg : Normalised Hubble expansion H/H_0
// HA1g : Normalised scale factor derivative of Hubble  aH dH/da / H0^2
// myfricF : Euler equation friction term (only present in Dark Scattering)
// mu : linear modification to the Poisson equation
// gamma2 : 2nd order modification to the Poisson equation
// gamma3 : 3rd order modification to the Poisson equation

// mymgF : nonlinear modification to the Poisson equation [Needed only for spherical collapse SCOL.cpp, not for 1-loop calculations ]
// WEFF : Effective dark energy fluid contribution to virial theorem [Needed only for spherical collapse SCOL.cpp, not for 1-loop calculations ]


/* B) Functions to edit for EFTofDE parametrisations needed for mu (cases 10-12): */

// riccibackgroundp : scale factor derivative of background Ricci scalar in FRLW
// alphai_eft : scale factor dependence of parameter : alpha_i(a)
// dalphai_eft : scale factor derivative :  d alpha_i(a)/ da
// ddalphai_eft : 2nd scale factor derivative : d^2 alpha_i(a)/ da^2
// HA3g : normalised 2nd scale factor derivative of Hubble :  d^2 H(a) / da^2  /H0


/* C) Functions to edit for custom background expansion: */
/* Note: to use these functions in your new model, you should initialise a spline with H(a) - run hubble_init : see SpecialFunctions.cpp */

// bespokehub : Normalised Hubble expansion: H(a)/H0 (this can involve a solution to some ODE)
// bespokehubd : aH dH/da / H0^2
// bespokehubdd : d^2 H(a) / da^2  /H0

// Once specified, you can use the following splines (initialised with hubble_init in SpecialFunctions class) in HAg and HA1g (and HA3g if EFTofDE is being considered)
// myhubble : H/H0
// myhubbled : aH dH/da / H0^2
// myhubbledd : d^2 H(a) / da^2  /H0

                                        ////
                                    /////////////
                            //////////////////////////////
                      /////////////////////////////////////////
                ///////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

// Some useful functions needed for EFTofDE computation of linear G_eff (mu-1)

// Ricci background function - see notebooks/GtoPT.nb
double riccibackground(double a, double omega0, double extpars[], int model){
	double h0 = 1./2997.92458;
	double a3 = pow3(a);
	// solutions currently assume LCDM
	switch(model) {
		case 10:
		/* EFTofDE with nPPF and full k-dependence in linear modification */
		return  3. * pow2(h0) * (omega0 + 4. * a3 * (1. - omega0)) / a3;

		case 11:
		/* EFTofDE unscreened and full k-dependence in linear modification */
		return  3. * pow2(h0) * (omega0 + 4. * a3 * (1. - omega0)) / a3;

		case 12:
		/* EFTofDE unscreened and full k-dependence in linear modification */
		return  3. * pow2(h0) * (omega0 + 4. * a3 * (1. - omega0)) / a3;


		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
}
}

// Ricci background scale factor derivtative - see notebooks/GtoPT.nb
double riccibackgroundp(double a, double omega0, double extpars[], int model){
		// solutions currently assume LCDM
		double h0 = 1./2997.92458;
		double a4 = pow(a,4);
		switch(model) {
			case 10:
			/* EFTofDE with nPPF and full k-dependence in linear modification */
			return  -9.*pow2(h0)*omega0/a4;

			case 11:
			/* EFTofDE unscreened and full k-dependence in linear modification */
			return  -9.*pow2(h0)*omega0/a4;

			case 12:
			/* EFTofDE unscreened and full k-dependence in linear modification */
			return  -9.*pow2(h0)*omega0/a4;

			default:
					warning("BeyondLCDM: invalid model choice, model = %d \n", model);
					return 0;
	}
}

// Ricci background 2nd scale factor derivtative
double riccibackgroundpp(double a, double omega0, double extpars[]){
	double h0 = 1./2997.92458;
	double a5 = pow(a,5);
	return 36.*pow2(h0)*omega0/a5;
}


/* Hu-Sawicki f(R) function (we assume Ricci as given by model 11 - was used for tests)*/
double fofR_hs(double a, double omega0, double fr0){
	double extra[20];
	double ricci = riccibackground(a,omega0,extra,11); // background Ricci
	double ricci0 = riccibackground(1.,omega0,extra,11); // background Ricci today
	return fr0*pow2(ricci0/ricci);
}

// f_R'
double fofRd_hs(double a, double omega0, double fr0){
	double extra[20];
	double ricci = riccibackground(a,omega0,extra,11); // background Ricci
	double riccid = riccibackgroundp(a,omega0,extra,11); // background Ricci sf derivative
	return -2.*riccid/ricci * fofR_hs(a,omega0,fr0);
}

// f_R''
double fofRdd_hs(double a, double omega0, double fr0){
	double extra[20];
	double ricci = riccibackground(a,omega0,extra,11); // background Ricci
	double riccid = riccibackgroundp(a,omega0,extra,11); // background Ricci sf derivative
	double riccidd = riccibackgroundpp(a,omega0,extra); // background Ricci 2nd sf derivative
	return 2.*fofR_hs(a,omega0,fr0) * (3.*pow2(riccid/ricci) - riccidd/ricci);
}

// alpha_M(a) in Hu-Sawicki theory (see GtoPT.nb)
double alpha_m_hs(double a, double omega0, double fr0){
	double fR = fofR_hs(a,omega0,fr0);
	double fRp = fofRd_hs(a,omega0,fr0);
	return a*fRp/(1.+fR);
}

// alpha_M'(a) in Hu-Sawicki theory (see GtoPT.nb)
double alpha_md_hs(double a, double omega0, double fr0){
	double fR = fofR_hs(a,omega0,fr0);
	double fRp = fofRd_hs(a,omega0,fr0);
	double fRpp = fofRdd_hs(a,omega0,fr0);
	return (fRp*(1.+fR - a*fRp)+a*(1.+fR)*fRpp)/pow2(1.+fR);
}


/* BACKGROUND QUANTITIES */

/*  Scale factor evolution of alpha_i for EFTofDE */
// can add in different scale factor dependencies for different alphas
// edit as necessary - but also edit the analytic scale factor derivative function [dalphai_eft] accordingly

// If Hubble is specified, then we should define M^2 as the integral of alpha_M

// integral to calculate time dependence of variation to planck mass from alpham generically
inline double M2integral(double alpha0, double omega0, double a){
	double alpham = alpha0 * a; // insert consistent expression for alpha_M as appearing in case 3 below
	return alpham / a ;
}

inline double alphai_eft(double a, double omega0, double alpha0, int model){
	switch(model) {
		case 1:
				return alpha0*a; // alpha_K(a)
		case 2:
				return alpha0*a; // alpha_B(a)
			//	 return -alpha_m_hs(a,omega0,alpha0); // Hu-Sawicki form
		case 3:
				return alpha0*a; // alpha_M(a)
			//	 return alpha_m_hs(a,omega0,alpha0); // Hu-Sawicki form
		case 4:
				return alpha0*a; // alpha_T(a)
		case 5:
				  return exp(alpha0*a); // M2(a)/m_planck = e^{Integrate[alpha_M / a]}  [We define alpha_M = a (dM^2/da) / M^2]
		//		return Integrate(bind(M2integral, alpha0, omega0, std::placeholders::_1), AMIN , a, 1e-3); // M2(a)/m_planck generic integral of alpha_M/a
		//		return (1.+fofR_hs(a,omega0,alpha0)); // Hu-Sawicki form
		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
			}
}


// scale factor derivatives of alphai_eft
inline double dalphai_eft(double a, double omega0, double alpha0, int model){
	switch(model) {
		case 1:
				return alpha0; // alpha_K'(a)
		case 2:
				return alpha0; // alpha_B'(a)
			//	return -alpha_md_hs(a,omega0,alpha0);  // Hu-Sawicki form
		case 3:
				return alpha0; // alpha_M'(a)
			//	return alpha_md_hs(a,omega0,alpha0); // Hu-Sawicki form
		case 4:
				return alpha0; // alpha_T'(a)
		case 5:
				return alpha0*exp(alpha0*a); // doesn't appear in equations (this is directly related to alpha_m)
		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
		}
}

// 2nd  scale factor derivatives of alphai_eft
inline double ddalphai_eft(double a, double omega0, double alpha0, int model){
	switch(model) {
		case 1:
				return 0.; // alpha_K''(a)
		case 2:
				return 0.; // alpha_B''(a)
		case 3:
				return 0.; // alpha_M''(a)
		case 4:
				return 0.; // alpha_T''(a)
		case 5:
				return 0.; // doesn't appear in equations
		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
		}
}



/* K-mouflage background equations - see 3.58,59,61 of 2110.00566, with various typos corrected (see notebooks/K-mouflage.nb) - derivative with respect to N=d/dlna  */

// Exponential conformal factor as in Eq. 2.25
double conf_fac(double phi, double beta0){
	return exp(phi*beta0);
//	return (1.+beta0*phi); // approximate exponential choice
}

/* Parameters passed to system of equations*/
// omega0 = Omega_{m,0}
// omeganu = Omega_{nu,0}
// maxpars - maximum number of extended parameters - specified in SpecialFunctions.h

struct param_type_km {
  double omega0;
	double omeganu;
	double extpars[maxpars];
};


// ODE solver for scalar field and its derivative - Goal is to create splined functions
Spline kmouflage_phi;
Spline kmouflage_phid;


//Jacobian of the system required when calling the system evolver, the below is not needed for solving
int jac_sub (double a, const double G[], double *dfdy, double dfdt[], void *params)
{
	return GSL_SUCCESS;
}

// analytic solution to Eq. 3.58 for n=2 (see k-mouflage.nb)
double kmouflage_hub2_n2(double a, double phid, double Omegam0, double lambda, double K0, double Aphi){
	double a2 = pow2(a);
	double a3 = a*a2;
	double lam2 = pow2(lambda);
	double phid2 = pow(phid,2.);
	double phid4 = phid2 * phid2;

	double denom = 3.*K0*phid4;
	// including radiation from glam = ~9e-5
	double Orad = 0.;
	double term1 = 36.*K0/(lam2*a3)*(Aphi*Omegam0 + Orad/a);

	double root = sqrt( 36. -12.*phid2 + phid4*(1. - 12.*K0 - term1) ) ;
	double term2 = -6. + phid2;

	return -a2 * lam2 * (term2 + root) / denom;
}


// Differential equations - all derivatives are with respect to ln[a]
int funcn_kmouflage_lna(double lna, const double G[], double F[], void *params)
{
  param_type_km p = *(param_type_km *)(params);

	// Cold dark matter + baryon fraction
	double omegacb = p.omega0 - p.omeganu;

	// K-mouflage parameters
	double n = p.extpars[0];
	double lambda = p.extpars[1];
	double K0 = p.extpars[2];
	double beta0 = p.extpars[3];

	double lam2 = pow2(lambda);

	// scale factor and square
	double a = exp(lna);
	double a2 = pow2(a);

// Legend:
	// G[0] = phi
	// G[1] = phi'
	// F[0] = phi'
	// F[1] = phi''

	// ' is a derivative with respect to ln[a]

	F[0] = G[1];

	// Background quantities:
	// Conformal factor
	double Aphi = conf_fac(G[0],beta0) ;
	// Conformal factor derivative wrt scalar field
	double Aphid = beta0*Aphi;

	// phi'^2
	double phid2 = pow2(G[1]);

	// Eq. 3.58 Friedmann equation: (a H/H0)^2, i.e. normalised conformal Hubble factor
	double myhub2 = kmouflage_hub2_n2(a, G[1], omegacb, lambda, K0, Aphi);

	// Eq. 3.59 - Raychaudhuri equation for given phi and phi' and H/H0: d (aH/H0)/ d tau
	// including radiation from glam = ~9e-5
	double Orad = 0.;
	double myhubd = - omegacb/2. * Aphi / a - Orad/a2
									+ 1./3.*a2*lam2 // Typo in 2110.00566 : 1/3 not 2/3
									- 1./3.*myhub2*phid2 // Typo in 2110.00566 : 1/3 not 2/3
									- 1./3.*(n + 1.)*K0*a2*lam2 * pow( myhub2 * phid2 / (2.*lam2*a2)  , n) ; // Typo in 2110.00566 : 1/6 not 1/3

	// Canonical kinetic term
	double X = myhub2 / (2.*lam2*a2) * phid2;
	// kinetic term derivative
	double KX = 1. + n * K0 * pow(X,n-1.);
	// kinetic term second derivative
	double KXX = n * (n-1.) * K0 * pow(X,n-2.);

	// Eq.3.61 - Klein-Gordon equation : phi''
	double kfactor = (KX + 2.*X*KXX);
	F[1] = ( (-2.*(KX - X * KXX)*myhub2*G[1] - 3.*Aphid * omegacb / a )/kfactor - myhubd*G[1] )/myhub2;

	return GSL_SUCCESS;
}

void init_kmouflage_lna(double pars[], double extpars[], int Na)
{
				double omega0 = pars[1]; // total matter fraction
				double omeganu = pars[2]; // massive neutrino fraction

			// Initial conditions for scalar field and first scale factor derivative (MG-GLAM set both to 1e-30 but this gives strange results )
				double phii =  1e-30;
				double phidi = -1e-6;

			// Our time variable is ln[a]
				double lnavals[Na];
			// Scale factor for final spline
				double avals[Na];

				double epsilon = 1e-4;
				double amax = 1.+0.5; // Go higher than 1. to have full coverage - much higher for Jordan frame transforms
				double amin = AMIN*(1.-epsilon); // Go a bit lower than AMIN to have full coverage

					for(int i = 0; i< Na; i++){
						 avals[i] = amin * exp(i*log(amax/amin)/(Na-1.));
						 lnavals[i] = log(avals[i]);
					}

        // Non-Eds ICs
			  double G[2] = {phii,phidi};

			/*Parameters passed to system of equations */
				struct param_type_km mypars;

					// set all parameters
					// total matter
					mypars.omega0 = omega0;
					// massive neutrinos
				  mypars.omeganu = omeganu;
					// model parameters
					for (int i=0; i<maxpars;i++){
						mypars.extpars[i] = 	extpars[i];
					}


				gsl_odeiv2_system sys = {funcn_kmouflage_lna, jac_sub, 2, &mypars};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-8, 1e-8, 1e-8);

				// evolve everything to redshift z=z_max and start sampling
				int status1 = gsl_odeiv2_driver_apply (d, &lnavals[0], lnavals[1], G);

				double phi[Na];
				double phid[Na];

				/*Allocation of array values */
				phi[0] = phii;
				phid[0] = phidi;
				phi[1] = G[0];
				phid[1] = G[1];

					for(int i = 2; i < Na; i++){
							status1 = gsl_odeiv2_driver_apply (d, &lnavals[i-1], lnavals[i], G);
							/*Allocation of array values */
							phi[i] = G[0];
							phid[i] = G[1];
						}
						// free memory
						gsl_odeiv2_driver_free(d);


// Create Splines for phi and its ln[scale factor] derivative as functions of the scale factor
 kmouflage_phi = LinearSpline(Na, avals , phi);
 kmouflage_phid = LinearSpline(Na, avals , phid);

}

// Planck mass normalised K-mouflage scalar field
double kmouflage_sf(double a){
	return kmouflage_phi(a);
}
// Its d/dlna derivative
double kmouflage_sf_der(double a){
	return kmouflage_phid(a);
}



/* Background normalised hubble function: general solution for initialisation with hub_init in SpecialFunctions.cpp */
// Useful if H(a) does not have a straight forward analytic form
// e.g. if we need to solve some equation for H(a):
// we will create a splined function using hub_init instead of solving the equation at each call
double bespokehub(double a, double omega0, double extpars[], int model){
	double counter, limit1, limit2, flimit1, flimit2, dmean, soluction, fsoluction, bolean;
	double omegaL, Aphi;
	/* Insert solver */
	switch(model) {
		case 10:
		/* EFTofDE with nPPF and full k-dependence in linear modification */
		//can change to some solution for H in terms of extra parameters
			return HA(a,omega0);

		case 11:
		/* EFTofDE unscreened and full k-dependence in linear modification */
		//can change to some solution for H in terms of extra parameters
			return HA(a,omega0);

		case 14:
		/* Generalised Cubic Galileon background solver */
			  omegaL = 1.-omega0;
				// Give range for solution
				limit1=0.;
				limit2=pow(10,9);

				// Set extreme values
				flimit1=omegaL+omega0/pow3(a)*pow(limit1,extpars[0])-pow(limit1,2.+extpars[0]);
				flimit2=omegaL+omega0/pow3(a)*pow(limit2,extpars[0])-pow(limit2,2.+extpars[0]);

				dmean=(limit2-limit1)/2.;
				soluction=limit2-dmean;
				fsoluction=1.;
				counter=0;
				while (sqrt(fsoluction*fsoluction)>1e-5 && counter < 1e5) {
					fsoluction=omegaL+omega0/pow3(a)*pow(soluction,extpars[0])-pow(soluction,2.+extpars[0]);
					bolean=fsoluction*flimit1;
						 if (bolean > 0.){
									 limit1=soluction;
									 flimit1=fsoluction;
							 }
							 if (bolean < 0.){
									 limit2=soluction;
									 flimit2=fsoluction;
							 }
							 dmean=(limit2-limit1)/2.;
							 soluction=limit1+dmean;
							 counter=counter+1.;
					}
					return soluction;


			case 15:
			/* QCDM */
			omegaL = 1.-omega0;
				limit1=0;
				if (limit1<0){
					limit1=0;
				}
				limit2=pow(10,9);
				flimit1=omegaL+omega0/pow3(a)*pow(limit1,extpars[0])-pow(limit1,2.+extpars[0]);
				flimit2=omegaL+omega0/pow3(a)*pow(limit2,extpars[0])-pow(limit2,2.+extpars[0]);
				dmean=(limit2-limit1)/2.;
				soluction=limit2-dmean;
				fsoluction=1;
				counter=0;
				while (sqrt(fsoluction*fsoluction)>1e-5 && counter < 1e5) {
					fsoluction=omegaL+omega0/pow3(a)*pow(soluction,extpars[0])-pow(soluction,2.+extpars[0]);
					bolean=fsoluction*flimit1;
						 if (bolean > 0.){
									 limit1=soluction;
									 flimit1=fsoluction;
							 }
							 if (bolean < 0.){
									 limit2=soluction;
									 flimit2=fsoluction;
							 }
							 dmean=(limit2-limit1)/2;
							 soluction=limit1+dmean;
							 counter=counter+1;
					}
					return soluction;

			case 16:
			/* K-mouflage */
				Aphi = conf_fac(kmouflage_phi(a),extpars[3]) ;
				return sqrt(kmouflage_hub2_n2(a, kmouflage_phid(a), omega0, extpars[1], extpars[2], Aphi)/pow2(a));

			case 17:
			/* K-mouflage with nPPF */
				Aphi = conf_fac(kmouflage_phi(a),extpars[3]) ;
				return sqrt(kmouflage_hub2_n2(a, kmouflage_phid(a), omega0, extpars[1], extpars[2], Aphi)/pow2(a));


		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
}
}

//  (dH/dt / H0^2) = aH dH/da / H0^2
// HA1 is the LCDM solution (see SpecialFunctions.cpp)
double bespokehubd(double a, double omega0, double extpars[], int model){
	/* Insert solver */
	switch(model) {
		case 10:
		/* EFTofDE with nPPF and full k-dependence in linear modification */
		//can change to some solution for H in terms of extra parameters
			return HA1(a,omega0);

		case 11:
		/* EFTofDE unscreened and full k-dependence in linear modification */
		//can change to some solution for H in terms of extra parameters
		 return HA1(a,omega0);

	 case 14:
		/* CG  */
			return 0; // we use analytic form - see below

		case 15:
		/* QCDM  */
			return 0; // we use analytic form - see below

		case 16:
  		/* K-mouflage  */
 			return 0; // we use analytic form - see below

		case 17:
  		/* K-mouflage with nPPF  */
 			return 0; // we use analytic form - see below


		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
}
}

//d^2 H(a) / da^2  /H0
// current solutions assume LCDM
double bespokehubdd(double a, double omega0, double extpars[], int model){
	double hub;
	switch(model) {
		case 10:
		/* EFTofDE with nPPF and full k-dependence in linear modification */
		   hub = HA(a,omega0); // H / H0
			return 3.*omega0*(-12.*(omega0-1.)/(4.*pow2(hub)) + 5.)/(4.*pow(a,5)*hub) ;

		case 11:
		/* EFTofDE unscreened and full k-dependence in linear modification */
			hub = HA(a,omega0); // H / H0
	 		return 3.*omega0*(-12.*(omega0-1.)/(4.*pow2(hub)) + 5.)/(4.*pow(a,5)*hub) ;

		case 14:
	 		return 0; // we use analytic form - see below

		case 15:
		 	return 0; // we use analytic form - see below

		case 16:
			return 0; // we use analytic form - see below

		case 17:
		 return 0; // we use analytic form - see below

		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
}
}


// H(a)/H0
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
		/* EFTofDE: nPPF - set to LCDM background for now  */
		return HA( a, omega0);

		case 8:
		/* EFTofDE: unscreened approx */
		return  HA( a, omega0);

		case 9:
		/* EFTofDE: superscreened approx */
		return  HA( a, omega0);

		case 10:
		/* EFTofDE with nPPF and full k-dependence in linear modification */
		return HA( a, omega0); //myhubble(a);

		case 11:
		/* EFTofDE: unscreened approx and full k-dependence in linear modification */
		return  HA( a, omega0);

		case 12:
		/* EFTofDE: Phenomenological non-linear model */
		return  HA( a, omega0);

		case 13:
		/* Minimal model independent : CPL + Gamma + Erf */
		A = -3.*(1.+extpars[0]+extpars[1]);
		omegaf = pow(a,A)*exp(3.*(-1.+a)*extpars[1]);
		omegaL= (1.-omega0)*omegaf;
		return  sqrt(omega0/pow(a,3)+omegaL);


		case 14:
		/* CG */
		return myhubble(a);

		case 15:
		/* QCDM */
		return myhubble(a);

		case 16:
		/* K-mouflage */
		return myhubble(a);

		case 17:
		/* K-mouflage with nPPF*/
		return myhubble(a);

		case 18:
		/*  Model independent Poisson equation parametrisation as in https://arxiv.org/pdf/2411.12026 with CPL background  */
		A = -3.*(1.+extpars[0]+extpars[1]);
		omegaf = pow(a,A)*exp(3.*(-1.+a)*extpars[1]);
		omegaL= (1.-omega0)*omegaf;
		return  sqrt(omega0/pow(a,3)+omegaL);

			default:
					warning("BeyondLCDM: invalid model choice, model = %d \n", model);
					return 0;
}
}

/* Normalised time derivative  */


//  (dH/dt / H0^2) = aH dH/da / H0^2
double HA1g(double a, double omega0, double extpars[], int model){
	double A, omegaf, omegaL,h2, myhub, Aphi, a2, a4, var1, Orad;
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
		/* EFTofDE: nPPF - set to LCDM background for now  */
		return HA1( a, omega0);

		case 8:
		/* EFTofDE: unscreened approximation */
		return HA1( a, omega0);

		case 9:
		/* EFTofDE: superscreened approximation */
		return HA1( a, omega0);

		case 10:
		/* EFTofDE with nPPF and full k-dependence in linear modification */
		return HA1( a, omega0);//myhubbled(a);

		case 11:
		/* EFTofDE unscreened and full k-dependence in linear modification */
		return  HA1( a, omega0);

		case 12:
		/* EFTofDE: Phenomenological non-linear model */
		return  HA1( a, omega0);

		case 13:
		/* Minimal model independent : CPL + Gamma + Erf */
		A = -3.*(1.+extpars[0]+extpars[1]);
		omegaf = pow(a,A)*exp(3*(-1.+a)*extpars[1]);
		omegaL= (1.-omega0)*omegaf;
		return -3.*omega0/(2.*pow(a,3)) + (A+3.*a*extpars[1])*omegaL/2.;

		case 14:
		/* CG */
		h2 = pow2(HAg(a,omega0,extpars,model));
		return -h2*(omega0/pow(a,3)*(-3.-extpars[0])+(extpars[0]+2.)*h2)/(extpars[0]*omega0/pow(a,3)-(extpars[0]+2.)*h2)-h2;

		case 15:
		/* QCDM */
		h2 = pow2(HAg(a,omega0,extpars,model));
		return -h2*(omega0/pow(a,3)*(-3.-extpars[0])+(extpars[0]+2.)*h2)/(extpars[0]*omega0/pow(a,3)-(extpars[0]+2.)*h2)-h2;


		case 16:
		/* K-mouflage */
		myhub = HAg(a,omega0,extpars,model);
		// conformal hubble
		h2 = pow2(a*myhub);
		Aphi = conf_fac(kmouflage_phi(a),extpars[3]) ;
		a2 = pow2(a);
		a4 = pow2(a2);
		// conformal Hubble conformal time derivative - See Eq. 3.59 of 2110.00566
		Orad = 0.;//9e-5;

		var1 = - omega0/2. * Aphi / a - Orad/a2
					+ 1./3.*a2*pow2(extpars[1])
					- 1./3.*h2*pow2(kmouflage_phid(a))
					- 1./3.*(extpars[0] + 1.)*extpars[2]*a2*pow2(extpars[1]) * pow(h2*pow2(kmouflage_phid(a))/(2.*pow2(extpars[1])*a2) , extpars[0]);

		return 1./ (a2) * (var1 - h2);

		case 17:
		/* K-mouflage with nPPF */
		myhub = HAg(a,omega0,extpars,model);
		// conformal hubble
		h2 = pow2(a*myhub);
		Aphi = conf_fac(kmouflage_phi(a),extpars[3]) ;
		a2 = pow2(a);
		a4 = pow2(a2);
		// conformal Hubble conformal time derivative - See Eq. 3.59 of 2110.00566
		var1 = - omega0/2. * Aphi / a
					+ 1./3.*a2*pow2(extpars[1])
					- 1./3.*h2*pow2(kmouflage_phid(a))
					- 1./3.*(extpars[0] + 1.)*extpars[2]*a2*pow2(extpars[1]) * pow(h2*pow2(kmouflage_phid(a))/(2.*pow2(extpars[1])*a2) , extpars[0]);

		return 1./ (a2) * (var1 - h2);

		case 18:
		/*  Model independent Poisson equation parametrisation as in https://arxiv.org/pdf/2411.12026 with CPL background  */
		A = -3.*(1.+extpars[0]+extpars[1]);
		omegaf = pow(a,A)*exp(3*(-1.+a)*extpars[1]);
		  omegaL= (1.-omega0)*omegaf;
		return -3.*omega0/(2.*pow(a,3)) + (A+3.*a*extpars[1])*omegaL/2.;
  

		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
}
}

/* Normalised 2nd time derivative  */

//d^2 H(a) / da^2  /H0
double HA3g(double a, double omega0, double extpars[], int model){
	double hub;
	switch(model) {
		case 10:
		/* EFTofDE with nPPF and full k-dependence in linear modification */
		hub = HA(a,omega0); // H / H0
		return 3.*omega0*(-12.*(omega0-1.)/(4.*pow2(hub)) + 5.)/(4.*pow(a,5)*hub) ;

		case 11:
		/* EFTofDE unscreened and full k-dependence in linear modification */
		hub = HA(a,omega0); // H / H0
	  return 3.*omega0*(-12.*(omega0-1.)/(4.*pow2(hub)) + 5.)/(4.*pow(a,5)*hub) ;

		case 12:
		/* EFTofDE Phenomenological model and full k-dependence in linear modification */
		hub = HA(a,omega0); // H / H0
	  return 3.*omega0*(-12.*(omega0-1.)/(4.*pow2(hub)) + 5.)/(4.*pow(a,5)*hub) ;

		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
}
}


// some useful functions of Hubble

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
		/* EFTofDE: nPPF */
					return 0.;
		case 8:
		/* EFTofDE: unscreened approx */
					return 0.;
		case 9:
		/* EFTofDE: superscreened approx */
					return 0.;
		case 10:
		/* EFTofDE with nPPF and full k-dependence in linear modification */
					return 0.;
		case 11:
		/* EFTofDE: unscreened approx and full k-dependence in linear modification */
					return 0.;

		case 12:
		/* EFTofDE: Phenomenological non-linear model */
		return  0.;

		case 13:
		/* Minimal model independent : CPL + Gamma + Erf */
		return  0.;

		case 14:
		/* CG */
		return  0.;

		case 15:
		/* QCDM */
		return  0.;

		case 16:
		/* K-mouflage */
		return   extpars[3] * kmouflage_phid(a) * HAg(a,omega0,extpars,model);

		case 17:
		/* K-mouflage with nPPF */
		return   extpars[3] * kmouflage_phid(a) * HAg(a,omega0,extpars,model);

		case 18:
		/*  Model independent Poisson equation parametrisation as in https://arxiv.org/pdf/2411.12026 with CPL background  */
		return 0.; 

		default:
					warning("BeyondLCDM: invalid model choice, model = %d \n", model);
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
	double var1, var2, var3, alphaofa[5],dalphaofa[5],ddalphaofa[5];
	double myA[3],myB[3],myC[4],myf[4];
	double hub, hubd, hubdd; // hubble and its derivatives

	// parameters for cases 7-12:
	double ct2; // speed of grav waves
	double ca; // background function
	double R0d; // background Ricci scale factor derivative
	double betaxi;
	double tol = 1e-15; // tolerance to avoid singularity in cases 7-12
	double qsa_test; // condition given in Eq.39 of 1712.00444 for QSA to hold
	double aqsa = 0.01; // scale factor at which we want condition to hold
	double kqsa = 0.01; // scale at which we want condition to hold

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
		// 7-9: nPPF, unscreened, superscreened assuming LCDM background and k->infinity limit
		case 7:
						// 0: alpha_K(a)
						// 1: alpha_B(a)
						// 2: alpha_M(a)
						// 3: alpha_T(a)
						// 4: M^2/M_planck^2

						// Initialise alpha_i (a) and their derivatives
						for(int i=0; i<5; i++){
							alphaofa[i] = alphai_eft(a,omega0,extpars[i],i+1);
							dalphaofa[i] = dalphai_eft(a,omega0,extpars[i],i+1);
							ddalphaofa[i] = ddalphai_eft(a,omega0,extpars[i],i+1);
						}
						ct2 = 1. + alphaofa[3];

					  var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha

						// protect against divergences when var1=0
						if (fabs(var1) <tol) {
							var1 = tol;
						}

						var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[1]/2.*ct2 + HA2g(a,omega0,extpars,model)
																					 + alphaofa[2] - alphaofa[3])
																					 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)/alphaofa[4]); //cs^2

						 // protect against divergences when var2=0
 						if (fabs(var2) <tol) {
 							var2 = tol;
 						}

						 // QSA check for late times
 						qsa_test = k0/(HAg(a,omega0,extpars,model)*h0*a) -var2;
 						if (qsa_test<0 && a>aqsa && k0>kqsa) {
							warning("BeyondLCDM: QSA violated at a: %e and k: %e \n", a,k0);
 						}

						betaxi = 2./(var2*var1) * pow2(alphaofa[1]/2. * ct2 + alphaofa[2] - alphaofa[3]);


					  return (1. + alphaofa[3] + betaxi) / (alphaofa[4]);

	  case 8:
						// 0: alpha_K(a)
						// 1: alpha_B(a)
						// 2: alpha_M(a)
						// 3: alpha_T(a)
						// 4: M^2/M_planck^2
						// Initialise alpha_i (a) and their derivatives
						for(int i=0; i<5; i++){
							alphaofa[i] = alphai_eft(a,omega0,extpars[i],i+1);
							dalphaofa[i] = dalphai_eft(a,omega0,extpars[i],i+1);
							ddalphaofa[i] = ddalphai_eft(a,omega0,extpars[i],i+1);
						}

						ct2 = 1. + alphaofa[3];

						var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha

						// protect against divergences when var1=0
						if (fabs(var1) <tol) {
							var1 = tol;
						}

						var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[1]/2.*ct2 + HA2g(a,omega0,extpars,model)
																					 + alphaofa[2] - alphaofa[3])
																					 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)/alphaofa[4]); //cs^2

						 // protect against divergences when var2=0
 						if (fabs(var2) <tol) {
 							var2 = tol;
 						}

						 // QSA check for late times
 						qsa_test = k0/(HAg(a,omega0,extpars,model)*h0*a) -var2;
						if (qsa_test<0 && a>aqsa && k0>kqsa) {
 							warning("BeyondLCDM: QSA violated at a: %e and k: %e \n", a,k0);
 						}

						betaxi = 2./(var2*var1) * pow2(alphaofa[1]/2. * ct2 + alphaofa[2] - alphaofa[3]);

						return (1. + alphaofa[3] + betaxi) / (alphaofa[4]);

		case 9:
						// 0: alpha_K(a)
						// 1: alpha_B(a)
						// 2: alpha_M(a)
						// 3: alpha_T(a)
						// 4: M^2/M_planck^2

						// Initialise alpha_i (a) and their derivatives
						for(int i=0; i<5; i++){
							alphaofa[i] = alphai_eft(a,omega0,extpars[i],i+1);
							dalphaofa[i] = dalphai_eft(a,omega0,extpars[i],i+1);
							ddalphaofa[i] = ddalphai_eft(a,omega0,extpars[i],i+1);
						}

						ct2 = 1. + alphaofa[3];

						var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha

						// protect against divergences when var1=0
						if (fabs(var1) <tol) {
							var1 = tol;
						}

						var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[1]/2.*ct2 + HA2g(a,omega0,extpars,model)
																					 + alphaofa[2] - alphaofa[3])
																					 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)/alphaofa[4]); //cs^2

						 // protect against divergences when var2=0
 						if (fabs(var2) <tol) {
 							var2 = tol;
 						}

						 // QSA check for late times
 						qsa_test = k0/(HAg(a,omega0,extpars,model)*h0*a) -var2;
						if (qsa_test<0 && a>aqsa && k0>kqsa) {
							warning("BeyondLCDM: QSA violated at a: %e and k: %e \n", a,k0);
 						}

						betaxi = 2./(var2*var1) * pow2(alphaofa[1]/2. * ct2 + alphaofa[2] - alphaofa[3]);

						return (1. + alphaofa[3] + betaxi) / (alphaofa[4]);

		case 10:
					/* EFTofDE with nPPF and full k-dependence in linear modification */
						// 0: alpha_K(a)
						// 1: alpha_B(a)
						// 2: alpha_M(a)
						// 3: alpha_T(a)
						// 4: M^2/M_planck^2
						// 19: c(a) prefactor [needed to set it to 0 manually when assuming a LCDM background]

						// Initialise alpha_i (a) and their derivatives
						for(int i=0; i<5; i++){
							alphaofa[i] = alphai_eft(a,omega0,extpars[i],i+1); // alpha(a)
							dalphaofa[i] = dalphai_eft(a,omega0,extpars[i],i+1); // alpha'(a)
							ddalphaofa[i] = ddalphai_eft(a,omega0,extpars[i],i+1); // alpha''(a)
						}

						// We take expressions from Chapter B of GtoPT.nb
						// Hubble
						hub = HAg(a,omega0,extpars,model)*h0;
						// Hubble scale factor derivative
						hubd = HA1g(a,omega0,extpars,model)*pow2(h0)/(a*hub);
						// Hubble 2nd scale factor derivative
						hubdd = HA3g(a,omega0,extpars,model)*h0;
						// Background Ricci scalar scale factor derivative
						R0d = riccibackgroundp(a,omega0,extpars,model);
						// speed of gravitational waves
						ct2 = 1. + alphaofa[3];

						var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha

						// protect against divergences when var1=0
						if (fabs(var1) <tol) {
							var1 = tol;
						}

						var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[1]/2.*ct2 + HA2g(a,omega0,extpars,model)
																					 + alphaofa[2] - alphaofa[3])
																					 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)/alphaofa[4]); //cs^2

						// QSA check for late times
						qsa_test = k0/(HAg(a,omega0,extpars,model)*h0*a) -var2;
						if (qsa_test<0 && a>aqsa && k0>kqsa) {
							warning("BeyondLCDM: QSA violated at a: %e and k: %e \n", a,k0);
						}

						// Note we have factorised M^2/G_N out of the following expressions
						// It then appears as 1/alphaofa[4] in the final result

						// c(a)  background function
						ca = 	extpars[19]*(- 3.*pow2(h0)*omega0/(2.*pow3(a)*alphaofa[4])
									- 1./2.*hub * (a*(ct2*(2.+alphaofa[2]) + a*dalphaofa[3])*hubd
					  			+ hub*(ct2*((alphaofa[2]-1.)*alphaofa[2] + a*dalphaofa[2]) + 2.*a*alphaofa[2]*dalphaofa[3]
									+ a*a*ddalphaofa[3])));

						// A terms
						myA[0] = 2.;
						myA[1] = alphaofa[1]*hub;
						myA[2] = 0.;
						// B terms
						myB[0] = - 1./ct2;
						myB[1] = 1.;
						myB[2] = (-	alphaofa[2] + alphaofa[3])*hub / ct2;
						// C terms
						myC[0] = - ct2 * myB[2];
						myC[1] = myA[1]/2.;

						myC[2] = ca + hub/2. * (-2.*alphaofa[3]*hub + pow2(alphaofa[2])*ct2*hub
						  			 + a*hub*dalphaofa[1] + a*hub*dalphaofa[2] + a*hub*alphaofa[3]*dalphaofa[2]
										 + 2.*a*alphaofa[3]*hubd + pow2(a)*dalphaofa[3]*hubd
										 + alphaofa[1]*((1.+alphaofa[2])*hub + a*hubd)
										 + alphaofa[2]*(hub*(1.-alphaofa[3] + 2.*a*dalphaofa[3]) + a*ct2*hubd)
									   + pow2(a)*hub*ddalphaofa[3]);


						myC[3] = -a*hub/4. * (12.*ca*hubd + hub*(6.*pow2(alphaofa[2])*ct2*hub*hubd
										+ 6.*alphaofa[1]*(2.*a*pow2(hubd) + hub*((4.+alphaofa[2])*hubd + a*hubdd))
									  + alphaofa[2]*(ct2*(12.*a*pow2(hubd) - R0d) + 6.*hub*(2.*(2.*ct2 + a*dalphaofa[3])*hubd + a*ct2*hubdd))
										+ a*(12.*(alphaofa[3] + a*dalphaofa[3])*pow2(hubd) - dalphaofa[3]*R0d
										+ 6.*hub*(hubd*(dalphaofa[1] + ct2*dalphaofa[2] + 5.*dalphaofa[3] + a*ddalphaofa[3]) + a*dalphaofa[3]*hubdd))));


						// f terms
						myf[0] = myB[1]*myC[2] - myC[0]*myB[2];
						myf[1] = myB[1]*myC[3];
						myf[2] = myA[0]*(myB[2]*myC[1]-myB[0]*myC[2]) + myA[1]*(myB[0]*myC[0]-myB[1]*myC[1]) + myA[2]*(myB[1]*myC[2]-myB[2]*myC[0]);
						myf[3] = myC[3]*(myA[2]*myB[1] - myA[0]*myB[0]);

						var3 = pow2(k0/a);
						// avoid singularity by shifting denominator by small amount
						tol = 1e-20;
						if (fabs(myf[2]*var3 + myf[3])<tol) {
							myf[3]*=0.001;
						}
						return 2./alphaofa[4]  *  (myf[0]*var3 + myf[1] ) / (myf[2]*var3 + myf[3]);


			case 11:
					/* EFTofDE unscreened and full k-dependence in linear modification */
						// 0: alpha_K(a)
						// 1: alpha_B(a)
						// 2: alpha_M(a)
						// 3: alpha_T(a)
						// 4: M^2/M_planck^2
						// 19: c(a) prefactor [needed to set it to 0 manually when assuming a LCDM background]

						// Initialise alpha_i (a) and their derivatives
						for(int i=0; i<5; i++){
							alphaofa[i] = alphai_eft(a,omega0,extpars[i],i+1);
							dalphaofa[i] = dalphai_eft(a,omega0,extpars[i],i+1);
							ddalphaofa[i] = ddalphai_eft(a,omega0,extpars[i],i+1);
						}

						// We take expressions from Chapter B of GtoPT.nb
						// Hubble
						hub = HAg(a,omega0,extpars,model)*h0;
						// Hubble scale factor derivative
						hubd = HA1g(a,omega0,extpars,model)*pow2(h0)/(a*hub);
						// Hubble 2nd scale factor derivative
						hubdd = HA3g(a,omega0,extpars,model)*h0;
						// Background Ricci scalar scale factor derivative
						R0d = riccibackgroundp(a,omega0,extpars,model);
						// speed of gravitational waves
						ct2 = 1. + alphaofa[3];

						var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha

						// protect against divergences when var1=0
						if (fabs(var1) <tol) {
							var1 = tol;
						}

						var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[1]/2.*ct2 + HA2g(a,omega0,extpars,model)
																					 + alphaofa[2] - alphaofa[3])
																					 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)/alphaofa[4]); //cs^2

						 // protect against divergences when var2=0
 						if (fabs(var2) <tol) {
 							var2 = tol;
 						}

						// QSA check for late times
						qsa_test = k0/(HAg(a,omega0,extpars,model)*h0*a) -var2;
						if (qsa_test<0 && a>aqsa && k0>kqsa) {
							warning("BeyondLCDM: QSA violated at a: %e and k: %e \n", a,k0);
						}


						// Note we have factorised M^2/G_N out of the following expressions
						// It then appears as 1/alphaofa[4] in the final result

						// c(a) background function
						ca = extpars[19]*(- 3.*pow2(h0)*omega0/(2.*pow3(a)*alphaofa[4])
									- 1./2.*hub * (a*( ct2*(2.+alphaofa[2]) + a*dalphaofa[3])*hubd
					  			+ hub*(ct2*((alphaofa[2]-1.)*alphaofa[2] + a*dalphaofa[2]) + 2.*a*alphaofa[2]*dalphaofa[3]
									+ a*a*ddalphaofa[3])));


						// A terms
						myA[0] = 2.;
						myA[1] = alphaofa[1]*hub;
						myA[2] = 0.;
						// B terms
						myB[0] = - 1./ct2;
						myB[1] = 1.;
						myB[2] = (-	alphaofa[2] + alphaofa[3])*hub / ct2;
						// C terms
						myC[0] = -ct2 * myB[2];
						myC[1] = myA[1]/2.;

						myC[2] = ca + hub/2. * (-2.*alphaofa[3]*hub + pow2(alphaofa[2])*ct2*hub
						  			 + a*hub*dalphaofa[1] + a*hub*dalphaofa[2] + a*hub*alphaofa[3]*dalphaofa[2]
										 + 2.*a*alphaofa[3]*hubd + pow2(a)*dalphaofa[3]*hubd
										 + alphaofa[1]*((1.+alphaofa[2])*hub + a*hubd)
										 + alphaofa[2]*(hub*(1.-alphaofa[3] + 2.*a*dalphaofa[3]) + a*ct2*hubd)
									   + pow2(a)*hub*ddalphaofa[3]);


						myC[3] = -a*hub/4. *  (12.*ca*hubd + hub*(6.*pow2(alphaofa[2])*ct2*hub*hubd
										+ 6.*alphaofa[1]*(2.*a*pow2(hubd) + hub*((4.+alphaofa[2])*hubd + a*hubdd))
									  + alphaofa[2]*(ct2*(12.*a*pow2(hubd) - R0d) + 6.*hub*(2.*(2.*ct2 + a*dalphaofa[3])*hubd + a*ct2*hubdd))
										+ a*(12.*(alphaofa[3] + a*dalphaofa[3])*pow2(hubd) - dalphaofa[3]*R0d
										+ 6.*hub*(hubd*(dalphaofa[1] + ct2*dalphaofa[2] + 5.*dalphaofa[3] + a*ddalphaofa[3]) + a*dalphaofa[3]*hubdd))));


						// f terms
						myf[0] = myB[1]*myC[2] - myC[0]*myB[2];
						myf[1] = myB[1]*myC[3];
						myf[2] = myA[0]*(myB[2]*myC[1]-myB[0]*myC[2]) + myA[1]*(myB[0]*myC[0]-myB[1]*myC[1]) + myA[2]*(myB[1]*myC[2]-myB[2]*myC[0]);
						myf[3] = myC[3]*(myA[2]*myB[1] - myA[0]*myB[0]);

						var3 = pow2(k0/a);
						tol = 1e-20;
						if (fabs(myf[2]*var3 + myf[3])<tol) {
							myf[3]*=0.001;
						}

						return 2./alphaofa[4]  *  (myf[0]*var3 + myf[1]) / (myf[2]*var3 + myf[3]);


				case 12:
				/* EFTofDE Phenomenological non-linear model  */

					// 0: alpha_K(a)
					// 1: alpha_B(a)
					// 2: alpha_M(a)
					// 3: alpha_T(a)
					// 4: M^2/M_planck^2
					// 19: c(a) prefactor [needed to set it to 0 manually when assuming a LCDM background]
					//
					// // Initialise alpha_i (a) and their derivatives
					for(int i=0; i<5; i++){
						alphaofa[i] = alphai_eft(a,omega0,extpars[i],i+1);
						dalphaofa[i] = dalphai_eft(a,omega0,extpars[i],i+1);
						ddalphaofa[i] = ddalphai_eft(a,omega0,extpars[i],i+1);
					}

					// We take expressions from Chapter B of GtoPT.nb
					// Hubble
					hub = HAg(a,omega0,extpars,model)*h0;
					// Hubble scale factor derivative
					hubd = HA1g(a,omega0,extpars,model)*pow2(h0)/(a*hub);
					// Hubble 2nd scale factor derivative
					hubdd = HA3g(a,omega0,extpars,model)*h0;
					// Background Ricci scalar scale factor derivative
					R0d = riccibackgroundp(a,omega0,extpars,model);
					// speed of gravitational waves
					ct2 = 1. + alphaofa[3];

					var1 = alphaofa[0] + 3./2.*pow2(alphaofa[1]); // alpha

					// protect against divergences when var1=0
					if (fabs(var1) <tol) {
						var1 = tol;
					}

					var2 = 2./var1*( (1.-	alphaofa[1]/2.) * (alphaofa[1]/2.*ct2 + HA2g(a,omega0,extpars,model)
																				 + alphaofa[2] - alphaofa[3])
																				 + a*dalphaofa[1]/2.- HA2g2(a,omega0,extpars,model)/alphaofa[4]); //cs^2

				 // protect against divergences when var2=0
					if (fabs(var2) <tol) {
						var2 = tol;
					}

					// QSA check for late times
					qsa_test = k0/(HAg(a,omega0,extpars,model)*h0*a) -var2;
					if (qsa_test<0 && a>aqsa && k0>kqsa) {
						warning("BeyondLCDM: QSA violated at a: %e and k: %e \n", a,k0);
					}


					// Note we have factorised M^2/G_N out of the following expressions
					// It then appears as 1/alphaofa[4] in the final result

					// c(a) background function
					ca = extpars[19]*(- 3.*pow2(h0)*omega0/(2.*pow3(a)*alphaofa[4])
								- 1./2.*hub * (a*( ct2*(2.+alphaofa[2]) + a*dalphaofa[3])*hubd
								+ hub*(ct2*((alphaofa[2]-1.)*alphaofa[2] + a*dalphaofa[2]) + 2.*a*alphaofa[2]*dalphaofa[3]
								+ a*a*ddalphaofa[3])));


					// A terms
					myA[0] = 2.;
					myA[1] = alphaofa[1]*hub;
					myA[2] = 0.;
					// B terms
					myB[0] = - 1./ct2;
					myB[1] = 1.;
					myB[2] = (-	alphaofa[2] + alphaofa[3])*hub / ct2;
					// C terms
					myC[0] = -ct2 * myB[2];
					myC[1] = myA[1]/2.;

					myC[2] = ca + hub/2. * (-2.*alphaofa[3]*hub + pow2(alphaofa[2])*ct2*hub
									 + a*hub*dalphaofa[1] + a*hub*dalphaofa[2] + a*hub*alphaofa[3]*dalphaofa[2]
									 + 2.*a*alphaofa[3]*hubd + pow2(a)*dalphaofa[3]*hubd
									 + alphaofa[1]*((1.+alphaofa[2])*hub + a*hubd)
									 + alphaofa[2]*(hub*(1.-alphaofa[3] + 2.*a*dalphaofa[3]) + a*ct2*hubd)
									 + pow2(a)*hub*ddalphaofa[3]);


					myC[3] = -a*hub/4. *  (12.*ca*hubd + hub*(6.*pow2(alphaofa[2])*ct2*hub*hubd
									+ 6.*alphaofa[1]*(2.*a*pow2(hubd) + hub*((4.+alphaofa[2])*hubd + a*hubdd))
									+ alphaofa[2]*(ct2*(12.*a*pow2(hubd) - R0d) + 6.*hub*(2.*(2.*ct2 + a*dalphaofa[3])*hubd + a*ct2*hubdd))
									+ a*(12.*(alphaofa[3] + a*dalphaofa[3])*pow2(hubd) - dalphaofa[3]*R0d
									+ 6.*hub*(hubd*(dalphaofa[1] + ct2*dalphaofa[2] + 5.*dalphaofa[3] + a*ddalphaofa[3]) + a*dalphaofa[3]*hubdd))));


					// f terms
					myf[0] = myB[1]*myC[2] - myC[0]*myB[2];
					myf[1] = myB[1]*myC[3];
					myf[2] = myA[0]*(myB[2]*myC[1]-myB[0]*myC[2]) + myA[1]*(myB[0]*myC[0]-myB[1]*myC[1]) + myA[2]*(myB[1]*myC[2]-myB[2]*myC[0]);
					myf[3] = myC[3]*(myA[2]*myB[1] - myA[0]*myB[0]);

					var3 = pow2(k0/a);
					tol = 1e-20;
					if (fabs(myf[2]*var3 + myf[3])<tol) {
						myf[3]*=0.001;
					}

					return 2./alphaofa[4]  *  (myf[0]*var3 + myf[1]) / (myf[2]*var3 + myf[3]);


				case 13:
					/* Minimal model independent : CPL + Gamma + Erf */
					// see eq.47 of PhysRevD.98.044051
					hub = HAg(a,omega0,extpars,model); // H/H0
					hubd = HA1g(a,omega0,extpars,model)/pow2(hub); // a H' / H
					var1 = omega0/pow3(a) / pow2(hub); // Omega_m
					var2 = pow(var1,extpars[2]); // Omega_m ^gamma
					var3 = -3. - 2. * hubd; // a Omega_m ' / Omega_m

					return 2./3. * var2/var1 * (var2 + 2. + hubd + extpars[2]*var3);


				case 14:
				/* CG */
					// Hubble and its derivative
					hub = HAg(a,omega0,extpars,model);
					////HA2g = -dH/dt/H^2 = -a dH/da / H
					hubd = HA2g(a,omega0,extpars,model);

					//alphaB
					alphaofa[0] = -(1.-omega0)*pow(1./hub,2.+extpars[0])*extpars[1]*extpars[0];
					// d/dt alphaB
					dalphaofa[0] = -alphaofa[0]*(2.+extpars[0])*h0*HA1g(a,omega0,extpars,model)/hub;

					// alphacs2
					var1 = 2.*hubd-3.*omega0/(pow3(a)*pow2(hub));
					var2 = -2.*alphaofa[0] + 2.*hubd*alphaofa[0] -2.*pow2(alphaofa[0])-2.*dalphaofa[0]/(h0*hub);
					var3 = var1 + var2;

				 return (1.+2.*pow2(alphaofa[0])/var3);

			 case 15:
	 		/* QCDM */
	 			 return  1. ;

			case 16:
	 			/* K-mouflage */
	 			// See Eq.61 of 1403.5424
	 			// normalise conformal hubble =  a H/H0
	 			hub = HAg(a,omega0,extpars,model)*a;
	 			// Canonical kinetic term
	 			var1 = pow2(hub*kmouflage_phid(a)/extpars[1]/a) / 2. ;
	 			// kinetic term derivative
	 			var2 = 1. + extpars[0] * extpars[2] * pow(var1,extpars[0]-1.);

	 			return  conf_fac(kmouflage_phi(a),extpars[3]) * (1. + 2. * pow2(extpars[3])/var2) ;

	 		case 17:
	 		/* K-mouflage with nPPF */
	 		// See Eq.61 of 1403.5424
	 		// normalise conformal hubble =  a H/H0
	 		hub = HAg(a,omega0,extpars,model)*a;
	 		// Canonical kinetic term
	 		var1 = pow2(hub*kmouflage_phid(a)/extpars[1]/a) / 2. ;
	 		// kinetic term derivative
	 		var2 = 1. + extpars[0] * extpars[2] * pow(var1,extpars[0]-1.);

	 		return  conf_fac(kmouflage_phi(a),extpars[3]) * (1. + 2. * pow2(extpars[3])/var2);

		case 18:
		/*  Model independent Poisson equation parametrisation as in Eq. 3.6 of https://arxiv.org/pdf/2411.12026 with CPL background  */
	
		// time dependency
		hub = HAg(a,omega0,extpars,model); // H/H0 
		var1 = -3.*(1.+extpars[0]+extpars[1]);
		var2 = pow(a,var1)*exp(3.*(-1.+a)*extpars[1]) / pow2(hub); // Omega_DE / Omega_DE(z=0)  

		// scale dependency
		var3  = (1. + extpars[3] * pow2(extpars[4]* hub / k0 ) ) /  (1. +  pow2(extpars[4]* hub / k0 ) );  

			return 1. + extpars[2] * var2 * var3; 

		default:
				warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
}
}


// Gamma 2
double gamma2(double a, double omega0, double k0, double k1, double k2, double u1, double extpars[], int model){
	double h0 = 1./2997.92458;
	double var1,var2,var3,var4,var5,var6,hub,hubd;
	double alphaofa[5],dalphaofa[5];
	switch(model) {
		case 1:
		/* GR */
						return  0.;
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
						return  0.;
		case 5:
		/* CPL */
						return  0.;
		case 6:
		/* Dark Scattering with CPL  */
						return  0.;
		case 7:
		/* EFTofDE nPPF*/
						return  0.;
		case 8:
		/* EFTofDE - unscreened  approximation  */
						return  0.;
		case 9:
		/* EFTofDE - superscreened approximation  */
						return  0.;
		case 10:
		/* EFTofDE with nPPF and full k-dependence in linear modification */
					return 0.;
		case 11:
		/* EFTofDE - unscreened approximation with full k-dependence in linear modification */
					return  0.;

		case 12:
		/* EFTofDE: Phenomenological non-linear model */
					return  0.;

		case 13:
		/* Minimal model independent : CPL + Gamma + Erf */
					return 0.;

		case 14:
			/* CG */
			// Hubble and its derivative
			hub = HAg(a,omega0,extpars,model);
			hubd = HA2g(a,omega0,extpars,model);

			//alphaB
			alphaofa[0] = -(1.-omega0)*pow(1./hub,2.+extpars[0])*extpars[1]*extpars[0];

			// d/dt alphaB
			dalphaofa[0] = -alphaofa[0]*(2.+extpars[0])*h0*HA1g(a,omega0,extpars,model)/hub;
			// alphacs2
			var1 = 2.*hubd-3.*omega0/(pow3(a)*pow2(hub));
			var2 = -2.*alphaofa[0] + 2.*hubd*alphaofa[0] -2.*pow2(alphaofa[0])-2.*dalphaofa[0]/(h0*hub);
			var3 = var1 + var2;

			return  -18.*pow4(alphaofa[0])/(pow4(hub)*pow3(var3))*pow2(omega0/pow3(a))*ker1(u1);

		 case 15:
			/* QCDM */
			return  0. ;

		case 16:
			/* K-mouflage */
			return 0;

		case 17:
			/* K-mouflage with nPPF */
			return 0;

		case 18:
		/*  Model independent Poisson equation parametrisation as in Eq. 3.6 of https://arxiv.org/pdf/2411.12026 with CPL background  */
			return 0.; 

		default:
		warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
}
}

// Partially symmetric gamma3 (in last 2 arguments)
double gamma3(double a, double omega0, double k0, double k1, double k2, double k3, double u1,double u2, double u3, double extpars[], int model){
	double h0 = 1./2997.92458;
	double var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,k23,hub,hubd;
	double alphaofa[5],dalphaofa[5];

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

						  return  1./3.*(var8*pow2(1./HA(a,omega0))*pow3(omega0/var1)/(36.
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
		/* EFTofDE nPPF  */
						return  0. ;
		case 8:
		/* EFTofDE - unscreened  approximation  */
						return  0. ;
		case 9:
		/* EFTofDE - superscreened approximation  */
						return  0. ;
		case 10:
		/* EFTofDE: full k dependency assuming LCDM background with Phenomenological approximation  */
					return 0.;

		case 11:
		/* EFTofDE: full k dependency assuming LCDM background with unscreened approximation  */
					return 0.;

		case 12:
		/* EFTofDE: full k dependency assuming LCDM background with Phenomenological approximation  */
					return  0.;

		case 13:
		/* Minimal model independent : CPL + Gamma + Erf */
					return 0.;

		case 14:
			/* CG */
			// Hubble and its derivative
			hub = HAg(a,omega0,extpars,model);
			hubd = HA2g(a,omega0,extpars,model);

			//alphaB
			alphaofa[0] = -(1.-omega0)*pow(1./hub,2.+extpars[0])*extpars[1]*extpars[0];
			// d/dt alphaB
			dalphaofa[0] = -alphaofa[0]*(2.+extpars[0])*h0*HA1g(a,omega0,extpars,model)/hub;
			// alphacs2
			var1 = 2.*hubd-3.*omega0/(pow3(a)*pow2(hub));
			var2 = -2.*alphaofa[0] + 2.*hubd*alphaofa[0] -2.*pow2(alphaofa[0])-2.*dalphaofa[0]/(h0*hub);
			var3 = var1 + var2;

			return  72.*pow(alphaofa[0],6)/(pow(hub,6)*pow(var3,5))*pow3(omega0/pow3(a))*ker1(u1)*ker1((k2*u2+k3*u3)/sqrt(k2*k2+2.*k2*k3*u1+k3*k3));

	 case 15:
		/* QCDM */
		return  0. ;

	case 16:
		/* K-mouflage */
		// See Eq.78 of 1403.5424

		// normalise conformal hubble =  a H/H0
		hub = HAg(a,omega0,extpars,model)*a;

		// Canonical kinetic term
		var1 = pow2(hub*kmouflage_phid(a)/extpars[1]/a) / 2. ;
		// kinetic term derivative
		var2 = 1. + extpars[0] * extpars[2] * pow(var1,extpars[0]-1.);
		// kinetic term 2nd derivative
		var3 = extpars[0] * (extpars[0]-1.) * extpars[2] * pow(var1,extpars[0]-2.);

		// Omega_m(a)
		var4 = omega0/pow2(hub)/a;
		// conformal factor
		var5 = conf_fac(kmouflage_phi(a),extpars[3]);

		// Symmetrised scale dependant piece x 3 k^2
		var6 = (u2 + 2.*u3*u1)/k1/k2 + (u3 + 2.*u1*u2)/k1/k3 + (u1 + 2.*u2*u3)/k2/k3;


		return -9./2. * var3 * pow3(var5 * var4) * pow(extpars[3]/var2,4.) * pow2(hub/extpars[1]) * pow2(hub*h0/a) *var6 ;


	case 17:
		/* K-mouflage with nPPF */
		// See Eq.78 of 1403.5424

		// normalise conformal hubble =  a H/H0
		hub = HAg(a,omega0,extpars,model)*a;
		// Canonical kinetic term
		var1 = pow2(hub*kmouflage_phid(a)/extpars[1]/a) / 2. ;
		// kinetic term derivative
		var2 = 1. + extpars[0] * extpars[2] * pow(var1,extpars[0]-1.);
		// kinetic term 2nd derivative
		var3 = extpars[0] * (extpars[0]-1.) * extpars[2] * pow(var1,extpars[0]-2.);

		// Omega_m(a)
		var4 = omega0/pow2(hub)/a;
		// conformal factor
		var5 = conf_fac(kmouflage_phi(a),extpars[3]);

		// Symmetrised scale dependant piece x 3 k^2
		var6 = (u2 + 2.*u3*u1)/k1/k2 + (u3 + 2.*u1*u2)/k1/k3 + (u1 + 2.*u2*u3)/k2/k3;


		return -9./2. * var3 * pow3(var5 * var4) * pow(extpars[3]/var2,4.) * pow2(hub/extpars[1]) * pow2(hub*h0/a) *var6 ;

	case 18:
		/*  Model independent Poisson equation parametrisation as in Eq. 3.6 of https://arxiv.org/pdf/2411.12026 with CPL background  */
		return 0.; 
			

		default:
		warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				return 0;
}
}

/* Modified non-linear Poisson equation function (F) for spherical collapse (see for example A.1 of appendix of 1812.05594) */

double mymgF(double a, double yh, double yenv, double Rth, double omega0, double extpars[], double delta_initial, int model){
	double h0 = 1./2997.92;
	double dod, dod2, dRRth, fr0, var1, var2, var3, term1, term2, term3,hub,hubd;
	double betadgp,xm3,xterm,delta,Mvir;
	double Aphi;
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

						dod = fr0 * yh * a * var1 * (term1 - term2) / Rth/ Rth / omega0 ;
						dod2 = pow2(dod);

						dRRth = dod - dod2 + dod2*dod/3.;

						return gsl_min(dRRth,1./3.);

	   case 3:
	 	 /*DGP normal branch*/
		      betadgp = beta(a,omega0,extpars[0]);

		      delta = (1.+delta_initial)/pow3(yh) - 1.;

		      xterm = 9.*pow2(betadgp)*extpars[0]/(2.*omega0/pow3(a)*(delta));

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
			/* EFTofDE - nPPF */
			// extpars[i]:
			// 0: alpha_K(a)
			// 1: alpha_B(a)
			// 2: alpha_M(a)
			// 3: alpha_T(a)
			// 4: M^2/M_planck^2
			// 5-12 : p1-p8

					var1 = extpars[5]/(extpars[5]-1.)*extpars[7]; // a
					Mvir = pow3(Rth/0.1)*5.*omega0/Gnewton; // virial mass x Gnewton - see definition in scol_init in HALO.cpp

					term1 = extpars[8]*pow(a,extpars[9]) * pow(2.*Gnewton*h0*Mvir/pow2(sofl),extpars[10])*pow(yenv/yh,extpars[11]); // y0

					xm3 = pow(term1/yh,var1); // (y0/yh)^a

					// Linear modification
					var2 = pow(10.,extpars[12])/ (yh * pow2(a) * Rth); // Fourier transform of Rth with some calibration
					term3 = mu(a,var2,omega0,extpars,model)-1.; // Linear G_effective


					return extpars[5]*extpars[6]*(pow(1.+ xm3, 1./extpars[5])-1.) / xm3 * term3 ;

			case 8:
			/* EFTofDE - unscreened approximation  */
					return mu(a,0.001,omega0,extpars,model)-1.;

			case 9:
			/* EFTofDE - superscreened approximation  */
					return 0.;

			case 10:
			/* EFTofDE with nPPF and full k-dependence in linear modification */
			// extpars[i]:
			// 0: alpha_K(a)
			// 1: alpha_B(a)
			// 2: alpha_M(a)
			// 3: alpha_T(a)
			// 4: M^2/M_planck^2
			// 5-12 : p1-p8

					var1 = extpars[5]/(extpars[5]-1.)*extpars[7]; // a
					Mvir = pow3(Rth/0.1)*5.*omega0/Gnewton; // virial mass x Gnewton - see definition in scol_init in HALO.cpp

					term1 = extpars[8]*pow(a,extpars[9]) * pow(2.*Gnewton*h0*Mvir/pow2(sofl),extpars[10])*pow(yenv/yh,extpars[11]); // y0

					xm3 = pow(term1/yh,var1); // (y0/yh)^a

					// Linear modification
					var2 = pow(10.,extpars[12])/ (yh * pow2(a) * Rth); // Fourier transform of Rth with some calibration
					term3 = mu(a,var2,omega0,extpars,model)-1.; // Linear G_effective

			    return extpars[5]*extpars[6]*(pow(1.+ xm3, 1./extpars[5])-1.) / xm3 * term3;

			case 11:
			/* EFTofDE with full k-dependence and unscreened approximation  */
					return mu(a,1.,omega0,extpars,model)-1.;


			case 12:
			/* EFTofDE: Phenomenological non-linear model */
					var1 = extpars[5] - extpars[6]*log(Rth); // Screening scale and mass dependence : Erf [ yh (R_scr/R_th)^p2 ] with R_scr = 10^(p1/p2)
					term1 = yh * a * pow(10,var1); // Argument for mass dependant error function

					if (yenv<=1.) {
					term2 = pow(yenv * a, extpars[7]);  // yenv dependence - power law
					}
					else{
						term2 = 1.;
					}

					var2 = pow(10.,extpars[8])/ (yh * pow2(a) * Rth); // Fourier transform of Rth with some calibration
					term3 = mu(a,var2,omega0,extpars,model)-1.; // Linear G_effective

					return term3 * erf(term1*term2);

			case 13:
			/* Minimal model independent : CPL + Gamma + Erf */
					var1 = extpars[3] - extpars[4]*log(Rth); // Screening scale and mass dependence : Erf [ yh (R_scr/R_th)^p2 ] with R_scr = 10^(p1/p2)
					term1 = yh * a * pow(10,var1); // Argument for mass dependant error function

					if (yenv<=1.) {
					term2 = pow(yenv * a, extpars[5]);  // yenv dependence - power law
					}
					else{
						term2 = 1.;
					}

					var2 = pow(10.,extpars[6])/ (yh * pow2(a) * Rth); // Fourier transform of Rth with some calibration
					term3 = mu(a,var2,omega0,extpars,model)-1.; // Linear G_effective

					return term3 * erf(term1*term2);

			case 14:
			/* CG */

				// Hubble and its derivative
				hub = HAg(a,omega0,extpars,model);
				// HA2g = -dH/dt/H^2 = -a dH/da / H
				hubd = HA2g(a,omega0,extpars,model);

				//alphaB
				alphaofa[0] = -(1.-omega0)*pow(1/hub,2.+extpars[0])*extpars[1]*extpars[0];
				// d/dt alphaB
				dalphaofa[0] = -alphaofa[0]*(2.+extpars[0])*h0*HA1g(a,omega0,extpars,model)/hub;

				// alphacs2
				var1 = 2.*hubd-3.*omega0/(pow3(a)*pow2(hub));
				var2 = -2.*alphaofa[0] + 2.*hubd*alphaofa[0] -2.*pow2(alphaofa[0])-2.*dalphaofa[0]/(h0*hub);
				var3 = var1 + var2;


				delta = (1.+delta_initial)/pow3(yh) - 1.;
				xterm = -2.*alphaofa[0]/(hub*var3);
				xm3 = 4.*omega0 * pow2(xterm)/pow3(a) * delta ;

				return	2.*(2.*pow2(alphaofa[0])/var3)*(sqrt(1. + xm3)-1.)/xm3;


		 case 15:
			/* QCDM */
			 return  0. ;


			 case 16:
				/* K-mouflage */
				// See Eq. 115
				 return  mu(a,1.,omega0,extpars,model)-1.;

			case 17:
			/* K-mouflage with nPPF  */
			// extpars[i]:
			// 0: n
			// 1: lambda
			// 2: K0
			// 3: betaK
			// 5-12 : p1-p8

					var1 = extpars[5]/(extpars[5]-1.)*extpars[7]; // a
					Mvir = pow3(Rth/0.1)*5.*omega0/Gnewton; // virial mass x Gnewton - see definition in scol_init in HALO.cpp

					term1 = extpars[8]*pow(a,extpars[9]) * pow(2.*Gnewton*h0*Mvir/pow2(sofl),extpars[10])*pow(yenv/yh,extpars[11]); // y0

					xm3 = pow(term1/yh,var1); // (y0/yh)^a

					// Linear modification
					var2 = pow(10.,extpars[12])/ (yh * pow2(a) * Rth); // Fourier transform of Rth with some calibration
					Aphi = conf_fac(kmouflage_phi(a),extpars[3]) ;

					term3 = mu(a,var2,omega0,extpars,model)-Aphi; // Linear G_effective/G_N corrected with conformal factor included

					return extpars[5]*extpars[6]*(pow(1.+ xm3, 1./extpars[5])-1.) / xm3 * term3 - 1.;


			case 18:
			/*  Model independent Poisson equation parametrisation as in Eq. 3.6 of https://arxiv.org/pdf/2411.12026 with CPL background  */
				var1 = extpars[5] - extpars[6]*log(Rth); // Screening scale and mass dependence : Erf [ yh (R_scr/R_th)^p2 ] with R_scr = 10^(p1/p2)
				term1 = yh * a * pow(10,var1); // Argument for mass dependant error function

				if (yenv<=1.) {
				term2 = pow(yenv * a, extpars[7]);  // yenv dependence - power law
				}
				else{
					term2 = 1.;
				}

				var2 = pow(10.,extpars[8])/ (yh * pow2(a) * Rth); // Fourier transform of Rth with some calibration
				term3 = mu(a,var2,omega0,extpars,model)-1.; // Linear G_effective

				return term3 * erf(term1*term2);

		default:
					warning("BeyondLCDM: invalid model choice, model = %d \n", model);
				  return 0;
  }
}

// Dark energy contribution to virial theorem - see Eq.A9 of 1812.05594
double  WEFF(double a, double omega0, double extpars[], int model){
	double h2,hubd,var1,var2,var3;
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
		 /* EFTofDE: k->infinity limit with nPPF */
		 return 2.*(1.-omega0);

		 case 8:
		 /* EFTofDE:  k->infinity limit with  unscreened approximation */
		 return 2.*(1.-omega0);

		 case 9:
		 /* EFTofDE:  k->infinity limit with superscreened approximation */
		 return 2.*(1.-omega0);

		 case 10:
		 /* EFTofDE: full k dependency assuming LCDM background with nPPF approximation */
		 return 2.*(1.-omega0);

		 case 11:
		 /* EFTofDE: full k dependency assuming LCDM background with unscreened approximation  */
		 return 2.*(1.-omega0);

		 case 12:
 		/* EFTofDE: full k dependency assuming LCDM background with Phenomenological approximation  */
 		 return  2.*(1.-omega0);

		 case 13:
		 /* Minimal model independent : CPL + Gamma + Erf */
 		 h2 = pow2(HAg(a,omega0,extpars,model));
 		 return -(1.+3.*(extpars[0]+(1.-a)*extpars[1]))*(h2-omega0/pow3(a));

		 case 14:
		 /* CG */
		 // Hubble and its derivative
		 // HA2g = -dH/dt/H^2 = -a dH/da / H
		 h2 = pow2(HAg(a,omega0,extpars,model));
		 hubd = HA2g(a,omega0,extpars,model);

		return (2.+hubd*extpars[0])*(h2-omega0/pow3(a));


		 case 15:
		 /* QCDM */
		 // Hubble and its derivative
		 h2 = pow2(HAg(a,omega0,extpars,model));
		 hubd = HA2g(a,omega0,extpars,model);

		return (2.+hubd*extpars[0])*(h2-omega0/pow3(a));


		case 16:
	 /* K-mouflage */
		 h2 = pow2(HAg(a,omega0,extpars,model));
		 // Canonical kinetic term
		 var1 = h2*pow2(kmouflage_phid(a)/extpars[1]) / 2. ;
		 // non-canonical kinetic term
		 var2 = -1. + var1 + extpars[2]*pow(var1,extpars[0]);
		 // kinetic term derivative
		 var3 = 1. + extpars[0] * extpars[2] * pow(var1,extpars[0]-1.);

		 return -pow2(extpars[1])/3. * (2.*var2 + h2*pow2(kmouflage_phid(a)/extpars[1])*var3);


		case 17:
	 /* K-mouflage with nPPF   */
		 h2 = pow2(HAg(a,omega0,extpars,model));
		 // Canonical kinetic term
		 var1 = h2*pow2(kmouflage_phid(a)/extpars[1]) / 2. ;
		 // non-canonical kinetic term
		 var2 = -1. + var1 + extpars[2]*pow(var1,extpars[0]);
		 // kinetic term derivative
		 var3 = 1. + extpars[0] * extpars[2] * pow(var1,extpars[0]-1.);

		 return -pow2(extpars[1])/3. * (2.*var2 + h2*pow2(kmouflage_phid(a)/extpars[1])*var3);


		 case 18:
	 	/*  Model independent Poisson equation parametrisation as in Eq. 3.6 of https://arxiv.org/pdf/2411.12026 with CPL background  */		
		 h2 = pow2(HAg(a,omega0,extpars,model));
		return -(1.+3.*(extpars[0]+(1.-a)*extpars[1]))*(h2-omega0/pow3(a));
 

		 default:
				 warning("BeyondLCDM: invalid model choice, model = %d \n", model);
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

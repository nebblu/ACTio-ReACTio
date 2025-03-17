#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cmath>


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <cmath>

#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SPT.h>
#include <Copter/GrowthFunction.h>
#include <Copter/SpecialFunctions.h>
#include <Copter/NoWigglePS.h>
#include <Copter/BeyondLCDM.h>

using namespace std;

/* Example code to output the EFTofLSS power spectrum as given in Eq.41 of 2311.13529 (should be equivalent to PBJ output)*/
/* Note that the higher order bias terms have no scale dependent corrections from MG and the growth are all calculated at the external momentum for these terms */

int main(int argc, char* argv[]) {

  // Which gravity or dark energy model?
  // 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL, 6: Dark scattering with CPL
  // 7 -9 : EFTofDE k->infinity limit,  w PPF, unscreened, superscreened respectively
  // 10-12: EFTofDE with PPF, unscreened and Error function Phenomenological model respectively.
      int mymodel = 2;

    // chosen redshift
      double myz = 0.0;
    // chosen omega_matter (total)
      double omega0 = 0.3132;
     // MG parameter (for f(R) this is |fr0| and for nDGP this is Omega_rc)
      double mgpar = 1e-4;

      //output file name
      const char* output = "EFTofLSS_f4_P024_cdm_z0.dat";
//	const char* output = "PL_z0.dat";

      // number of bins between kmin and kmax
      int Nk =50;
      double kmin = 0.01;
      double kmax = 0.4;

      // choose RSD model :
      // 0: Kaiser
      // 1: TNS with q-bias (see 1507.01592 for example)
      // 2: TNS with Lagrangian bias (incomplete!)
      // 3: SPT 1-loop RSD with linear bias only
      // 4: SPT 1-loop RSD with higher order bias
      // 5: shot noise and EFTofLSS counter terms
      double rsd_model = 4;

      // Lagrangian bias params (See 1607.03149 for example) OR Q-bias params (see PRSD_mg in src/SPT.cpp).
      double bias[8];
    /*  bias[0] = 1.5; // b1
      bias[1] = 0.7; // b2
      bias[2] = 0.3; // bG2
      bias[3] = -0.7; // bG3
      bias[4] = 200.; // N
      bias[5] = 5.; // epsilon_0
      bias[6] = 200.; // epsilon_2
      bias[7] = 0.001.; // nbar
      */
      bias[0] = 1.; // b1
      bias[1] = 0.; // b2
      bias[2] = 0.; // bG2
      bias[3] = 0.; // bG3
      bias[4] = 0.; // N
      bias[5] = 0.; // epsilon_0
      bias[6] = 0.; // epsilon_2
      bias[7] = 0.001; // nbar

      // absolute error in RSD integrals
      double err = 1e-3;

    // cosmo file for P_L(k)
    const char* cstr = "transfers/planck";

    // Keep it z=0 to keep Copter's Growth @ 1
    real z =0;
    // Relative error in magnitude integrations
    real epsrel = 1e-3;

    // Initialise all classes
    Cosmology C(cstr);
    LinearPS P_l(C, z);
    SPT spt(C, P_l, epsrel);
    IOW iow;

    // base parameter values
    double pars[10];
    pars[0] = 1./(1.+myz); // scale factor
    pars[1] =  omega0; // Omega_{m,0}
    pars[2] =  0.; // Omega_{nu,0}

    // ODE solver initial step
    pars[7] = 1e-4;
    // ODE solver relative error
    pars[8] = 1e-3;
    // ODE solver absolute error
    pars[9] = 1e-3;

    // extended model parameters
    double extpars[maxpars]; // Currently maxed out at 20 extra params
    extpars[0] = mgpar;

    // RSD parameters
    // 0: sigma_v in linear theory
    // 1: c0
    // 2: c2
    // 3: c4
    double rsdpars[4];

    // initialise LCDM or wCDM lin growth for PS normalisation
    iow.initnorm(pars,extpars,mymodel);
    rsdpars[0] = spt.sigmav_init(pars, extpars, mymodel); // calculates linear theory dispersion
    rsdpars[1] = 5.;
    rsdpars[2] = 10.;
    rsdpars[3] = 20.;


    /* Output section */
    /* Open output file */
    FILE* fp = fopen(output, "w");

    real p1,p2,p3,p4,p5,p6,p7,p8,p9;


for(int i =0; i <Nk;  i ++) {

     real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));

    /* for more quantities check out the SPT.h file in the src directory */
      p1 = spt.PRSD_mg(rsd_model, 1, pars, extpars, rsdpars, bias, k, err, mymodel); // P0 SPT
      p2 = spt.PRSD_mg(rsd_model, 2, pars, extpars, rsdpars, bias, k, err, mymodel); // P2 SPT
      p3 = spt.PRSD_mg(rsd_model, 3, pars, extpars, rsdpars, bias, k, err, mymodel); // P4 SPT

      p4 = spt.PRSD_mg(5, 1, pars, extpars, rsdpars, bias, k, err, mymodel); // P0 ct+shot
      p5 = spt.PRSD_mg(5, 2, pars, extpars, rsdpars, bias, k, err, mymodel); // P2 ct+shot
      p6 = spt.PRSD_mg(5, 3, pars, extpars, rsdpars, bias, k, err, mymodel); // P4 ct+shot


     printf("%e %e %e %e %e %e %e  \n", k, p1,p2,p3,p4,p5,p6); // print to terminal
     fprintf(fp,"%e %e %e %e %e %e %e  \n", k, p1,p2,p3,p4,p5,p6); // print to file

//	p1 = P_l(k);
//	 printf("%e %e  \n", k, p1);
//	 fprintf(fp,"%e %e  \n", k, p1);

}


	/*close output file*/
    fclose(fp);
    return 0;
}

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

#include <omp.h>

#include <Copter/HALO.h>
#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SpecialFunctions.h>
#include <Copter/SPT.h>

using namespace std;

using std::ifstream;
using std::string;
using std::istringstream;


/* Example code to output the halo model powerspectrum for modified gravity */

int main(int argc, char* argv[]) {

  // Which gravity or dark energy model?
  // 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL
  int mymodel = 2;

  // Modified gravity active?
  bool modg = true;
  // Is the transfer being fed to ReACT of the target cosmology?
  // If false, the transfer should be LCDM at z=0.
  bool mgcamb = false;


   // cosmo file for P_L(k)
   const char* cstr = "transfers/Matteo_fr";
   // Keep it z=0 to keep Copter's Growth @ 1
   real z = 0;
   // Relative error in magnitude integrations
   real epsrel = 1e-3;
   // Initialise all classes
   Cosmology C(cstr);
   LinearPS P_l(C, z);
   HALO halo(C, P_l,P_l,P_l,P_l, epsrel);
   SPT spt(C, P_l, epsrel);
   IOW iow;

   // chosen redshift
      double myz = 0.;
   // chosen omega_matter (total)
     double omegam0 = 0.3072;
   // chosen omega_nu (total)
     double omeganu0 = 0.;

   // extended model parameters
   double extpars[maxpars]; // Currently maxed out at 20 extra params
   extpars[0] = 1e-5; // fr0 for f(R)

   // Base parameters
    double pars[7];
    pars[0] = 1./(1.+myz); // scale factor
    pars[1] = omegam0;  // chosen omega_matter (total)
    pars[2] = omeganu0;  // chosen omega_matter (total)
    pars[5] = 50.; // number of mass bins for spherical collapse solver
    pars[6] = 0.0; // omega_neutrinos

//initialise spherical collapse quantities and reaction quantities
halo.initialise(pars,extpars,mgcamb,modg,mymodel);

/* Output section */

//output file name
const char* output = "halo_ps.dat";

/* Open output file */
FILE* fp = fopen(output, "w");

real p1,p2,p3,p4,p5,p6;

int Nk =10;
double kmin = 0.001;
double kmax = 10.;

 for(int i =0; i < Nk;  i ++) {

  real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));

      p1 =  halo.one_halo(k, pars); // one halo term : real
      p2 =  halo.one_halop(k, pars); // one halo term : pseudo
      p3 =  halo.reaction_nu(k, pars, mgcamb); // halo model reaction

     printf("%d %e %e %e %e \n", i, k, p1,p2,p3); // print to terminal
     fprintf(fp,"%e %e %e %e \n", k, p1,p2, p3); // print to file

}

	/*close output file*/
    fclose(fp);
    return 0;
}

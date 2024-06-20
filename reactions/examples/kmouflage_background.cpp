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

#include <Copter/HALO.h>
#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SpecialFunctions.h>
#include <Copter/SPT.h>
#include <Copter/BeyondLCDM.h>

//using namespace std;
using std::ifstream;
using std::string;
using std::istringstream;

vector<vector<double> > mytrans;
vector<vector<double> > mytransl;
vector<vector<double> > mypk;


/* Example code to output the reaction and halo spectra for minimal parametrisation */
int main(int argc, char* argv[]) {


// Which gravity or dark energy model?
// 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL, 6: Dark scattering with CPL
// 7 -9 : EFTofDE k->infinity limit,  w PPF, unscreened, superscreened respectively
// 10-12: EFTofDE with PPF, unscreened and Error function Phenomenological model respectively.
// 13: Minimal model independent parametrisation : w0,wa for background, gamma for linear, Phenomenological for nonlinear
// 14: Generalised cubic galileon
// 15: QCDM
// 16: K-mouflage

int mymodel = 16;

// K-mouflage parameters
double n = 2.;
double beta0 = 0.1; // glam
double K0 = 1.;
double lambda = 1.4518; //K0=1; b=0.1

// Omega_m today
double Omega_m = 0.3089;

// store beyond-LCDM parameters
double extpars[maxpars];
    // Linear and background
    extpars[0] = n;
    extpars[1] = lambda;
    extpars[2] = K0;
    extpars[3] = beta0;

double pars[8];
pars[0] = 1.;// scale factor
pars[1] = Omega_m;
pars[2] = 0.;

IOW iow;

// initialise K-mouflage scalar field and derivative
init_kmouflage_lna(pars,extpars,1000);
printf("%s \n ", "Initialised scalar field");

// initialise Hubble
iow.hubble_init(Omega_m,extpars,1000,mymodel);
printf("%s \n ", "Initialised Hubble");

/* Output section */

real p1,p2,p3,p4,p5,p6;

const char* output = "kmouflage_background_z0.dat";


FILE* fp = fopen(output, "w");

double Na = 100;
double amin = 0.001;
double amax = 1.;

 for(int i =0; i <Na;  i ++) {

      real a = amin * exp(i*log(amax/amin)/(Na-1.));

      p1 = HAg(a, Omega_m, extpars, mymodel); // Kmouflage Hubble/H0
      p2 = HA(a , Omega_m); // LCDM Hubble/H0
      p3 = HA1g(a, Omega_m, extpars, mymodel)/p1/a; // dH/da / H0^2
      p4 = kmouflage_sf(a); // Planck mass normalised scalar field
      p5 = kmouflage_sf_der(a); // its derivative wrt ln a

      printf("%e %e  %e %e %e %e %e  \n", a, p1, p2 ,p3, p4,p5, p1/p2);
      fprintf(fp,"%e %e %e  %e  %e %e  \n", a, p1, p2, p3, p4, p5);


}
    fclose(fp);
    return 0;
}

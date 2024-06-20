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

// output redshift
double myz = 0.;

// integration error
real epsrel = 1e-3;
// number of mass bins between 5<Log10[M]<18 for the spherical collapse solver to iterate over
double massb = 50.;
// perform 1-loop corrections? This prompts the calculation of k_star and \mathcal{E}.
bool modg = true;
// Modified - true - or LCDM (z=0) - false - input spectrum/transfer
bool mgcamb = false;


// K-mouflage parameters
double n = 2.;
double beta0 = 0.1; // glam
double K0 = 1.;
double lambda = 1.4518; //K0=1; b=0.1


// Base cosmology associated with the transfer function
double h  = 0.6774; // hubble constant
double n_s = 0.9667; // spectral index
double Omega_m = 0.3089; // total matter GLAM
double Omega_b  = 0.0486; //  baryon fraction
double pscale = 0.05; // pivot scale
double As = 2.064e-9; // As

// store beyond-LCDM parameters
double extpars[maxpars];
    // Linear and background
    extpars[0] = n;
    extpars[1] = lambda;
    extpars[2] = K0;
    extpars[3] = beta0;

double pars[8];
pars[0] = 1./(1.+myz); // scale factor
pars[1] = Omega_m;
pars[2] = 0.;
pars[5] = massb;


// Load  LCDM transfer function at z=0
ifstream flcdm("transfers/KM/glam_km_transfer_00.dat");

// Load in the transfer data
string linelcdm;
    while (getline(flcdm, linelcdm)) {      // for each line
            vector<double> lineData;           // create a new row
            double val;
            istringstream lineStream(linelcdm);
            while (lineStream >> val) {          // for each value in line
                    lineData.push_back(val);           // add to the current row
            }
            mytransl.push_back(lineData);         // add row to allData
    }

int Nkrl = mytransl.size();
int* Nktl = &Nkrl;

array Tml(*Nktl);
array kil(*Nktl);

  for(int i = 0; i< Nkrl; i++){
          kil[i] = mytransl[i][0];
          Tml[i] = mytransl[i][6];
      }



// Load cosmology classes
Cosmology Cm(h, n_s, Omega_m, Omega_b, As, pscale, kil, Tml);
// Get linear P(k) from input transfer
LinearPS_as P_l(Cm, 0.);
IOW iow;

// Load halo class witth all linear P(k)
HALO halo(Cm, P_l, P_l, P_l, P_l, epsrel);
SPT spt(Cm, P_l, epsrel);


// initialise K-mouflage scalar field and derivative
init_kmouflage_lna(pars,extpars,1000);
printf("%s \n ", "Initialised scalar field");

// initialise Hubble
iow.hubble_init(Omega_m,extpars,1000,mymodel);
printf("%s \n ", "Initialised Hubble");

//initialise spherical collapse quantities and reaction quantities
halo.initialise(pars,extpars,mgcamb,modg,mymodel);
printf("%s \n ", "Initialised collapse quantities ");

// initialise halofit parameters
halo.phinit_pseudo(pars,mgcamb);
printf("%s \n ", "Initialised halofit parameters");


real p1,p2,p3,p4,p5;

const char* output = "kmouflage_nonlinear_z0.dat";

FILE* fp = fopen(output, "w");

// number of bins between kmin and kmax
int Nk = 100;
double kmin = 0.01;
double kmax = 5.;

for(int i =0; i <Nk;  i ++) {

   real k = kmin* exp(i*log(kmax/kmin)/(Nk-1));

   p1 = halo.plinear_cosmosis(k); // Linear spectrum
   p2 = halo.reaction_nu(k,pars, mgcamb); // halo model reaction
   p3 = halo.PHALO_pseudo(k,mgcamb);//  halofit pseudo spectrum
  // p4 = spt.PLOOPn2(1, k, pars, extpars, 1e-3, mymodel); // 1-loop matter spectrum

   printf("%e %e %e %e %e  \n", k, p1,p2, p3,p3*p2); // output to terminal
   fprintf(fp,"%e %e %e %e %e \n", k, p1,p2,p3,p3*p2); // output to file : k , P_linear, R, P_pseudo, P_nl = RxP_pseudo


}
    fclose(fp);
    return 0;
}

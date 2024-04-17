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
// 14:  Cubic galileon
// 15:  QCDM

int mymodel = 14;

// target redshift
double myz = 0.;

// (sum of) neutrino mass
double mnu = 0.00;


// CG or QCDM parameters
double s = 2.;
double q = 0.5;

//output file name
const char* output = "cg_z0.dat";


// number of bins between kmin and kmax
int Nk = 300;
double kmin = 0.01;
double kmax = 5.;

// perform 1-loop corrections? This prompts the calculation of k_star and \mathcal{E}.
bool modg = true;

// Is the transfer being fed to ReACT of the target cosmology?
// If false, the transfer should be LCDM at z=0 and ReACT will rescale P_L using internally computed modified growth - see README.
// If true, the transfer function should be that of the target cosmology at the target redshift (with MG or/and massive neutrinos)
// Note that ReACT does not calculate growth factors for massive neutrino cosmologies and so the target transfer function should be supplied.
bool mgcamb = false;

// Load  LCDM transfer function at z=0
ifstream fin("transfers/cg/test_nl_transfer_00.dat");

// Base cosmology associated with the transfer function
double h  = 0.6732; // hubble constant
double n_s = 0.9655; // spectral index
double Omega_c = 0.264201; // cold dark matter
double Omega_b  = 0.0490218; //  baryon fraction
double Omega_nu = mnu/93.14/pow2(h);
double pscale = 0.05;
double As = 2.0095e-9; // initial amplitude of fluctuations

double Omega_m = Omega_c + Omega_b + Omega_nu;

// Load in the transfer data
string line;
    while (getline(fin, line)) {      // for each line
            vector<double> lineData;           // create a new row
            double val;
            istringstream lineStream(line);
            while (lineStream >> val) {          // for each value in line
                    lineData.push_back(val);           // add to the current row
            }
            mytrans.push_back(lineData);         // add row to allData
    }


int Nkr = mytrans.size();
int* Nkt = &Nkr;
array Tm(*Nkt);
array ki(*Nkt);

  for(int i = 0; i< Nkr; i++){
           ki[i] = mytrans[i][0];
           Tm[i] =  mytrans[i][6]; // total
      }

// integration error
real epsrel = 1e-3;

// number of mass bins between 5<Log10[M]<18 for the spherical collapse solver to iterate over
double massb = 50.;

// store params for passing into React functions
double pars[7];
    pars[0] = 1./(1.+myz); //  scale factor
    pars[1] = Omega_m;
    pars[2] = Omega_nu ; //  extra param
    pars[5] = massb; // number of mass bins between 5<Log10[M]<18

// store beyond-LCDM parameters
double extpars[maxpars];
    extpars[0] = s;
    extpars[1] = q;

    // Load cosmology classes
    Cosmology Cm(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tm);
    // Get linear P(k) from input transfer
    LinearPS_as P_l(Cm, 0.);
    // Load special function class (used for solving ODEs and initialising quantities)
    IOW iow;

    // Load halo class with all linear P(k)
    HALO halo(Cm, P_l, P_l, P_l, P_l, epsrel);

// Number of logarithmically spaced steps from a=3e-5 to a=10 over which to spline the Hubble rate
int hubble_steps = 2000;

// initialise Hubble
printf("%s \n", "Starting initialisation");
iow.hubble_init(Omega_m,extpars,hubble_steps,mymodel);
printf("%s \n", "Initialised background ");
//initialise spherical collapse quantities and reaction quantities
halo.initialise(pars,extpars,mgcamb,modg,mymodel);
printf("%s \n", "Initialised spherical collapse quantities and linear growth");
// initialise halofit parameters
halo.phinit_pseudo(pars,mgcamb);
printf("%s \n", "Initialised halofit parameters");

/* Output section */

/* Open output file */
FILE* fp = fopen(output, "w");

real p1,p2,p3,p4,p5,p6;

 for(int i =0; i <Nk;  i ++) {

      real k = kmin* exp(i*log(kmax/kmin)/(Nk-1));

      p1 = halo.plinear_cosmosis(k);// mypk[i][1]; // Linear spectrum
      p2 = halo.reaction_nu(k,pars, mgcamb);// mypk[i][2]; // halo model reaction
      p3 = halo.PHALO_pseudo(k,mgcamb);// mypk[i][3]; // halofit pseudo spectrum

      printf("%e %e %e %e %e  \n", k, p1,p2, p3,p3*p2); // output to terminal
      fprintf(fp,"%e %e %e %e %e \n", k, p1,p2,p3,p3*p2); // output to file : k , P_linear, R, P_pseudo, P_nl = RxP_pseudo

}

	/*close output file*/
    fclose(fp);
    return 0;
}

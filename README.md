  
# ReACT v2: ACT-io et ReACT-io 
=======

# [ReACT v2: ACT-io et ReACT-io](reactions/data/doggie.jpeg)

## This repository contains version 2 of the [ReACT](https://github.com/nebblu/ReACT) code. 

## Table of contents: 

1. [Introduction](https://github.com/nebblu/ACTio-ReACTio#introduction)
2. [Requirements](https://github.com/nebblu/ACTio-ReACTio#requirements)
3. [Installation](https://github.com/nebblu/ACTio-ReACTio#installation)
4. [Running ReACT](https://github.com/nebblu/ACTio-ReACTio#running-react)
5. [Models of gravity and dark energy](https://github.com/nebblu/ACTio-ReACTio#adding-in-models)
6. [Parameter description](https://github.com/nebblu/ACTio-ReACTio#parameters)
7. [Citation](https://github.com/nebblu/ACTio-ReACTio#citation)
8. [Notes](https://github.com/nebblu/ACTio-ReACTio#notes)
9. [What's new?](https://github.com/nebblu/ACTio-ReACTio#what-is-new)


## Introduction

ReACT is a halo model and standard perturbation theory code based on the software packages [Copter](http://mwhite.berkeley.edu/Copter/) (0905.0479) and MG-Copter (1606.02520). It allows the efficient computation of many large scale structure observables for a wide class of gravity and dark energy models. 

The package is made up of 5 folders: 

* `reactions` : contains all the c++ source code for the calculations. Specifically, this is in `reactions/src`. 
* `pyreact` : contains the wrapper functions calling the c++ code from python. 
* `cosmosis` : contains the cosmosis module for ReACT. 
* `notebooks` : example python notebooks for ReACT computations and a Mathematica file performing translations between covariant theories and modifications to the Poisson equation. 


ReACT can perform the follwing calculations for general theories beyond LCDM

* Spherical collapse  (1812.05594): `reactions/src/SCOL.cpp`.

* Halo model power spectrum (1812.05594):  `reactions/src/HALO.cpp`.

* Real and redshift space LSS 2 point statistics (1606.02520): `reactions/src/SPT.cpp`.

* Numerical perturbation theory kernels up to 3rd order (1606.02168): `reactions/src/SpericalFunctions.cpp`.

* Real space bispectrum at 1-loop level** (1808.01120): `reactions/src/BSPT.cpp`.

* Exact perturbation theory kernels up to 4th order** (1808.01120): `reactions/src/BSPTN.cpp`.

** for LCDM, nDGP and f(R) gravity only. 


## Requirements

### C++ Compiler and automake
ReACT is written in C++, so you'll need a modern C++ compiler.

### GSL
The ReACT extension makes use of a number of [gsl](http://www.gnu.org/software/gsl/) packages such as odeiv2 package which is part of the gsl library.  You will need to have a version of gsl installed (gsl 2. or later). 

### SUNDIALS

For Spherical Collapse module (`reactions/src/SCOL.cpp`) you will also need the [SUNDIALS](https://computing.llnl.gov/projects/sundials) package version 4.1.0 (also tested to work with version 5.0) **Note** that later versions of sundials may produce errors as they update their functions so please work with version 4.1.0. 

Get these packages using your package manager of choice (e.g., `homebrew` on mac OS).

One should also make sure that sundials is on your library path, for example 

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/sundials/install_dir/lib64:${LD_LIBRARY_PATH} 
```

### Python
A recent version of python should be installed. I have worked with Python 3.6.8 without issue. 

## Installation 

ReACT can be installed by following the proceeding instructions. If you experience any issues, please have a look at the [Issues tab of version 1](https://github.com/nebblu/ReACT/issues?q=is%3Aissue+is%3Aclosed). Leave an issue in this git if you experience any problems and we'll get back to you as soon as possible, or contact ben.bose@ed.ac.uk. 

We have also included a Docker file to run ReACT within a container. This is explained below. 

### General installation
The Python interface can then be installed with
```
$ python setup.py install
```
or 
```
$ pip install .
```
If you want to work on the library, install it in "developer" mode:
```
$ python setup.py develop
```
or
```
$ pip install -e .
```

### Linux specific installation 

**In all the following change the paths accordingly**

1) Make sure the following are installed: python, sundials, g++, gsl. The versions I've tested with are:

Python 3.6.8

g++ (GCC) 8.5.0

sundials-4.1.0

2) Clone ACTio-ReACTio:

```
$ git clone https://github.com/nebblu/ACTio-ReACTio.git
```

3) Add in sundials directory to pyreact/Makefile :

**Example:** 

```
LDFLAGS += -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -L/home/bose/sundials/instdir/lib64
```
```
CPPFLAGS += -I/home/bose/sundials/instdir/include
```

4) Export sundials library to LD_LIBRARY_PATH:

```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/sundials/instdir/lib64:${LD_LIBRARY_PATH}
```

5) install react

```
$ python3 setup.py develop --user
```

6) Check that everything is working. Go to the reactions/examples/ directory and try to run one of the example files. For example:

```
$ g++ -I/home/bose/react_tutorial/ACTio-ReACTio/reactions/include -L/home/bose/react_tutorial/ACTio-ReACTio/reactions/lib -lcopter -lgsl -lstdc++ halo_ps.cpp -o test
```

```
$ ./test
```

Example output:

```
$ 0 1.000000e-03 4.039658e+02 4.213977e+02 9.959108e-01 ...
```

**Note** you may also need to export the copter library path to run the example files: 

```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/bose/react_tutorial/ACTio-ReACTio/reactions/lib:${LD_LIBRARY_PATH}
```


### Docker

We have also included a Docker file, `DockerFile'. Docker allows you to run the code within what is called a container which is a dedicated environment built according to specifications given in the Docker file. Once you've installed [docker](https://www.docker.com/) you can just place this file into a folder and build the images into a container

```
docker build -t mybuild . 
```

Once this builds, you can jump into the container (which has ReACT installed) using the following command

```
docker run -v /Users/bbose/Desktop/myfolder:/home -i -t mybuild
```

This will also automatically take all the local files in `/Users/bbose/Desktop/myfolder` to the `/home` folder within the container. Anything placed in the home folder from within the container will then automatically show up locally. This lets you transfer ReACT output produced within the container to the local system. 

Alternatively, if you would like to use jupyter notebooks, then you should assign a port to your container in order to allow the connection from a host browser
```
docker run -v /Users/bbose/Desktop/myfolder:/home -it -p 8888:8888 mybuild
```

I've tested this works on the tests and example files within the `reactions/examples` and `reactions/tests` directories. You can use the following command to run them

```
g++ -I/ACTio-ReACTio/reactions/include -L/ACTio-ReACTio/reactions/lib  spt.cpp -lgsl -lcopter
```
```
./a.out
```

To run jupyter notebooks you can jump to the 'ACTio-ReACTio/notebooks'-folder and run the following command:
```
jupyter notebook --ip 0.0.0.0 --allow-root
```
You can access the notebooks through your desktops browser on http://localhost:8888 , i.e. copy-paste the second automatically generated hyperlink.  

To exit the container simply use the exit command 

```
exit 
```

## Running ReACT

### Python
The pyreact wrapper allows the C++ code (the native language of ReACT) to be called by Python code. Example jupyter notebooks that demonstrate the usage of ReACT can be found in the `notebooks` folder. Specifically we have included : 

* pyreact_demo.ipynb : Demonstrates the basic halo model reaction [2005.12184](https://arxiv.org/abs/2005.12184).
*  pyreact_demo-ext.ipynb : Demonstrates the implementation of the EFTofDE + PPF halo model reaction and compares with the DGP implementation. 
* pyreact_demo-neutrinos.ipynb : Demonstrates the implementation of the massive neutrinos halo model reaction [2105.12114](https://arxiv.org/abs/2105.12114).
* pyreact_demo-rsd.ipynb : Demonstrates the implementation of the redshift space power spectrum multipoles [1606.02520](https://arxiv.org/abs/1606.02520).
* pyreact_demo-bspt.ipynb : Demonstrates the implementation of the real space bispectrum [1808.01120](https://arxiv.org/abs/1808.01120).


###  C++
One can also run ReACT and MG-Copter in their native C++. Again, a number of example output C++ scripts have been included in `reactions/examples` as well as a number of cosmologies in `reactions/examples/transfers`. Specifically we have included : 

* actio_reactio.cpp : Example code to output the reaction and halo spectra for EFTofDE & PPF parametrisations.
* reaction_mnu.cpp : Example code to output the reaction and halo spectra for beyond LCDM physics & massive neutrinos.
* halo_ps.cpp : Example code to output the halo model powerspectrum.
* spt_rsd.cpp : Example code to output the 1-loop powerspectrum in real and redshift space.
* spt.cpp : Example code to output the 1-loop powerspectrum in real space.
* bs.cpp : Example code to output the real space bispectrum.


We can compile these examples with a command similar to : 

> gcc -I/Users/bbose/Desktop/ACTio-ReACTio/reactions/include -L/Users/bbose/Desktop/ACTio-ReACTio/reactions/lib -lcopter -lgsl -lstdc++ bs.cpp -o test

Then just run 

> ./test 

All example files require the specification of a transfer function or linear power spectrum. The given examples take transfer function inputs. You can specify the modified transfer function in two possible ways: 

1. Directly as produced using a Boltzmann solver (e.g. [MGCAMB](https://github.com/sfu-cosmo/MGCAMB)) as in reaction_mnu.cpp. **Note** that this file should be headerless. 
2. Otherwise, you must create your own cosmology file as in `reactions/examples/transfers' within which you should specify the associated z=0 LCDM transfer function file with the following two column format: {k[h/Mpc], T(k)} with the transfer function normalised to 1 at small k. This is then assumed by the code to rescale the power spectrum to the target, beyond LCDM cosmology  using internally computed LCDM and MG/DE growth factors.  

The internal flag **modcamb** tells ReACT whether or not to treat the input transfer function file as in option (1) - true value - or option (2) - false value. Default is false as in the original version of the code. 


## Adding in models

In Pyreact we currently have the following models and model parameters 

1. gr : general relativity. No extra parameters.  
2. f(r) : [Hu-Sawicki f(R)](https://arxiv.org/abs/0705.1158). **extpars[0]** = $f_{R0}$. 
3. dgp : normal branch of [DGP gravity](https://arxiv.org/abs/hep-th/0005016). **extpars[0]** = $\Omega_{rc}$. 
4. quintessence : Quintessence. **extpars[0]** = $w_0$.
5. cpl : [CPL evolving dark energy](https://arxiv.org/abs/gr-qc/0009008) w = w0 + (1-a)wa (). **extpars[0-1]**=  {$w_0,w_a$}.
6. ds : [Dark Scattering with CPL background](https://arxiv.org/abs/1605.05623). **extpars[0-2]** = {$w_0,w_a,\xi*h$}.
7. eftppf :  effective field theory of dark energy with a post parametrised friedmannian $G_{eff,non-linear}$ in spherical collapse equations. **extpars[0-9]** = {$\alpha_{k0},\alpha_{b0},\alpha_{m0},p1,...,p7$}. 
8. eftus :  effective field theory of dark energy without screening, i.e. $G_{eff, non-linear}$ = $G_{eff,linear}$. **extpars[0-2]** ={$\alpha_{k0},\alpha_{b0},\alpha_{m0}$}. 
9. eftss :  effective field theory of dark energy with superscreening , i.e. $G_{eff, non-linear}$ = $G_{N}$. **extpars[0-2]**  = {$\alpha_{k0},\alpha_{b0},\alpha_{m0}$}.  

where a `0` subscript means the value today. 

**Note** For EFTofDE models (7-9), the scale factor dependence and scale factor derivatives of the alpha functions must be specified in `reactions/src/BeyondLCDM.cpp` - see `alphai_eft` and `dalphai_eft` functions respectively. 

In C++ code this is specified as the integer of each model in the last argument of the functions , e.g. for the 1-loop real space power spectrum in f(R) gravity we would specify the following functional call 

```
PLOOPn2(1, k, pars, extpars, err, 2)
```
where the arguments are: 1 specifies a call to the 1-loop matter-matter calculation, k is the wave mode, pars hold base cosmological parameters, extpars[0] = fr0, err is the absolute error on the 1-loop integrals and 2 specifies f(R) gravity respectively. We forward to `reactions/src/examples` for more details.  

Note that all cases except 4,5,6 assume a LCDM background expansion. 

One can add in new models by simply going to the source code, `reactions/src/BeyondLCDM.cpp` and adding in a new `case n:` in all the functions with the required modified background and Poisson equation functions. The array extpars stores the theory parameters. The default size of this array is 20 and can be increased by changing the maxpars parameter in `reactions/src/BeyondLCDM.h`. 

Once you have added in the required background, linear, 1-loop and non-linear modifications you just need to re-compile the source code by going to the `reactions` directory and running 

```
make && make install
```

If you want to add in this model to the Pyreact wrapper, you will also need to add in the nth model in `pyreact/react.py`.  

## Parameters

**modcamb** tells ReACT whether or not to treat the input transfer function file as in option (1) - true value - or option (2) - false value. Default is false as in the original version of the code. 

**modg**: Tells ReACT to manually set $k_\star$ and $\mathcal{E}$ to LCDM values (1e-6 and 1. resp.). This is needed because of sensitivity of $\mathcal{E}$ to the ratio of 1-halo terms which may not be exactly equal at large scales for different cosmologies even when modified gravity isn't present. 


**Note** the cosmoSIS module has not yet been extended to include massive neutrinos. 


## Citation

When using ReACT in a publication, please acknowledge the code by citing the relevant papers from the following:

[arXiv:1812.05594](https://arxiv.org/abs/1812.05594) : "On the road to percent accuracy I: non-linear reaction of the matter power spectrum to dark energy and modified gravity"

[arXiv:1909.02561](https://arxiv.org/abs/1909.02561) : "On the road to percent accuracy III: non-linear reaction of the matter power spectrum to massive neutrinos"

[arXiv:2005.12184](https://arxiv.org/abs/2005.12184) : "On the road to per-cent accuracy IV:  ReACT -- computing the non-linear power spectrum beyond LCDM"

[arXiv:2105.12114](https://arxiv.org/abs/2105.12114) : "On the road to percent accuracy V: the non-linear power spectrum beyond LCDM with massive neutrinos and baryonic feedback"

[arXiv:2111.13598](https://arxiv.org/abs/2111.13598)  : "On the road to per cent accuracy VI: the non-linear power spectrum for interacting dark energy with baryonic feedback and massive neutrinos
"

[arXiv:1606.02520](https://arxiv.org/abs/1606.02520) : "A Perturbative Approach to the Redshift Space Power Spectrum: Beyond the Standard Model"

[arXiv:1808.01120](https://arxiv.org/abs/1808.01120) : "The one-loop matter bispectrum as a probe of gravity and dark energy"



Respective bibtex entries:

```
 @article{Cataneo:2018cic,
   author = "Cataneo, Matteo and Lombriser, Lucas and Heymans, Catherine and Mead, Alexander and Barreira, Alexandre and Bose, Sownak and Li, Baojiu",
    title = "{On the road to percent accuracy: non-linear reaction of the matter power spectrum to dark energy and modified gravity}",
     eprint = "1812.05594",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1093/mnras/stz1836",
    journal = "Mon. Not. Roy. Astron. Soc.",
    volume = "488",
    number = "2",
    pages = "2121--2142",
    year = "2019"
}
```

```
@article{Cataneo:2019fjp,
    author = "Cataneo, Matteo and Emberson, J. D. and Inman, Derek and Harnois-Deraps, Joachim and Heymans, Catherine",
    title = "{On the road to per cent accuracy \textendash{} III. Non-linear reaction of the matter power spectrum to massive neutrinos}",
    eprint = "1909.02561",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1093/mnras/stz3189",
    journal = "Mon. Not. Roy. Astron. Soc.",
    volume = "491",
    number = "3",
    pages = "3101--3107",
    year = "2020"
}
```



```
@article{Bose:2020wch,
    author = "Bose, Benjamin and Cataneo, Matteo and Tröster, Tilman and Xia, Qianli and Heymans, Catherine and Lombriser, Lucas",
    title = "{On the road to per-cent accuracy IV: ReACT -- computing the non-linear power spectrum beyond $\Lambda$CDM}",
    eprint = "2005.12184",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    month = "5",
    year = "2020"
}
```

```
@article{Bose:2021mkz,
    author = "Bose, Benjamin and Wright, Bill S. and Cataneo, Matteo and Pourtsidou, Alkistis and Giocoli, Carlo and Lombriser, Lucas and McCarthy, Ian G. and Baldi, Marco and Pfeifer, Simon and Xia, Qianli",
    title = "{On the road to percent accuracy V: the non-linear power spectrum beyond $\Lambda$CDM with massive neutrinos and baryonic feedback}",
    eprint = "2105.12114",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    month = "5",
    year = "2021"
}
```

```
@article{Carrilho:2021rqo,
    author = "Carrilho, Pedro and Carrion, Karim and Bose, Benjamin and Pourtsidou, Alkistis and Hidalgo, Juan Carlos and Lombriser, Lucas and Baldi, Marco",
    title = "{On the road to~per\,cent accuracy VI: the non-linear power spectrum for interacting dark energy with baryonic feedback and massive neutrinos}",
    eprint = "2111.13598",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1093/mnras/stac641",
    journal = "Mon. Not. Roy. Astron. Soc.",
    volume = "512",
    number = "3",
    pages = "3691--3702",
    year = "2022"
}
```

```
@article{Bose:2016qun,
    author = "Bose, Benjamin and Koyama, Kazuya",
    title = "{A Perturbative Approach to the Redshift Space Power Spectrum: Beyond the Standard Model}",
    eprint = "1606.02520",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    doi = "10.1088/1475-7516/2016/08/032",
    journal = "JCAP",
    volume = "08",
    pages = "032",
    year = "2016"
}
```

```
@article{Bose:2018zpk,
    author = "Bose, Benjamin and Taruya, Atsushi",
    title = "{The one-loop matter bispectrum as a probe of gravity and dark energy}",
    eprint = "1808.01120",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    reportNumber = "YITP-18-88",
    doi = "10.1088/1475-7516/2018/10/019",
    journal = "JCAP",
    volume = "10",
    pages = "019",
    year = "2018"
}
```

## Notes 

### Parameter ranges (Updated: 22/03/21)
* To optimise root finding within the spherical collapse portion of the code, the maximum redshift that one can solve the reaction for currently is z=2.5. 
* There are some current issues in the wCDM part of the code. Namely for very particular values of w0 and wa in the CPL evolving dark energy case, the spherical collapse library cannot solve the virial theorem. We advise sticking to the ranges 
-1.3<w0<-0.7 and -1.5<wa<0.6 to avoid these issues. 
* Spherical collapse will not solve for values of sigma_8(z=0)< 0.55 or 1.4<sigma_8(z=0).  
* We have tested the massive neutrino code for m_nu=0.48eV (Omega_nu = 0.0105) without issue in the absence of modified gravity or evolving dark energy. 
* Spherical collapse may encounter issues for high fr0 values in combination with high Omega_nu. We have tested the code for Log[fr0]=-4 with m_nu=0.3eV with success but do not recommend higher values than these in combination.  

### Miscellaneous (Updated: 28/05/20)
* If errors in spherical collapse are experienced for non-f(R) theories, try setting yenvf=0 in the scol_init function in reactions/src/HALO.cpp.
* Note if using the stand-alone version of ReACT, the reaction may have deviations away from unity of the order of ~0.1-0.3% for k<1e-3. Pyreact automatically sets it to unity at such large scales. 


## What is new 

We have implemented the following to v.2:

* Added in EFTofDE + PPF model for model independent predictions. 
* Cleaned up directories, added heavy commenting in source code and notebooks and re ordered some source code for clarity. 
* Merged [Dark Scattering model](https://github.com/PedroCarrilho/ReACT/tree/react_with_interact_baryons). 
* Split off beyond LCDM functions to BeyondLCDM.cpp in src directory for easily adding in new models. 
* Upgraded all beyond LCDM functions to use an n-dimensional array for parameters allowing for arbitrary number of theory parameters to be used. 
* Upgraded redshift space functions to use an n-dimensional array for rsd and bias parameters.
* Optimised RSD multipole computation. 
* Created python wrapper for RSD multipoles. 
* Created new example python notebooks for rsd multipoles and EFTofDE calculations. 

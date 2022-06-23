import os
import ctypes as ct
import numpy as np

def array_ctype(ndim, dtype=np.float64, flags="C"):
    return [ct.POINTER(ct.c_int)]*ndim + [np.ctypeslib.ndpointer(ndim=ndim, dtype=dtype, flags=flags)]

def array_arg(a):
    arr = a
    return (*(ct.c_int(s) for s in arr.shape), arr)

class ReACT:
    libname = "libreact_wrapper.so"
    module_name = "reaction_module"

    def __init__(self):
        self.load_lib()

    def load_lib(self, path=None):
        if path is None:
            path = os.path.dirname(__file__)
        libpath = os.path.abspath(os.path.join(path, self.libname))
        self.lib = ct.CDLL(libpath)

    def get_function(self, name, c_bind=True):
        if c_bind:
            return getattr(self.lib, name)
        else:
            return getattr(self.lib, f"__{self.module_name}_MOD_{name}")

    def test_func(self, a):
        f = self.get_function("test_func")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=2, dtype=np.float64)]

        r = f(*array_arg(a))
        return r

# Compute reaction without massive neutrinos
    def compute_reaction(self, h, n_s, omega_m, omega_b, sigma_8,
                               z, k, Pk,
                               model="f(r)", fR0=None, Omega_rc=None, w=None, wa=None,
                               is_transfer=False, mass_loop=30,
                               verbose=True):


        if max(z) > 2.5:
            raise ValueError("ReACT is unstable above z=2.5, try limiting the range of z values.")

        if len(z) > 1 and z[0] > z[-1]:
            raise ValueError("The z array needs to be ordered from low to high redshift.")

        if model.lower() == "f(r)":
            reaction_model = 2
            modified_gravity_param = fR0
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "dgp":
            reaction_model = 3
            modified_gravity_param = Omega_rc
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "gr":
            reaction_model = 1
            modified_gravity_param = 0.0
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "quintessence":
            reaction_model = 4
            modified_gravity_param = w
            modified_gravity_param2 = 0.0
            modified_gravity_param3 = 0.0
        elif model.lower() == "cpl":
            reaction_model = 5
            modified_gravity_param = w
            modified_gravity_param2 = wa
            modified_gravity_param3 = 0.0
        else:
            raise ValueError(f"model '{model}' not supported.")

        if modified_gravity_param is None:
            raise ValueError("fR0, Omega_rc or w0 need to be specified.")

        f = self.get_function("compute_reaction")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=1, dtype=np.float64), # P(k, z=0)
                      *array_ctype(ndim=1, dtype=np.float64), # k
                      *array_ctype(ndim=1, dtype=np.float64), # z
                      ct.POINTER(ct.c_bool),       # is_transfer
                      ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # n_s
                      ct.POINTER(ct.c_double),     # omega_m
                      ct.POINTER(ct.c_double),     # omega_b
                      ct.POINTER(ct.c_double),     # sigma_8
                      ct.POINTER(ct.c_double),     # model parameter 1 : for f(R) this is fr0, for dgp this is Omega_rc, for CPL or quintessence this is w0
                      ct.POINTER(ct.c_double),     # model parameter 2 : for CPL this is wa
                      ct.POINTER(ct.c_double),     # model parameter 3
                      ct.POINTER(ct.c_int),        # mass_loop
                      ct.POINTER(ct.c_int),        # model (1: GR, 2: f(R), 3; DGP, 4: quintessence, 5: CPL)
                      *array_ctype(ndim=2, dtype=np.float64), # reaction (output)
                      *array_ctype(ndim=2, dtype=np.float64), # linear MG power spectrum (output)
                      np.ctypeslib.ndpointer(ndim=1, dtype=np.float64, flags="C"),     # modified sigma_8 storage variable
                      ct.POINTER(ct.c_int),        # verbose
                     ]
        reaction = np.zeros((len(k), len(z)), dtype=Pk.dtype, order="C")
        p_lin = np.zeros((len(k), len(z)), dtype=Pk.dtype, order="C")
        sigma8 = np.zeros(1, dtype=Pk.dtype, order="C")

        r =   f(*array_arg(np.ascontiguousarray(Pk, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(k, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(z[::-1], dtype=np.float64)), # ReACT expect z order from high to low
                ct.c_bool(is_transfer),
                ct.c_double(h), ct.c_double(n_s), ct.c_double(omega_m), ct.c_double(omega_b), ct.c_double(sigma_8),
                ct.c_double(modified_gravity_param),
                ct.c_double(modified_gravity_param2),
                ct.c_double(modified_gravity_param3),
                ct.c_int(mass_loop),
                ct.c_int(reaction_model),
                *array_arg(reaction),
                *array_arg(p_lin),
                sigma8,
                ct.c_int(verbose),
                )
        if r != 0:
            string_type = ct.c_char * ct.c_int.in_dll(self.lib, "ERROR_MESSAGE_LEN").value
            error_message = string_type.in_dll(self.lib, "error_message").value.decode()
            raise RuntimeError(f"ReACT code terminated with an error: {error_message}")
        # Get into CAMB ordering (z, k), with increasing z
        #Tilman: to check output of modsig8
        return reaction[:,::-1].T, p_lin[:,::-1].T, sigma8[0]



  # Compute reaction without massive neutrinos with additional parameters
    def compute_reaction_ext(self, h, n_s, omega_m, omega_b, sigma_8,
                               z, k, Pk,
                               model, extpars,
                               compute_pseudo=False,
                               is_transfer=False, mass_loop=30,
                               verbose=True):


        if max(z) > 2.5:
            raise ValueError("ReACT is unstable above z=2.5, try limiting the range of z values.")

        if len(z) > 1 and z[0] > z[-1]:
            raise ValueError("The z array needs to be ordered from low to high redshift.")

        if model.lower() == "gr":
            reaction_model = 1
        elif model.lower() == "f(r)":
            reaction_model = 2
        elif model.lower() == "dgp":
            reaction_model = 3
        elif model.lower() == "quintessence":
            reaction_model = 4
        elif model.lower() == "cpl":
            reaction_model = 5
        elif model.lower() == "ds":
            reaction_model = 6
        elif model.lower() == "eftppf":
            reaction_model = 7
        elif model.lower() == "eftus":
            reaction_model = 8
        elif model.lower() == "eftss":
            reaction_model = 9
        else:
            raise ValueError(f"model '{model}' not supported.")

        if extpars[0] is None:
            raise ValueError("First extended parameter need to be specified.")

        f = self.get_function("compute_reaction_ext")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=1, dtype=np.float64), # P(k, z=0)
                      *array_ctype(ndim=1, dtype=np.float64), # k
                      *array_ctype(ndim=1, dtype=np.float64), # z
                      ct.POINTER(ct.c_bool),       # is_transfer
                      ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # n_s
                      ct.POINTER(ct.c_double),     # omega_m
                      ct.POINTER(ct.c_double),     # omega_b
                      ct.POINTER(ct.c_double),     # sigma_8
                      *array_ctype(ndim=1, dtype=np.float64), # extended parameters
                      ct.POINTER(ct.c_int),        # mass_loop
                      ct.POINTER(ct.c_int),        # model (1: GR, 2: f(R), 3; DGP, 4: quintessence, 5: CPL)
                      *array_ctype(ndim=2, dtype=np.float64), # reaction (output)
                      *array_ctype(ndim=2, dtype=np.float64), # linear MG power spectrum (output)
                      *array_ctype(ndim=2, dtype=np.float64), # pseudo spectrum
                      np.ctypeslib.ndpointer(ndim=1, dtype=np.float64, flags="C"),     # modified sigma_8 storage variable
                      ct.POINTER(ct.c_bool),       # compute_pseudo
                      ct.POINTER(ct.c_int),        # verbose
                     ]
        reaction = np.zeros((len(k), len(z)), dtype=Pk.dtype, order="C")
        p_lin = np.zeros((len(k), len(z)), dtype=Pk.dtype, order="C")
        pseudo = np.zeros((len(k), len(z)), dtype=Pk.dtype, order="C")
        sigma8 = np.zeros(1, dtype=Pk.dtype, order="C")

        r =   f(*array_arg(np.ascontiguousarray(Pk, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(k, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(z[::-1], dtype=np.float64)), # ReACT expect z order from high to low
                ct.c_bool(is_transfer),
                ct.c_double(h), ct.c_double(n_s), ct.c_double(omega_m), ct.c_double(omega_b), ct.c_double(sigma_8),
                *array_arg(np.ascontiguousarray(extpars, dtype=np.float64)),
                ct.c_int(mass_loop),
                ct.c_int(reaction_model),
                *array_arg(reaction),
                *array_arg(p_lin),
                *array_arg(pseudo),
                sigma8,
                ct.c_bool(compute_pseudo),
                ct.c_int(verbose),
                )
        if r != 0:
            string_type = ct.c_char * ct.c_int.in_dll(self.lib, "ERROR_MESSAGE_LEN").value
            error_message = string_type.in_dll(self.lib, "error_message").value.decode()
            raise RuntimeError(f"ReACT code terminated with an error: {error_message}")
        # Get into CAMB ordering (z, k), with increasing z
        #Tilman: to check output of modsig8
        return reaction[:,::-1].T, p_lin[:,::-1].T, sigma8[0], pseudo[:,::-1].T



    def compute_reaction_nu_ext(self, h, n_s, omega_m, omega_b, omega_nu, As,
                               z, k, Tm, Tcb, klcdm, Tcblcdm,
                               pscale,
                               model, extpars,
                               compute_pseudo=False,
                               is_transfer=False, mass_loop=30,
                               verbose=True):

        if max(z) > 2.5:
            raise ValueError("ReACT is unstable above z=2.5, try limiting the range of z values.")

        if len(z) > 1 and z[0] > z[-1]:
            raise ValueError("The z array needs to be ordered from low to high redshift.")

        if model.lower() == "gr":
            reaction_model = 1
        elif model.lower() == "f(r)":
            reaction_model = 2
        elif model.lower() == "dgp":
            reaction_model = 3
        elif model.lower() == "quintessence":
            reaction_model = 4
        elif model.lower() == "cpl":
            reaction_model = 5
        elif model.lower() == "ds":
            reaction_model = 6
        elif model.lower() == "eftppf":
            reaction_model = 7
        elif model.lower() == "eftus":
            reaction_model = 8
        elif model.lower() == "eftss":
            reaction_model = 9
        else:
            raise ValueError(f"model '{model}' not supported.")


        if extpars[0] is None:
            raise ValueError("First extended parameter need to be specified.")


        f = self.get_function("compute_reaction_nu_ext")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=1, dtype=np.float64), # full matter transfer
                      *array_ctype(ndim=1, dtype=np.float64), # k
                      *array_ctype(ndim=1, dtype=np.float64), # z
                      *array_ctype(ndim=1, dtype=np.float64), # CDM + baryon transfer
                      *array_ctype(ndim=1, dtype=np.float64), # LCDM CDM+baryon transfer
                      *array_ctype(ndim=1, dtype=np.float64), # k for LCDM transfer
                      ct.POINTER(ct.c_bool),       # is_transfer
                      ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # n_s
                      ct.POINTER(ct.c_double),     # omega_m
                      ct.POINTER(ct.c_double),     # omega_b
                      ct.POINTER(ct.c_double),     # omega_nu
                      ct.POINTER(ct.c_double),     # A_s
                      ct.POINTER(ct.c_double),     # pscale
                      *array_ctype(ndim=1, dtype=np.float64), # extended parameters
                      ct.POINTER(ct.c_int),        # mass_loop
                      ct.POINTER(ct.c_int),        # model (1: GR, 2: f(R), 3; DGP, 4: quintessence, 5: CPL)
                      *array_ctype(ndim=2, dtype=np.float64), # reaction (output)
                      *array_ctype(ndim=2, dtype=np.float64), # linear MG power spectrum (output)
                      *array_ctype(ndim=2, dtype=np.float64), # pseudo spectrum
                      np.ctypeslib.ndpointer(ndim=1, dtype=np.float64, flags="C"),     # modified sigma_8 storage variable
                      ct.POINTER(ct.c_bool), # compute_pseudo
                      ct.POINTER(ct.c_int),        # verbose
                     ]

        reaction = np.zeros((len(k), len(z)), dtype=k.dtype, order="C")
        p_lin = np.zeros((len(k), len(z)), dtype=k.dtype, order="C")
        pseudo = np.zeros((len(k), len(z)), dtype=k.dtype, order="C")
        sigma8 = np.zeros(1, dtype=k.dtype, order="C")


        r =   f(*array_arg(np.ascontiguousarray(Tm, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(k, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(z[::-1], dtype=np.float64)), # ReACT expect z order from high to low
                *array_arg(np.ascontiguousarray(Tcb, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(Tcblcdm, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(klcdm, dtype=np.float64)),
                ct.c_bool(is_transfer),
                ct.c_double(h), ct.c_double(n_s), ct.c_double(omega_m), ct.c_double(omega_b), ct.c_double(omega_nu),
                ct.c_double(As), ct.c_double(pscale),
                *array_arg(np.ascontiguousarray(extpars, dtype=np.float64)),
                ct.c_int(mass_loop),
                ct.c_int(reaction_model),
                *array_arg(reaction),
                *array_arg(p_lin),
                *array_arg(pseudo),
                sigma8,
                ct.c_bool(compute_pseudo),
                ct.c_int(verbose),
                )
        if r != 0:
            string_type = ct.c_char * ct.c_int.in_dll(self.lib, "ERROR_MESSAGE_LEN").value
            error_message = string_type.in_dll(self.lib, "error_message").value.decode()
            raise RuntimeError(f"ReACT code terminated with an error: {error_message}")

        return reaction[:,::-1].T, p_lin[:,::-1].T, sigma8[0], pseudo[:,::-1].T




  # Compute reaction without massive neutrinos with additional parameters
    def compute_multipoles(self, h, n_s, omega_m, omega_b, sigma_8,
                               z, k, Pk, kout,
                               model, rsd_model, whichmulti,
                               extpars, rsdpars, biaspars, errpars,
                               is_transfer=False,
                               verbose=True):


        if z > 2.5:
            raise ValueError("Code is unstable above z=2.5.")

        if model.lower() == "gr":
            reaction_model = 1
        elif model.lower() == "f(r)":
            reaction_model = 2
        elif model.lower() == "dgp":
            reaction_model = 3
        elif model.lower() == "quintessence":
            reaction_model = 4
        elif model.lower() == "cpl":
            reaction_model = 5
        elif model.lower() == "ds":
            reaction_model = 6
        elif model.lower() == "eftppf":
            reaction_model = 7
        elif model.lower() == "eftus":
            reaction_model = 8
        elif model.lower() == "eftss":
            reaction_model = 9
        else:
            raise ValueError(f"model '{model}' not supported.")

        if extpars[0] is None:
            raise ValueError("First extended parameter need to be specified.")

        f = self.get_function("compute_multipoles")
        f.restype = np.int
        f.argtypes = [*array_ctype(ndim=1, dtype=np.float64), # P(k, z=0)
                      *array_ctype(ndim=1, dtype=np.float64), # k
                      *array_ctype(ndim=1, dtype=np.float64), # kout
                      ct.POINTER(ct.c_bool),       # is_transfer
                      ct.POINTER(ct.c_double),     # z
                      ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # n_s
                      ct.POINTER(ct.c_double),     # omega_m
                      ct.POINTER(ct.c_double),     # omega_b
                      ct.POINTER(ct.c_double),     # sigma_8
                      *array_ctype(ndim=1, dtype=np.float64), # extended parameters
                      *array_ctype(ndim=1, dtype=np.float64), # bias parameters
                      *array_ctype(ndim=1, dtype=np.float64), # rsd parameters
                      *array_ctype(ndim=1, dtype=np.float64), # error parameters
                      ct.POINTER(ct.c_int),        # rsd model (0: kaiser, 1: TNS w qbias, 2: TNS w eul bias, 3: 1-loop SPT)
                      ct.POINTER(ct.c_int),        # model (1: GR, 2: f(R), 3; DGP, 4: quintessence, 5: CPL)
                      ct.POINTER(ct.c_int),        # how many multipoles? (1: p0, 2: p0 + p2, 3: p0 + p2 + p4)
                      *array_ctype(ndim=1, dtype=np.float64), # p0 (output)
                      *array_ctype(ndim=1, dtype=np.float64), # p2 (output)
                      *array_ctype(ndim=1, dtype=np.float64), # p4 (output)
                      *array_ctype(ndim=1, dtype=np.float64), # linear growth factor
                      *array_ctype(ndim=1, dtype=np.float64), # linear growth rate
                      ct.POINTER(ct.c_int),        # verbose
                     ]
        myp0 = np.zeros(len(kout), dtype=Pk.dtype, order="C")
        myp2 = np.zeros(len(kout), dtype=Pk.dtype, order="C")
        myp4 = np.zeros(len(kout), dtype=Pk.dtype, order="C")
        mygf = np.zeros(len(kout), dtype=Pk.dtype, order="C")
        mygr = np.zeros(len(kout), dtype=Pk.dtype, order="C")

        r =   f(*array_arg(np.ascontiguousarray(Pk, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(k, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(kout, dtype=np.float64)),
                ct.c_bool(is_transfer),
                ct.c_double(z),
                ct.c_double(h), ct.c_double(n_s), ct.c_double(omega_m), ct.c_double(omega_b), ct.c_double(sigma_8),
                *array_arg(np.ascontiguousarray(extpars, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(biaspars, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(rsdpars, dtype=np.float64)),
                *array_arg(np.ascontiguousarray(errpars, dtype=np.float64)),
                ct.c_int(rsd_model),
                ct.c_int(reaction_model),
                ct.c_int(whichmulti),
                *array_arg(myp0),
                *array_arg(myp2),
                *array_arg(myp4),
                *array_arg(mygf),
                *array_arg(mygr),
                ct.c_int(verbose),
                )
        if r != 0:
            string_type = ct.c_char * ct.c_int.in_dll(self.lib, "ERROR_MESSAGE_LEN").value
            error_message = string_type.in_dll(self.lib, "error_message").value.decode()
            raise RuntimeError(f"ReACT code terminated with an error: {error_message}")
            #Output multipoles and growth
        return myp0[:], myp2[:], myp4[:], mygf[:], mygr[:]

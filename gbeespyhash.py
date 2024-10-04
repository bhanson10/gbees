# gbeespy.py, https://github.com/bhanson10/gbees-hash
# Copyright 2024 by Benjamin Hanson, published under BSD 3-Clause License.

import ctypes as ct

lib = ct.CDLL("../../gbees-hash.so")
TOL = 1E-8

class Meas(ct.Structure):
    _fields_ = [
        ("dim", ct.POINTER(ct.c_int)),
        ("mean", ct.POINTER(ct.c_double)),
        ("cov", ct.POINTER(ct.POINTER(ct.c_double))),
        ("T", ct.c_double)
    ]

def Meas_create(dim, M_DIR, M_FILE):
    lib.Meas_create.argtypes = [ct.c_int, ct.c_char_p, ct.c_char_p]
    lib.Meas_create.restype = Meas

    M_DIR = ct.c_char_p(M_DIR.encode('utf-8'))
    M_FILE = ct.c_char_p(M_FILE.encode('utf-8'))
    return lib.Meas_create(dim, M_DIR, M_FILE)

class Grid(ct.Structure):
    _fields_ = [
        ("dim", ct.c_int),
        ("thresh", ct.c_double),
        ("t", ct.c_double),
        ("dt", ct.c_double),
        ("center", ct.POINTER(ct.c_double)),
        ("dx", ct.POINTER(ct.c_double)),
        ("hi_bound", ct.POINTER(ct.c_double)),
        ("lo_bound", ct.POINTER(ct.c_double))
    ]

def Grid_create(dim, thresh, center, dx):
    lib.Grid_create.argtypes = [ct.c_int, ct.c_double, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)]
    lib.Grid_create.restype = Grid

    dx = (ct.c_double * dim)(*dx)
    return lib.Grid_create(dim, thresh, center, dx)

class Traj(ct.Structure):
    _fields_ = [
        ("coef", ct.POINTER(ct.c_double))
    ]

def Traj_create(n, coef):
    lib.Traj_create.argtypes = [ct.c_int, ct.POINTER(ct.c_double)]
    lib.Traj_create.restype = Traj

    coef = (ct.c_double * n)(*coef)
    return lib.Traj_create(n, coef)

def run_gbees(f, h, BOUND_f, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, CAPACITY, DIM_h, OUTPUT, RECORD, MEASURE, BOUNDS, COLLISIONS):

    c_P_DIR = P_DIR.encode('utf-8')
    c_P_DIR = ct.create_string_buffer(c_P_DIR)
    c_P_DIR = ct.cast(c_P_DIR, ct.POINTER(ct.c_char))

    c_M_DIR = M_DIR.encode('utf-8')
    c_M_DIR = ct.create_string_buffer(c_M_DIR)
    c_M_DIR = ct.cast(c_M_DIR, ct.POINTER(ct.c_char))

    SYS_CALLBACK_FUNC = ct.CFUNCTYPE(None, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_double, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double))
    BOUND_CALLBACK_FUNC = ct.CFUNCTYPE(ct.c_double, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double))

    @SYS_CALLBACK_FUNC
    def c_f(xk, x, t, dx, coef):
        v = f(x, t, dx, coef)
        for i in range(len(v)):
            xk[i] = v[i]

    @SYS_CALLBACK_FUNC
    def c_h(y, x, t, dx, coef):
        v = h(x, t, dx, coef)
        for i in range(len(v)):
            y[i] = v[i]
    
    if(BOUNDS):
        @BOUND_CALLBACK_FUNC
        def c_BOUND_f(x, coef):
            J = BOUND_f(x, coef)
            return J
        
        lib.run_gbees.argtypes = [SYS_CALLBACK_FUNC, SYS_CALLBACK_FUNC, BOUND_CALLBACK_FUNC, Grid, Meas, Traj, ct.c_char_p, ct.c_char_p, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_bool, ct.c_bool, ct.c_bool, ct.c_bool, ct.c_bool]
        lib.run_gbees.restype = None
    else:
        c_BOUND_f = None
        lib.run_gbees.argtypes = [SYS_CALLBACK_FUNC, SYS_CALLBACK_FUNC, ct.c_char_p, Grid, Meas, Traj, ct.c_char_p, ct.c_char_p, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_bool, ct.c_bool, ct.c_bool, ct.c_bool, ct.c_bool]
        lib.run_gbees.restype = None

    lib.run_gbees(c_f, c_h, c_BOUND_f, G, M, T, c_P_DIR, c_M_DIR, ct.c_int(NUM_DIST), ct.c_int(NUM_MEAS), ct.c_int(DEL_STEP), ct.c_int(OUTPUT_FREQ), ct.c_int(CAPACITY), ct.c_int(DIM_h), ct.c_bool(OUTPUT), ct.c_bool(RECORD), ct.c_bool(MEASURE), ct.c_bool(BOUNDS), ct.c_bool(COLLISIONS))
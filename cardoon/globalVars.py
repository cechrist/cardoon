"""
:mod:`globalVars` -- Global simulator options and physical constants
--------------------------------------------------------------------

.. moduleauthor:: Carlos Christoffersen

Physical constants are from http://physics.nist.gov/cuu/Constants/

Usage example::

    from globalVars import const, glVar

    # Calculate thermal voltage (Vt) at default temperature
    vt = const.k * glVar.temp / const.q

**Important Note:** The ``glVar`` name really refers to the variables
in the ``.options`` line, **not** the netlist variables defined in
``.vars``. Those are kept in a dictionary in ParamSet.

"""

from paramset import ParamSet
import numpy as np

# Physical constants
constDict = dict(
    k = ('Boltzmann constant', 'J K^{-1}', float, 1.3806488e-23),
    q = ('Elementary charge', 'C', float, 1.602176565e-19),
    epsilon0 = ('Permittivity of free space', 'F m^{-1}', 
     float, 8.854187817e-12),
    mu0 = ('Permeability of free space', 'H m^{-1}', float,
     np.pi*4e-7),
    c0 = ('Speed of light in free space', 'm s^{-1}', float, 2.99792458e8),
    T0 = ('Zero degree Celsius temperature', 'K', float, 273.15),
    epSi = ('Permitivity of silicon', '', float, 104.5e-12),
    epOx = ('Permitivity of silicon oxide', '', float, 34.5e-12),
    Np2dB = ('Neper to dB conversion constant', 'dB/Np', 
             float, 8.6858896380650368)
    )

const = ParamSet(constDict)
const.set_attributes()

globDict = dict(
    temp = ('Ambient temperature', 'C', float, 27.),
    abstol = ('Absolute tolerance for nodal variables', 'nodal unit', 
              float, 1e-7),
    reltol = ('Relative tolerance for nodal variables', '', float, 1e-4),
    maxiter = ('Maximum number of Newton iterations', '', int, 40),
    softiter = ('Number of Newton iterations if no progress', '', int, 20),
    maxdelta = ('Maximum allowed change in one Newton iteration', 
                'nodal unit', float, 20.),
    errfunc = ("Enable additional test for error function in Newton's method", 
               'bool', bool, False),
    sparse = ("Use sparse matrices in analyses if possible", 'bool', 
              bool, True),
    verbose = ('Verbose output', '', bool, False),
    gyr = ('Default gain in internal gyrators', 'S', float, 1e-3)
)

glVar = ParamSet(globDict)
glVar.set_attributes()


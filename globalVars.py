"""
Global simulator variables, including global settings and physical
constants.

Physical constants are from http://physics.nist.gov/cuu/Constants/

Use const.mu0 to access constants and glVar.temp for global
variables.

--------------------------------------------------------------------

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
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
    epOx = ('Permitivity of silicon oxide', '', float, 34.5e-12)
    )

const = ParamSet(constDict)
const.set_attributes()

globDict = dict(
    temp = ('Ambient temperature', 'C', float, 27.),
    abstol = ('Absolute tolerance', '-', float, 1e-8),
    reltol = ('Relative tolerance', '-', float, 1e-5),
    gyr = ('Default gain in internal gyrators', 'S', float, 1e-2)
)

glVar = ParamSet(globDict)
glVar.set_attributes()


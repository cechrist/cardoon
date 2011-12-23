"""
:mod:`svbjt` -- State-Variable-Based Bipolar Transistor
-------------------------------------------------------

.. module:: svbjt
.. moduleauthor:: Carlos Christoffersen

"""

import svbjti
from extbjt import extrinsic_bjt

Device = extrinsic_bjt(svbjti.Device)

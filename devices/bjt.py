"""
:mod:`bjt` -- Bipolar Transistor
--------------------------------

.. module:: bjt
.. moduleauthor:: Carlos Christoffersen

"""

import bjti
from extbjt import extrinsic_bjt

Device = extrinsic_bjt(bjti.Device)

"""
This package contains all device models.  Device classes are imported
into a dictionary. Keys are the device types. So to create a new
device just use:

devices.devClass['devType']('instancename')

Examples: 

from devices import devClass
# Create device instance
m1 = devClass['mosacm']('m1n')

Refer to EmptyDev.py for a template how to build a device model.

----------------------------------------------------------------------
Copyright Carlos Christoffersen <c.christoffersen@ieee.org>

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""
import autoThermal

# Regular 'netlist' elements must be listed here
netElemList = ['mosACM', 'resistor', 'capacitor', 'inductor', 
               'idc', 'vdc', 'diode', 'svdiode', 'mosEKV', 'bjt']

# Add here any modules to be imported in addition to netElemList
__all__ = netElemList 

devClass = {}
parType = {}
for modname in netElemList:
    module = __import__(modname, globals(), locals())
    devClass[module.Device.devType] = module.Device
    # If module has defined the makeAutoThermal flag then attempt
    # to create electrothermal device
    if module.Device.makeAutoThermal:
        Device = autoThermal.thermal_device(module.Device)
        devClass[Device.devType] = Device

        
        




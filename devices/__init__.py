"""
devices package
---------------

This package contains a library with device models.  Device classes
are imported into a dictionary. Keys are the device types. So to
create a new device just use::

    devices.devClass['devType']('instancename')

Examples::

    from devices import devClass
    # Create device instance
    m1 = devClass['mosacm']('m1n')

"""
import autoThermal

# Regular 'netlist' elements must be listed here
netElemList = ['mosACM', 'resistor', 'capacitor', 'inductor', 
               'idc', 'vdc', 'diode', 'svdiode', 'mosEKV', 'bjt', 'svbjt']

# Add here any modules to be imported in addition to netElemList
__all__ = netElemList + ['cppaddev']

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

        
        




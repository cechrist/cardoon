"""
:mod:`devices` -- Device Library
--------------------------------

The ``devices`` package contains a library with device models.  Device
classes are imported into a dictionary. Keys are the device types. To
create a new device use the following::

    devices.devClass['devType']('instancename')

Example (from python)::

    from devices import devClass
    # Create device instance
    m1 = devClass['mosacm']('m1n')

Example (from netlist)::

    mosacm:m1n <nodes> <parameters>

Making a device model to be recognized by this package
++++++++++++++++++++++++++++++++++++++++++++++++++++++

Suppose the new model is implemented in a file named
``newmodel.py``. Save this file in the ``devices`` directory and edit
``devices/__init__.py``. Add your module name to ``netElemList`` as
shown below::

    # Regular 'netlist' elements must be listed here
    netElemList = ['mosACM', 'resistor', 'capacitor', 'inductor', 'idc', 'vdc', 
                   'diode', 'svdiode', 'mosEKV', 'bjt', 'svbjt', 'newelem']

That's all!

"""
import autoThermal

# Regular 'netlist' elements must be listed here
netElemList = ['resistor', 'capacitor', 'inductor', 'idc', 'vdc', 
               'diode', 'svdiode', 'mosEKV', 'mosACM', 'bjt', 'svbjt', 
               'tlinpy4', 'tlinps4', 'isin', 'vsin', 'vpulse']

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

        
        




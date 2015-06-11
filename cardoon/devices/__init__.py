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
netElemList = ['resistor', 'capacitor', 'nonidealcapacitor', 'inductor', 
               'nonidealinductor', 'memristor', 'memductor',
               'gyrator', 'vccs',
               'idc', 'vdc', 'isin', 'vsin', 'ipulse', 'vpulse',
               'diode', 'svdiode', 'bjt', 'svbjt', 
               'mosEKV', 'mosACM', 'mosACMs', 'mosBSIM3v3',
               'mesfetc',
               'tlinpy4', 'tlinps4']

# Add here any modules to be imported in addition to netElemList
__all__ = netElemList + ['cppaddev']

devClass = {}
parType = {}

def add_device(device):
    """
    Appends one device model to the global dictionary (devClass)
    """
    devClass[device.devType] = device
    # If device has defined the makeAutoThermal flag then attempt
    # to create electrothermal device
    if device.makeAutoThermal:
        tDevice = autoThermal.thermal_device(device)
        devClass[tDevice.devType] = tDevice

for modname in netElemList:
    module = __import__(modname, globals(), locals())
    try:
        for device in module.devList:
            add_device(device)
    except AttributeError:
        # If there is no deviceList then there must be a Device class
        add_device(module.Device)


        
        




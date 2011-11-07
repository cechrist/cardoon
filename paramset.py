"""
:mod:`paramset` -- Classes for parameter handling
-------------------------------------------------

.. moduleauthor:: Carlos Christoffersen

Handles sets of parameters used for Elements, Models and Analyses
"""


class ParamError(Exception):
    """
    Used for exceptions raised in this module
    """
    pass


class ParamSet:
    """
    Handles a set of parameters with values. 

    Useful for devices, analyses and anything that accepts
    parameters. 
    """

    def __init__(self, paramDict):
        """
        Parameter definitions given in a dictionary with the following
        format::

            paramDict = dict(
                w = ('Channel width', 'm', float, 10e-6),
                l = ('Channel length', 'm', float, 10e-6),
                vt0 = ('Threshold Voltage', 'V', float, 0.532),
                vsat = ('Velocity Saturation', 'm/s', float, 8e4),
                tox = ('Oxide Thickness', 'm', float, 7.5e-9)
                )

        """
        # 
        self.paramDict = paramDict
        # valueList is a dictionary with overriding parameter values
        self.valueDict = dict()

    # Printing/formatting stuff -----------------------------------------

    def __str__(self):
        """
        Convert to string
        """
        desc = ' Name    |   Value    | Unit  \n'
        desc += '------------------------------\n'
        for key in sorted(self.paramDict.iterkeys()):
            paraminfo = self.paramDict[key]
            try:
                value = getattr(self, key)
            except AttributeError:
                try:
                    value = self.valueDict[key]
                except KeyError:
                    value = paraminfo[3]
            desc += '{0:^8} | {1:^10} | {2:^5}\n'.format(key, 
                                                         value, paraminfo[1])
        return(desc)

    def netlist_string(self):
        """
        Output netlist-formatted string listing only non-default parameters
        """
        desc = ''
        for param, value in self.valueDict.iteritems():
            desc += '{0} = {1} '.format(param, value)
        return(desc)

    def describe_parameters(self, paramName = None):
        """
        Returns a string with parameter information 

        If no parameter is specified all parameters are listed
        """
        helpstring = ' ========= ============ ============ ===================================================== \n'
        helpstring += ' Name       Default      Unit         Description                                          \n'
        helpstring += ' ========= ============ ============ ===================================================== \n'
        if paramName:
            helpstring = self.format(paramName)
        else:
            # Format the whole list of parameters. This could be
            # achieved using functional programming but no need for
            # now.
            for param in self.paramDict.iterkeys():
                helpstring += self.format(param)
        helpstring += ' ========= ============ ============ ===================================================== \n'
        return helpstring

    def format(self, paramName):
        """
        Returns a string describing parameter named param
        """
        try:
            paraminfo = self.paramDict[paramName]
        except NameError:
            out = 'Parameter not found.\n'
        else:
            out = ' {0:<8}   {1:<10}   {2:<10}   {3:<52} \n'.format(paramName,
                                 paraminfo[3], paraminfo[1], paraminfo[0])
        return out

    def list_parameters(self):
        """
        Briefly list parameter names
        """
        desc = ''
        for key in self.paramDict.iterkeys():
            desc += ' ' + key
        return desc

    # Individual parameter action/information ----------------------------

    def set_param(self, paramName, value):
        """ 
        Set parameter given with paramName to value

        Mostly used by the parser.  Note that actual attribute is not
        set until set_attributes() is executed.  A check is made to
        ensure the parameter name is valid
        """
        if self.paramDict.has_key(paramName):
            self.valueDict[paramName] = value
        else:
            raise ParamError('{0}: not a valid parameter name'.format(paramName))


    def get_type(self, paramName):
        """
        Returns type if parameter exists
        """
        try:
            return self.paramDict[paramName][2]
        except KeyError:
            raise ParamError('{0}: not a valid parameter name'.format(paramName))

    def is_set(self, paramName):
        """
        Returns True if paramName is valid and manually set
        """
        if self.valueDict.has_key(paramName):
            return True
        else:
            return False

    def get_float_attributes(self):
        """
        Returns a list with the names and values of 'float' attributes

        This is handy for sensitivity calculation. Each item in the
        list is a tuple: (<name>, <value>)

        The set is assumed to be already initialized. Values are taken
        directly from instance attributes.
        """
        paramList = list()
        for key, value in self.paramDict.iteritems():
            if value[2] == float:
                paramList.append((key, getattr(self, key)))
        return paramList
                

    # Conversion to attributes ----------------------------------------

    def clean_attributes(self):
        """
        Delete parameter attributes
        """
        for key, value in self.paramDict.iteritems():
            try:
                delattr(self, key)
            except AttributeError:
                pass

    def set_attributes(self, useDefaults = True):
        """
        Set attributes named after parameters in self.valueDict

        If useDefaults is True, set defaults from self.paramDict 
        """
        for key, value in self.valueDict.iteritems():
            setattr(self, key, value)

        if useDefaults:
            for key, value in self.paramDict.iteritems():
                # Only add these if not already set
                if not hasattr(self, key):
                    setattr(self, key, value[3])


                
#--------------------------------------------------------------------
class Model(ParamSet):
    """
    Provides '.model' functionality (used for Elements)

    Provides some extra functionality in addition to ParamSet methods
    """

    def __init__(self, name, modelType, paramDict):
        """
        name is the name of this particular set of parameters

        modelType is a string describing the type of model that uses
        this parameter set.

        Parameter definitions given in a dictionary
        """
        # 
        self.name = name
        self.modelType = modelType
        ParamSet.__init__(self, paramDict)

    def __str__(self):
        """
        Convert to string
        """
        desc = 'Model: {0}, Type: {1}\n\n'.format(self.name, self.modelType)
        desc += ParamSet.__str__(self)
        return(desc)

    def netlist_string(self):
        """
        Output netlist-formatted string using pset definitions and
        attributes values from var
        """
        desc = '.model {0} {1} ('.format(self.name, self.modelType)
        desc += ParamSet.netlist_string(self)
        desc += ')'
        return(desc)

    def set_missing_attributes(self, target):
        """
        Set attributes not present in target with values in stored in self
        """
        for key, paraminfo in self.paramDict.iteritems():
            # if attribute not in target, set it to value
            if not hasattr(target, key):
                # value not in target
                try:
                    # value in model
                    value = self.valueDict[key]
                except KeyError:
                    # Use default
                    value = paraminfo[3]
                setattr(target, key, value)

        
    

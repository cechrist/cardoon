"""
:mod:`circuit` -- Classes for internal circuit representation
-------------------------------------------------------------

.. moduleauthor:: Carlos Christoffersen

Example of how to use this module
+++++++++++++++++++++++++++++++++

The following example was originally taken from
``analyses/nodal.py``. The objective of this code is to assign
row/column numbers to all non-reference nodes and -1 to all reference
nodes in a ``Circuit`` instance (``ckt``)::

    # make a list of all non-reference terminals in circuit 
    ckt.nD_termList = ckt.termDict.values() + ckt.get_internal_terms()
    # get reference node (if any)
    if ckt.has_term('gnd'):
        ckt.nD_ref = ckt.get_term('gnd')
        # remove ground node from terminal list
        ckt.nD_termList.remove(ckt.nD_ref)
    # remove ground node from terminal list
    ckt.nD_termList.remove(ckt.nD_ref)
    # For reference nodes assign -1
    ckt.nD_ref.nD_namRC = -1
    # Make a list of all elements
    ckt.nD_elemList = ckt.elemDict.values()
    # Set RC number of reference terminals to -1
    for elem in ckt.nD_elemList:
        if elem.localReference:
            elem.connection[elem.localReference].nD_namRC = -1
    # Assign a number (starting from 0) to all nodes.
    for i, term in enumerate(ckt.nD_termList):
        term.nD_namRC = i

The following example shows how to create and add devices from the
point of view of a parser::

    import circuit as cir
    from devices import devClass
    
    # Create circuit:
    mainckt = cir.Circuit()
    
    # Add elements with private model
    
    # creates element type 'ind'
    dev = devClass['ind']('Lchoke')
    
    # override 'l' parameter value (otherwise default value is used)
    dev.set_param('l', 1e-9) 
    
    # add to circuit 
    mainckt.add_elem(dev)
    
    # connect to terminals 'vdd' and and 'drain', automatically created
    mainckt.connect(dev, ['vdd', 'drain'])
    
    # add element with shared (public) model: model 'my_nmos' is created
    # if it is not already in circuit
    m1 = devClass['mosacm']('mos1n')
    mainckt.add_elem(m1, 'my_nmos')
    
    # To later retrieve model
    model = cir.get_model('my_nmos')
    # if model not found returns None, use cir.add_model() if needed
    
    # Subcircuits: create and add subcircuit instance
    x1 = cir.Xsubckt('x1', 'LM741')
    mainckt.add_subckt(x1)
    # Connect subcircuit instance with cir.connect()
    
    # Subcircuit definition: tl2 are the external connections
    tl2 = ['vdd', 'vee', 'in', 'out']
    amp1 = cir.SubCircuit('LM741', tl2)
    # now you can treat amp1 as a regular circuit instance
    dev = devClass['ind']('L2')
    dev.override = [('l', 2e-9)] 
    amp1.add_elem(dev)
    amp1.connect(dev, ['in', 'n1'])


Notes
+++++

* If anything goes wrong a CircuitError exception is thrown

* The device type is automatically added to Element names. This
  prevents name conflicts from devices of different type. Example:
  'mosacm:mos1n'

* When all circuits and subcircuits are entered, you can flatten() the
  main circuit to remove subcircuit hierarchy (optional).

* Before elements can be used, the circuit should be initialized with
  init(). This performs some simple checks and initializes all
  elements in circuit (and sub-circuits if not flattened).

* .model statements have global scope: they can be used and referred
  from any circuit/subcircuit.

* The same applies to subcircuit definitions. They are global and can
  be referred from any circuit/subcircuit.

"""

from __future__ import print_function
from paramset import ParamSet, Model
from globalVars import glVar
import copy

#---------------------------------------------------------------------
class CircuitError(Exception):
    """
    Used for exceptions raised in this module
    """
    pass


#---------------------------------------------------------------------
class GraphNode(object):
    """
    Simple graph node base class (not circuit node) 

    Used for circuit Elements, subcircuits and terminals.  Links are
    stored using lists.
    """
    
    def __init__(self, instanceName):
        """
        instanceName is the Name of the node. No check is made here that
        the name is not already in use.
        """
        self.instanceName = instanceName
        self.connection = []

    def __str__(self):
        """convert to string"""
        desc = '"' + self.instanceName + '"\n' 
        desc += 'Linked nodes: '
        for n in self.connection:
            desc += ' ' + n.instanceName
        return(desc)

#---------------------------------------------------------------------
class Terminal(GraphNode):
    """
    Represent circuit terminals (i.e., nodes)

    This class should used only for 'external' terminals (i.e.,
    terminals that appear in the netlist). See also InternalTerminal
    """
    # By default terminals have volts as units. Devices may change this
    unit = 'V'

    def __init__(self, instanceName):
        # Call base class constructors
        GraphNode.__init__(self, instanceName)

    def __str__(self):
        """convert to string"""
        desc = 'Terminal: {1}\nUnit: {0}'.format(self.unit, 
                                                 GraphNode.__str__(self))
        return(desc)
    
    def get_label(self):
        """Return label for plots/tables in a formatted string"""
        return 'T: "{0}"'.format(self.instanceName)

#---------------------------------------------------------------------
class InternalTerminal(Terminal):
    """
    Represent terminals that are internal to one Element instance

    They only have one connection (the parent Element instance)
    """

    def __init__(self, element, name):
        """
        Creates and connects internal terminal

        element: parent Element instance
        name: internal name (only unique to parent Element)
        """
        # Call base class constructor
        Terminal.__init__(self, name)
        # Connect to parent element
        self.connection.append(element)
        element.connection.append(self)
        
    def __str__(self):
        """convert to string"""
        desc = 'Internal Terminal: "{1}:{2}", Unit: {0}'.format(
            self.unit,
            self.connection[0].instanceName,
            self.instanceName)
        return(desc)

    def get_label(self):
        """Return label for plots/tables in a formatted string"""
        return 'IT: "{0}:{1}"'.format(self.connection[0].instanceName,
                                      self.instanceName)

#---------------------------------------------------------------------
class Element(GraphNode, ParamSet):
    """
    Base class for circuit Elements. 

    Default flags are set here. 
    """
    # numTerms = 0 means do not (automatically) check number of connections
    numTerms = 0

    # temperature parameter: must be called 'temp' (if needed at all)
    tempItem = (('temp', ('Device temperature (None: use global temp.)', 
                          'C', float, None)), )
    
    # Has electro-thermal implementation (using autothermal.py)
    # Element must implement all nonlinear functions to use this
    makeAutoThermal = False

    # Flags for nonlinear devices. 
    isNonlinear = False
    # Number of time-delayed control ports
    nDelays = 0
    # Used for transmission lines, etc.
    isFreqDefined = False
    # Element contributes to source vector (DC, time-domain and freq-domain)
    isDCSource = False
    isTDSource = False
    isFDSource = False

    # General attributes
    linearVCCS = []
    linearVCQS = []
    noisePorts = ()

    # Local reference: if not zero, reference (internal) terminal number
    localReference = 0

    def __init__(self, instanceName):
        """
        The name of the element is formed by combining the given
        instanceName and the device type. 

        Example: diode:d1
        """
        # Call base class constructors
        GraphNode.__init__(self, self.devType + ':' + instanceName)
        # Note: paramDict must be defined by the derived class
        ParamSet.__init__(self, self.paramDict)
        # Default is not to have a separate model
        self.dotModel = False

    # Printing and info-related functions ----------------------------------
    def __str__(self):
        """convert to string"""
        desc = 'Element: ' + GraphNode.__str__(self)
        desc += '\nDevice type: ' + self.devType + '\n' 
        if self.dotModel:
            desc += 'Model: {0}\n'.format(self.dotModel.name)
        desc += 'Overridden parameters: ' + ParamSet.netlist_string(self)
        return(desc)

    def netlist_string(self):
        """
        Output netlist-formatted string in netlist format
        """
        desc = '{0} '.format(self.instanceName)
        # Add terminals
        for i, term in enumerate(self.connection):
            # Do not include internal terminals. The following works
            # even when numTerms is not set.
            if issubclass(type(term), InternalTerminal):
                break
            desc += term.instanceName + ' '
        # Model (if any)
        if self.dotModel:
            desc += 'model = {0} '.format(self.dotModel.name)
        # Parameters
        desc += ParamSet.netlist_string(self)
        return(desc)

    def print_vars(self):
        """
        Nicely print parameter list and OP (if any) with values and units
        """
        print('Parameter values:\n')
        print(ParamSet.__str__(self))

    # General initialization --------------------------------------------
    def init(self):
        """
        General initialization function

        Set the attribute values, check basic terminal connectivity
        and process parameters.
        """
        self.set_attributes()
        self.check_terms()
        self.process_params()

    # Parameter-related functions ----------------------------------------
    def is_set(self, paramName):
        """
        Returns True if paramName is valid and manually set
        """
        # check if parameter in model first
        if self.dotModel:
            answer = self.dotModel.is_set(paramName)
        else:
            answer = False
        # now check in element
        if not answer:
            answer = ParamSet.is_set(self, paramName)
        return answer

    def set_attributes(self):
        """
        Set parameters as attributes. 

        Priority is as follows: first manually set parameters, then
        manually set parameters in model (if any) and finally default
        values
        """
        # Ambient temperature (temp) initially set to global
        # temperature, but may be overriden by the model or the
        # element line attributes
        self.temp = glVar.temp
        # Set overrides first
        if self.dotModel:
            ParamSet.set_attributes(self, useDefaults = False)
            # Get other attributes from model
            self.dotModel.set_missing_attributes(self)
        else:
            ParamSet.set_attributes(self, useDefaults = True)

    # Connectivity-related functions -----------------------------------------
    def check_terms(self):
        """
        Checks terminal connections

        If numTerms is not zero, checks that the number of connected
        terminals is equal to numTerms. Raises an exception if they do
        not match. 

        Else sets numTerms = number of connected terminals
        """
        if self.numTerms:
            if (len(self.connection) != self.numTerms):
                raise CircuitError(self.instanceName + 
                                   ': must have ' + str(self.numTerms) 
                                   + ' terminals.')
        else:
            # Set numterms to number of external connections (useful
            # to quickly find internal terminals)
            self.numTerms = len(self.connection)

    def disconnect(self, terminal):
        """
        Disconnect a terminal from an element. Arguments are Element and
        Terminal instances. 
    
        Assumption is that if terminal is in element's list, then the nodes are
        linked and element should be on terminal's list. If nodes not connected
        an exception is raised. 
        """
        try:
            self.connection.remove(terminal)
            terminal.connection.remove(self)
        except ValueError:
            raise CircuitError('Nodes not linked')

    def add_internal_term(self, name, unit):
        """
        Create and connect one internal terminal

        name: internal terminal name
        unit: internal variable unit

        Returns internal terminal index
        """
        # Create internal term (connects automatically)
        term = InternalTerminal(self, name)
        term.unit = unit
        return len(self.connection) - 1

    def add_reference_term(self):
        """
        Create and connect one internal terminal to be used as local reference
        """
        # There is no need for more than one
        assert not self.localReference
        # Create internal term (connects automatically)
        term = InternalTerminal(self, 'lref')
        term.unit = '-'
        # Set to reference terminal number
        self.localReference = len(self.connection) - 1
        return self.localReference

    def get_internal_terms(self):
        """
        Returns a list of internal terms (if any) excluding local references
        """
        intTerms = self.connection[self.numTerms:]
        if self.localReference:
            intTerms.pop(self.localReference - self.numTerms)
        return intTerms
    
    def clean_internal_terms(self):
        """
        Disconnect any internal terms, 

        Normally used before calling process_params() for a second
        time or when an element is removed from circuit.
        """
        for term in self.connection[self.numTerms:]:
            self.disconnect(term)
        # Chop adjacency list
        self.connection = self.connection[:self.numTerms]
        # Clean local reference
        self.localReference = 0

#---------------------------------------------------------------------
class Xsubckt(GraphNode):
    """
    Represent subcircuit instances (not definitions, use SubCircuit
    for those)
    """

    def __init__(self, instanceName, cktName):
        """
        instanceName: self-documenting

        cktName: name of the circuit that contains the definition for
        this instance.
        """
        # Call base class constructors
        assert instanceName[0] == 'x'
        GraphNode.__init__(self, instanceName)
        self.cktName = cktName
        

    def __str__(self):
        """convert to string"""
        desc = 'xsubckt ' + GraphNode.__str__(self)
        desc += '\nSubcircuit definition: ' + self.cktName
        return(desc)


    def netlist_string(self):
        """ 
        Convert to string in netlist format
        """
        desc = self.instanceName + ' '
        for term in self.connection:
            desc += term.instanceName + ' '
        desc += self.cktName 
        return desc

    def check_terms(self):
        # First make sure the definition exists
        try:
            cktDef = Circuit.cktDict[self.cktName]
        except KeyError:
            raise CircuitError(self.instanceName + \
                                ': subcircuit definition "'\
                                + self.cktName + '" not found')
        else:
            if  len(self.connection) != len(cktDef.extConnectionList):
                raise CircuitError(self.instanceName + \
                                    ': xsubckt connections do not match '\
                                    + 'definition in "' + self.cktName \
                                    + '"')


#---------------------------------------------------------------------
class Circuit(object):
    """
    Holds a circuit. 

    There are 2 global dictionaries defined at the class level:

    * cktDict: Contains references to all circuit/subcircuit
      definitions

    * modelDict: References to all models in any circuit. Thus .model
      statements are global and can be defined and referred anywhere

    Element, Xsubckt and (external) Terminal references are stored in
    dictionaries, empty by default:

    * elemDict
    * subcktDict
    * termDict

    Internal terminals must be accessed directly from the parent
    Element instance.

    Ground node: terminals '0' and 'gnd' are considered to be the same
    reference node. If a circuit does not contain a ground node then
    it is up to the user to set a reference.

    Plot/Save requests are stored in lists: plotReqList/saveReqList

    In the future we may implement topology checking utilities here.
    """

    cktDict = dict()
    modelDict = dict()

    def __init__(self, name):
        """
        name: circuit name.  
        """
        self.name = name
        self._initialized = False
        self._flattened = False
        # adds itself to main dictionary
        if self.cktDict.has_key(name):
            raise CircuitError('Circuit "' + name + '" already exists')
        self.cktDict[name] = self
        # Create (empty) dictionaries
        self.termDict = dict()
        self.elemDict = dict()
        self.subcktDict = dict()
        self.plotReqList = list()
        self.saveReqList = list()

    # Printing stuff ---------------------------------------------

    def __str__(self):
        """
        Generates a short descriptive string
        """
        desc = 'Circuit: {0}\n'.format(self.name)
        if hasattr(self, 'title'):
            desc = self.title + '\n'
        return desc

    def netlist_string(self):
        """
        Generates a 'netlist-like' description of the circuit
        """
        desc = ''
        if hasattr(self, 'title'):
            desc = self.title + '\n'
        desc += '# Circuit: ' + self.name + '\n'

        desc += '#                         *** Elements ***\n'
        for elem in self.elemDict.itervalues():
            desc += elem.netlist_string() + '\n'
            
        usedSubCKTs = set()
        if not self._flattened and self.subcktDict:
            desc += '#                     *** Subcircuit instances ***\n'
            for subckt in self.subcktDict.itervalues():
                desc += subckt.netlist_string() + '\n'
                # Save used subckt name
                usedSubCKTs.add(subckt.cktName)
            desc += '#                     *** Subcircuit Definitions ***\n'
            for cktName in usedSubCKTs:
                desc += Circuit.cktDict[cktName].netlist_string()
        
        for outreq in self.plotReqList:
            desc += '\n.plot ' + outreq.type 
            for name in outreq.varlist:
                desc += ' ' + name
        for outreq in self.saveReqList:
            desc += '\n.save ' + outreq.type 
            for name in outreq.varlist:
                desc += ' ' + name
        desc += '\n'

        return desc

    def globals_to_str(self):
        """
        Convert all global stuff to netlist format

        This includes: models, netlist variables and .options variables
        """
        desc = '.vars '
        for param, value in ParamSet.netVar.iteritems():
            desc += '{0} = {1} '.format(param, value)
        desc += '\n\n'
        desc += '.options ' + glVar.netlist_string()
        desc += '\n\n'
        desc += '#                  *** Models *** \n'
        for model in Circuit.modelDict.itervalues():
            desc += model.netlist_string() + '\n\n'
        return desc

    # Actions on the whole circuit --------------------------------------

    def copy(self, newName):
        """ 
        Returns a copy of an uninitialized circuit

        newName is the circuit name for the copy

        The elements in the copy point to the ParamSet of the original
        element (models are not cloned). Elements are shallow-copied,
        but terminal connections are re-generated.

        Example::

          ckt2 = ckt1.copy('ckt1copy')

        If the circuit has been flattened, the copy does not include
        subcircuit instances, otherwise a shallow copy of subcircuit
        instances are generated.
        """
        assert not self._initialized 
        # create new circuit with given name
        cktCopy = Circuit(newName)
        # Get all elements, add and connect them to this circuit
        for elem in self.elemDict.itervalues():
            # We need a copy of the element without the terminal
            # connections. Make shallow copy:
            elemCopy = copy.copy(elem)
            # Clean terminal list
            elemCopy.connection = []
            # Add to circuit 
            cktCopy.add_elem(elemCopy)
            # Create connection list
            termList = [term.instanceName for term in elem.connection]
            cktCopy.connect(elemCopy, termList)
        if not self._flattened:
            for xsubckt in self.subcktDict.itervalues():
                xsubcktCopy = copy.copy(xsubckt)
                # Clean terminal list
                xsubcktCopy.connection = []
                # Add to circuit 
                cktCopy.add_subckt(xsubcktCopy)
                # Create connection list
                termList = [term.instanceName for term in xsubckt.connection]
                cktCopy.connect(xsubcktCopy, termList)
        return cktCopy

    def init(self):
        """
        To be used after all elements/terminals have been created. Can
        be called multiple times but only has an effect if
        ``self._initialized`` is False

        1. Initialize all elements. This includes a check for terminal
           connections and parameter processing. Elements may add
           internal terminals at this point.

        2. If not flattened, checks that subcircuit definitions exist
           and checks number of terminal connections.

        3. Checks that output requests contain valid terminal names

        """
        # Can only initialize once
        if self._initialized:
            return
        for elem in self.elemDict.itervalues():
            elem.init()

        if not self._flattened:
            for xsubckt in self.subcktDict.itervalues():
                xsubckt.check_terms()
                # Initialize subcircuit definition
                Circuit.cktDict[xsubckt.cktName].init()

        # Check and initialize output requests
        for outreq in self.plotReqList + self.saveReqList:
            outreq.termlist = []
            for termname in outreq.varlist:
                # Distinguish external/internal terminals
                if type(termname) == str:
                    # External terminal
                    if termname == '0':
                        # Special treatment for ground terminal
                        termname = 'gnd'
                    try:
                        outreq.termlist.append(self.termDict[termname])
                    except KeyError:
                        raise CircuitError(
                            'Output request for nonexistent terminal: ' 
                            + termname)
                else:
                    # Internal terminal
                    try:
                        elem = self.elemDict[termname[0]]
                    except KeyError:
                        raise CircuitError(
                            'Output request for nonexistent element: ' 
                            + termname[0] + ':' + termname[1])
                    iterm = None
                    for term in elem.get_internal_terms():
                        if term.instanceName == termname[1]:
                            iterm = term
                            break
                    if iterm:
                        outreq.termlist.append(iterm)
                    else:
                        raise CircuitError(
                            'Output request for nonexistent terminal: '
                            + termname[0] + ':' + termname[1])


        self._initialized = True


    def get_internal_terms(self):
        """
        Returns a list with all non-reference internal terminals

        Circuit must be initialized first. Note that the same effect
        can be achieved by directly polling the elements.
        """
        assert self._initialized
        intTermList = []
        for elem in self.elemDict.itervalues():
            intTermList += elem.get_internal_terms()
        return intTermList


    def flatten(self):
        """ 
        Expand all subcircuit instances in place, assuming the circuit
        and subcircuits have not been initialized. 

        The flattened elements point to the ParamSet of the original
        element (models are not cloned). Elements are shallow-copied,
        but terminal connections are re-generated.
        """
        assert not self._initialized 
        # Can only flatten once
        if self._flattened:
            return
        for xsubckt in self.subcktDict.itervalues():
            cktDef = self.cktDict[xsubckt.cktName]
            # Recursively flatten nested subcircuits
            if not cktDef._flattened:
                cktDef.flatten()
            # Get circuit from subcircuit definition
            cktDef.get_circuit(xsubckt, self)
        self._flattened = True

    def check_sanity(self):
        """
        This is slow for large circuits.  Use for debugging to
        make sure that all elements and terminals have been added to
        the circuit. Also checks for floating terminals.

        For now only works for unflattened circuits
        """
        for elem in self.elemDict.itervalues():
            for term in elem.connection:
                assert self.termDict[term.instanceName] == term

        if not self._flattened:
            # If flattened Xsubckt instances are ignored
            for subckt in self.subcktDict.itervalues():
                for term in subckt.connection:
                    assert self.termDict[term.instanceName] == term

        for term in self.termDict.itervalues():
            # Make sure there are no floating terminals
            assert term.connection
            for node in term.connection:
                if isinstance(node, Element):
                    assert self.elemDict[node.instanceName] == node
                else:
                    assert self.subcktDict[node.instanceName] == node


    # Actions on individual elements/terminals -------------------------
            
    def connect(self, element, termList):
        """
        Connect an Element (or subckt) instance to terminals 

        Terminals are specified by a list of terminal names
        (termList). If a terminal does not exist in the circuit it is
        created and added.

        **Important notes**: 

            * Element/subckt should be in circuit dictionary. No
              attempt to check this is made (use ``add_elem()``).
        
            * Order in the adjacency list is important for Elements
        
              For example, for a MOSFET the first node corresponds to
              the drain, the second to the gate, etc.

        """
        # Some sanity checking. This function should be used once per
        # element (we could change this)
        assert not element.connection

        for termName in termList:
            terminal = self.get_term(termName)
            terminal.connection.append(element)
            element.connection.append(terminal)
    

    def add_elem(self, elem, modelName = None):
        """ 
        Adds an element to a circuit.
    
        If modelName is not given it is assumed that no global model
        will be used
    
        Otherwise the model is retrieved (or created if necessary) and
        assigned to elem.dotModel. This model can be shared with other
        elements.
    
        A check is made to make sure the instance name is unique
        """
        if self.elemDict.has_key(elem.instanceName):
            raise CircuitError(elem.instanceName + ': Element already exists')
        else:
            self.elemDict[elem.instanceName] = elem
    
        if modelName:
            if elem.dotModel:
                raise CircuitError('{0} already has an ' +
                                   'assigned model'.format(elem.instanceName))
            # Public model: Check if model already in circuit
            if Circuit.modelDict.has_key(modelName):
                model = Circuit.modelDict[modelName]
                if model.modelType != elem.devType:
                    raise CircuitError(
                        'Incorrect model type "{0}" for element "{1}"'.format(
                            model.modelType, elem.instanceName))
                elem.dotModel = model
            else:
                elem.dotModel = Model(modelName, elem.devType, elem.paramDict)
                Circuit.modelDict[modelName] = elem.dotModel
    
    def remove_elem(self, elemName):
        """
        Disconnect and remove an element from a circuit. Internal
        terminals are also removed.
        """
        # Make sure the node is in circuit (and take it out)
        try:
            elem = self.elemDict.pop(elemName)
        except KeyError:
            raise CircuitError(elemName + ': Element not found' )
        else:
            # Remove internal terminals first
            elem.clean_internal_terms(self)
            # unlink connections (may be a bit slow)
            for n1 in elem.connection:
                # Disconnect any terminals from elem
                n1.connection.remove(elem)
                    
    def add_subckt(self, xsubckt):
        """
        Adds a Xsubckt instance to circuit
        """
        if self.subcktDict.has_key(xsubckt.instanceName):
            raise CircuitError('add_subckt: Subcircuit instance name "'\
                                + xsubckt.instanceName + '" already in use')
        else:
            self.subcktDict[xsubckt.instanceName] = xsubckt 
    
    def get_term(self, termName):
        """ 
        Returns an external terminal instance with the given name. 

        A new instance is created if necessary
        """
        if termName == '0':
            # Special treatment for ground terminal
            termName = 'gnd'
        if self.termDict.has_key(termName):
            term = self.termDict[termName]
        else:
            term = self.termDict[termName] = Terminal(termName)
    
        return term

    def has_term(self, termName):
        """ 
        Returns True if terminal present in circuit
        """
        if termName == '0':
            # Special treatment for ground terminal
            termName = 'gnd'
        if self.termDict.has_key(termName):
            return True
        else:
            return False
            
    def find_term(self, termName):
        """
        Find a terminal (external/internal) by name

        termName: full name, including container Element name for
                  internal terminals

        Returns terminal instance if found or raises CircuitError
        exception if not
        """
        result = None
        # May be it is an internal terminal. Search for
        # containing element
        token = termName.rpartition(':')
        if token[0]:
            # looks like an internal terminal
            t1 = token[0].rpartition(':')
            if t1[0]:
                try:
                    terms = self.elemDict[token[0]].get_internal_terms()
                    for t in terms:
                        if t.instanceName == token[2]:
                            result = t
                            break
                    if result == None:
                        raise CircuitError('Terminal "' + token[2] + 
                                           '" not found in element "' + 
                                           token[0] + '"')
                except KeyError:
                    raise CircuitError('Element "' + token[0] + '" not found')
            else:
                 raise CircuitError('Terminal name syntax is wrong: "' 
                                    + termName + '"')
        else:
            # It can be an external terminal
            try:
                result = self.termDict[termName]
            except KeyError:
                raise CircuitError('Terminal "' + termName + '" not found')

        return result


    def remove_term(self, termName):
        """
        Removes a terminal from a circuit. Links are *not* removed, you
        must take care of that.
        """
        # Make sure the node is in circuit (and take it out)
        try:
            term = self.termDict.pop(termName)
        except KeyError:
            raise CircuitError(termName + ': Terminal not found' )
    
    def add_plot_request(self, outreq):
        self.plotReqList.append(outreq)

    def add_save_request(self, outreq):
        self.saveReqList.append(outreq)

    def get_requested_terms(self, reqtype):
        """
        Returns a set with terminals to be plotted or saved 

        Set generated from all plot and save requests of the specified
        type.
        """
        # Uses a set to avoid repeated terminals
        termSet = set()
        for outreq in self.plotReqList + self.saveReqList:
            if outreq.type == reqtype:
                for term in outreq.termlist:
                    termSet.add(term)
        return termSet


#---------------------------------------------------------------------
class SubCircuit(Circuit):
    """
    Almost identical to Circuit but adds support for external connections.

    External subcircuit terminals have an additional attribute set:
    ``subCKTconnection``. This attribute is set to the subcircuit's
    connection number. Example::

        .subckt  inverter in out vplus vminus
        
        in.subCKTconnection = 0
        out.subCKTconnection = 1
        vplus.subCKTconnection = 2
        vminus.subCKTconnection = 3

    The reference node ('gnd' or '0') is global , *i.e.*, a terminal
    named 'gnd' in a subcircuit is assumed to be connected to the same
    ground node in all other circuits/subcircuits. To avoid problems,
    the ground terminal can not be included in the external terminal
    list. One of the possible problems is short-circuiting a node to
    ground when a subcircuit is connected. Example of invalid code::

        vdc:vcc 1 0 vdc=10V
        vdc:vee 2 0 vdc=-10V
        xamp1 1 2 in out amplifier
        
        .subckt amplifier 1 gnd in out
        res:rc 1 out r=5k
        cap:cin in 2 c=1uF
        bjt:q1 out 2 gnd type=npn
        .ends

    Possible solutions: 

      1. Rename 'gnd' node in subcircuit to something else 

      2. Remove 'gnd' node (implicitly connected to terminal '0' in
         main circuit) from subcircuit external terminal list

    """
    def __init__(self, name, termList):
        """
        termList is a list of terminal names that connect the
        subcircuit to the external world. These are added to the
        subcircuit instance.
        """
        Circuit.__init__(self, name)
        self.extConnectionList = termList
        for i, termName in enumerate(termList):
            terminal = self.get_term(termName)
            if terminal.instanceName == 'gnd':
                raise CircuitError(name +
""": gnd is the global reference. It is implicitly connected and can 
not be included in the external terminal list in subckt""")
            # Save terminal connection number here
            terminal.subCKTconnection = i

    def copy(self, newName):
        """ 
        Returns a copy of an uninitialized circuit

        newName is the circuit name for the copy

        Similar as Circuit.copy() but also takes care of subcircuit
        connections
        """
        assert not self._initialized 
        # create new circuit with the given name
        cktCopy = SubCircuit(newName, self.extConnectionList)
        # Get all elements, add and connect them to this circuit
        for elem in self.elemDict.itervalues():
            # We need a copy of the element without the terminal
            # connections. Make shallow copy:
            elemCopy = copy.copy(elem)
            # Clean terminal list
            elemCopy.connection = []
            # Add to circuit 
            cktCopy.add_elem(elemCopy)
            # Create connection list
            termList = [term.instanceName for term in elem.connection]
            cktCopy.connect(elemCopy, termList)
        if not self._flattened:
            for xsubckt in self.subcktDict.itervalues():
                xsubcktCopy = copy.copy(xsubckt)
                # Clean terminal list
                xsubcktCopy.connection = []
                # Add to circuit 
                cktCopy.add_subckt(xsubcktCopy)
                # Create connection list
                termList = [term.instanceName for term in xsubckt.connection]
                cktCopy.connect(xsubcktCopy, termList)


        return cktCopy

    def get_connections(self):
        """
        Returns a list of subcircuit terminals
        """
        return [self.get_term(termName) 
                for termName in self.extConnectionList]

    def get_circuit(self, xsubckt, target):
        """
        Dump copy of circuit into target (for hierarchy flattening)

        Subcircuit instances are ignored
        """
        # This code (originally in Circuit class) belongs here because
        # it uses attributes known only by the SubCircuit class.
        assert self._flattened
        # Get all elements, add and connect them to this circuit
        for elem in self.elemDict.itervalues():
            # We need a copy of the element without the terminal
            # connections. Make shallow copy:
            elemCopy = copy.copy(elem)
            # Clean terminal list
            elemCopy.connection = []
            # Change instance name
            elemCopy.instanceName = xsubckt.instanceName + \
                ':' + elem.instanceName
            # Add to circuit 
            target.add_elem(elemCopy)
            # Create connection list
            termList = []
            for term in elem.connection:
                if hasattr(term, 'subCKTconnection'):
                    # Must get terminal name from xsubckt
                    termName = \
                        xsubckt.connection[term.subCKTconnection].instanceName
                elif term.instanceName == 'gnd':
                    # Special treatment for global reference node
                    termName = 'gnd'
                else:
                    termName = xsubckt.instanceName + ':' + term.instanceName
                termList.append(termName)
            target.connect(elemCopy, termList)

    def netlist_string(self):
        """
        Ganerates a 'netlist-like' description of subcircuit.

        Produces the same output as the Circuit class plus extra
        .subckt and .ends lines.
        """
        desc = '#-----------------------------------------------------\n'
        desc += '.subckt ' + self.name 
        for termName in self.extConnectionList:
            desc += ' ' + termName
        desc += '\n'
        desc += Circuit.netlist_string(self)
        desc += '.ends\n\n'
        return desc


#---------------------------------------------------------------------
# Utility functions
#
# To avoid the overhead of using virtual functions, utility functions
# are defined outside the circuit class.
#
#---------------------------------------------------------------------
def get_mainckt():
    """
    Used to get the circuit called 'main' weather it has been already
    created or not
    """
    if not Circuit.cktDict.has_key('main'):
        return Circuit('main')
    else:
        return Circuit.cktDict['main']

def reset_allckt():
    """
    Erases all existing circuits

    Used mostly for debugging purposes
    """
    Circuit.cktDict = dict()
    Circuit.modelDict = dict()

#---------------------------------------------------------------------
def add_model(model):
    """ 
    Adds a public model (parameter set) to a circuit.

    A check is made to make sure the instance name is unique
    """
    if Circuit.modelDict.has_key(model.name):
        raise CircuitError(model.name + ': Model already exists')
    else:
        Circuit.modelDict[model.name] = model

#---------------------------------------------------------------------
def get_model(modelName):
    """
    Returns model reference if present in circuit, otherwise returns
    None.
    """
    if Circuit.modelDict.has_key(modelName):
        return  Circuit.modelDict[modelName]
    else:
        return None




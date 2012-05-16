"""
:mod:`circuit` -- Classes for internal circuit representation
-------------------------------------------------------------

.. moduleauthor:: Carlos Christoffersen

Example of how to use this module
+++++++++++++++++++++++++++++++++

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
    tl2 = ['vdd', 'gnd', 'in', 'out']
    amp1 = cir.SubCircuit('LM741', tl2)
    # now you can treat amp1 as a regular circuit instance
    dev = devClass['ind']('L2')
    dev.override = [('l', 2e-9)] 
    amp1.add_elem(dev)
    amp1.connect(dev, ['in', 'n1'])

The following example taken from ``analyses/nodal.py``. Here ``ckt``
is a reference to a ``Circuit`` object::

    # get ground node
    ckt.nD_ref = ckt.get_term(reference)

    # make a list of all non-reference terminals in circuit 
    ckt.nD_termList = ckt.termDict.values() + ckt.get_internal_terms()
    # remove ground node from terminal list
    ckt.nD_termList.remove(ckt.nD_ref)
    # Assign a number (0-inf) to all nodes. For reference nodes
    # assign -1 
    ckt.nD_ref.nD_namRC = -1
    # Make a list of all elements
    ckt.nD_elemList = ckt.elemDict.values()
    # Set RC number of reference terminals to -1
    for elem in ckt.nD_elemList:
        if elem.localReference:
            elem.neighbour[elem.localReference].nD_namRC = -1
    # For the future: use graph techniques to find the optimum
    # terminal order
    for i, term in enumerate(ckt.nD_termList):
        term.nD_namRC = i

    # Store internal RC numbers for later use
    for elem in ckt.nD_elemList:
        elem.nD_intRC = [term.nD_namRC for term in 
                         elem.neighbour[elem.numTerms:]]

    # Dimension is the number of unknowns to solve for
    ckt.nD_dimension = len(ckt.nD_termList)
    # Number of external terminals excluding reference
    ckt.nd_nterms = len(ckt.termDict.values()) - 1

    # Create specialized element lists
    ckt.nD_nlinElem = filter(lambda x: x.isNonlinear, ckt.nD_elemList)
    ckt.nD_freqDefinedElem = filter(lambda x: x.isFreqDefined, ckt.nD_elemList)
    ckt.nD_sourceDCElem = filter(lambda x: x.isDCSource, ckt.nD_elemList)
    ckt.nD_sourceTDElem = filter(lambda x: x.isTDSource, ckt.nD_elemList)
    ckt.nD_sourceFDElem = filter(lambda x: x.isFDSource, ckt.nD_elemList)


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
class GraphNode:
    """
    Simple graph node base class (not circuit node) 

    Used for circuit Elements, subcircuits and terminals.  Links are
    stored using lists.
    """
    
    def __init__(self, nodeName):
        """
        nodeName is the Name of the node. No check is made here that
        the name is not already in use.
        """
        self.nodeName = nodeName
        self.neighbour = []

    def __str__(self):
        """convert to string"""
        desc = self.nodeName + '\n' 
        desc += 'Linked nodes: '
        for n in self.neighbour:
            desc += ' ' + n.nodeName
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
        desc = 'Terminal({0}): {1}'.format(self.unit, GraphNode.__str__(self))
        return(desc)

#---------------------------------------------------------------------
class InternalTerminal(Terminal):
    """
    Represent terminals that are internal to one Element instance

    They only have one neighbour (the parent Element instance)
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
        self.neighbour.append(element)
        element.neighbour.append(self)
        

    def __str__(self):
        """convert to string"""
        desc = 'Internal Terminal ({0}): {1}:{2}'.format(
            self.unit,
            self.neighbour[0].nodeName,
            self.nodeName)
        return(desc)


#---------------------------------------------------------------------
class Element(GraphNode, ParamSet):
    """
    Base class for circuit Elements. 

    Default flags are set here. 
    """
    # numTerms = 0 means do not (automatically) check number of connections
    numTerms = 0

    # temperature parameter: must be called 'temp' (if needed at all)
    tempItem = (('temp', ('Device temperature', 'C', float, None)), )
    
    # Has electro-thermal implementation (using autothermal.py)
    # Element must implement all nonlinear functions to use this
    makeAutoThermal = False

    # Flags for nonlinear devices. 
    isNonlinear = False
    needsDelays = False
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
        desc = 'Element ' + GraphNode.__str__(self)
        desc += '\nDevice type: ' + self.devType + '\n' 
        if self.dotModel:
            desc += 'Model: {0}\n'.format(self.dotModel.name)
        desc += 'Overridden parameters: ' + ParamSet.netlist_string(self)
        return(desc)

    def netlist_string(self):
        """
        Output netlist-formatted string in netlist format
        """
        desc = '{0} '.format(self.nodeName)
        # Add terminals
        for i, term in enumerate(self.neighbour):
            # Do not include internal terminals. The following works
            # even when numTerms is not set.
            if issubclass(type(term), InternalTerminal):
                break
            desc += term.nodeName + ' '
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
        if hasattr(self, 'OP'):
            print('Operating point information:\n')
            print(' Variable  |  Value ')
            print('-------------------------')
            print(self.format_OP())

    def format_OP(self):
        """ 
        Return OP information in a formatted string
        """
        try:
            # Format operating point information
            s = ''
            for key in sorted(self.OP.iterkeys()):
                s += '{0:10} | {1}\n'.format(key, self.OP[key])
            return s
        except AttributeError:
            return ''

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
            # Get othe attributes from model
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
            if (len(self.neighbour) != self.numTerms):
                raise CircuitError(self.nodeName + 
                                   ': must have ' + str(self.numTerms) 
                                   + ' terminals.')
        else:
            # Set numterms to number of external connections (useful
            # to quickly find internal terminals)
            self.numTerms = len(self.neighbour)

    def disconnect(self, terminal):
        """
        Disconnect a terminal from an element. Arguments are Element and
        Terminal instances. 
    
        Assumption is that if terminal is in element's list, then the nodes are
        linked and element should be on terminal's list. If nodes not connected
        an exception is raised. 
        """
        try:
            self.neighbour.remove(terminal)
            terminal.neighbour.remove(self)
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
        return len(self.neighbour) - 1

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
        self.localReference = len(self.neighbour) - 1
        return self.localReference

    def get_internal_terms(self):
        """
        Returns a list of internal terms (if any) excluding local references
        """
        intTerms = self.neighbour[self.numTerms:]
        if self.localReference:
            intTerms.pop(self.localReference - self.numTerms)
        return intTerms
    
    def clean_internal_terms(self):
        """
        Disconnect any internal terms, 

        Normally used before calling process_params() for a second
        time or when an element is removed from circuit.
        """
        for term in self.neighbour[self.numTerms:]:
            self.disconnect(term)
        # Chop adjacency list
        self.neighbour = self.neighbour[:self.numTerms]
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
        desc = self.nodeName + ' '
        for term in self.neighbour:
            desc += term.nodeName + ' '
        desc += self.cktName 
        return desc

    def check_terms(self):
        # First make sure the definition exists
        try:
            cktDef = Circuit.cktDict[self.cktName]
        except KeyError:
            raise CircuitError(self.nodeName + \
                                ': subcircuit definition "'\
                                + self.cktName + '" not found')
        else:
            if  len(self.neighbour) != len(cktDef.extConnectionList):
                raise CircuitError(self.nodeName + \
                                    ': xsubckt connections do not match '\
                                    + 'definition in "' + self.cktName \
                                    + '"')


#---------------------------------------------------------------------
class OutRequest:
    """
    Holds a set of output variables requests

    Output request consist in: 

      1. Type of of request (``type``): dc, ac_*, tran, etc.

      2. List of variables (``varlist``): these are terminal names.

    After initialization the circuit adds a list of terminals in the
    ``termlist`` attribute.
    """
    validTypes = ['dc', 'ac_mag', 'ac_phase', 'ac_dB', 'tran', 'hb']

    def __init__(self, reqtype, varlist):
        if reqtype not in self.validTypes:
            raise CircuitError(
                'Not a valid output request type: {0}'.format(reqtype)
                + '\nValid types: {0}'.format(self.validTypes))
        self.type = reqtype
        self.varlist = varlist


#---------------------------------------------------------------------
class Circuit:
    """
    Holds a circuit. 

    There are 2 global dictionaries defined at the class level:

    * cktDict: Contains references to all circuit/subcircuit instances

    * modelDict: References to all models in any circuit. Thus .model
      statements are global and can be defined and referred anywhere

    Element, Xsubckt and (external) Terminal references are stored in
    dictionaries, empty by default:

    * elemDict
    * subcktDict
    * termDict

    Internal terminals must be accessed directly from the parent
    Element instance.

    Ground node: treated specially. Nodes '0' and 'gnd' are considered
    to be the same. If a circuit does not contain a ground node then
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
                if termname == '0':
                    # Special treatment for ground terminal
                    termname = 'gnd'
                try:
                    outreq.termlist.append(self.termDict[termname])
                except KeyError:
                    raise CircuitError(
                        'Output request for nonexistent terminal: ' + termname)

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
            # Get all elements, add and connect them to this circuit
            for elem in cktDef.elemDict.itervalues():
                # We need a copy of the element without the terminal
                # connections. Make shallow copy:
                newelem = copy.copy(elem)
                # Clean terminal list
                newelem.neighbour = []
                # Change instance name
                newelem.nodeName = xsubckt.nodeName + ':' + elem.nodeName
                # Add to circuit 
                self.add_elem(newelem)
                # Create connection list
                termList = []
                for term in elem.neighbour:
                    if hasattr(term, 'subCKTconnection'):
                        # Must get terminal name from xsubckt
                        termName = \
                            xsubckt.neighbour[term.subCKTconnection].nodeName
                    else:
                        termName = xsubckt.nodeName + ':' + term.nodeName
                    termList.append(termName)
                self.connect(newelem, termList)

        self._flattened = True

    def check_sanity(self):
        """
        This is slow for large circuits.  Use for debugging to
        make sure that all elements and terminals have been added to
        the circuit. Also checks for floating terminals.

        For now only works for unflattened circuits
        """
        for elem in self.elemDict.itervalues():
            for term in elem.neighbour:
                assert self.termDict[term.nodeName] == term

        if not self._flattened:
            # If flattened Xsubckt instances are ignored
            for subckt in self.subcktDict.itervalues():
                for term in subckt.neighbour:
                    assert self.termDict[term.nodeName] == term

        for term in self.termDict.itervalues():
            # Make sure there are no floating terminals
            assert term.neighbour
            for node in term.neighbour:
                if isinstance(node, Element):
                    assert self.elemDict[node.nodeName] == node
                else:
                    assert self.subcktDict[node.nodeName] == node


    # Actions on individual elements/terminals -------------------------
            
    def connect(self, element, termList):
        """
        Connect an Element (or subckt) instance to terminals specified by
        a list of terminal names (termList). If any terminal does not
        exist in the circuit it is created and added.
    
        **Order in the adjacency list is important for Elements**
    
        For example, for a MOSFET the first node corresponds to the drain,
        the second to the gate, etc.
        """
        # Some sanity checking. This function should be used once per element
        assert not element.neighbour
    
        for termName in termList:
            terminal = self.get_term(termName)
            terminal.neighbour.append(element)
            element.neighbour.append(terminal)
    

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
        if self.elemDict.has_key(elem.nodeName):
            raise CircuitError(elem.nodeName + ': Element already exists')
        else:
            self.elemDict[elem.nodeName] = elem
    
        if modelName:
            if elem.dotModel:
                raise CircuitError('{0} already has an ' +
                                   'assigned model'.format(elem.nodeName))
            # Public model: Check if model already in circuit
            if Circuit.modelDict.has_key(modelName):
                model = Circuit.modelDict[modelName]
                if model.modelType != elem.devType:
                    raise CircuitError(
                        'Incorrect model type "{0}" for element "{1}"'.format(
                            model.modelType, elem.nodeName))
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
            for n1 in elem.neighbour:
                # Disconnect any terminals from elem
                n1.neighbour.remove(elem)
                    
    def add_subckt(self, xsubckt):
        """
        Adds a Xsubckt instance to circuit
        """
        try:
            self.subcktDict[xsubckt.nodeName] = xsubckt 
        except KeyError:
            raise CircuitError('add_subckt: Subcircuit instance name "'\
                                + xsubckt.nodeName + '" already in use')
    
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
    Almost identical to Circuit but adds support for external
    connections.
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
            # Save terminal connection number here
            terminal.subCKTconnection = i

    def get_connections(self):
        """
        Returns a list of subcircuit terminals
        """
        return [self.get_term(termName) 
                for termName in enumerate(self.extConnectionList)]

    def netlist_string(self):
        """
        Just adds an extra .subckt line
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




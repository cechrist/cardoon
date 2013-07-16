"""
The main function available in this module is parse_file(). Also
analysisQueue is a list with all analyses to be performed.

Everything else should not be needed outside of the module.

If something goes wrong throws ParseError

--------------------------------------------------------------------
Copyright Carlos Christoffersen <c.christoffersen@ieee.org>

This file is part of the cardoon electronic circuit simulator.

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

http://www.gnu.org/licenses/gpl.html
"""

import pyparsing as pp

from paramset import ParamSet, ParamError
import circuit as cir
from globalVars import glVar
from devices import devClass
from analyses import anClass


class ParseError(Exception):
    """
    Used for exceptions raised in this module
    """
    pass

class NetVarException(Exception):
    """
    Used when a netlist variable is detected
    """
    pass

# Stack to hold main circuit while processing subcircuits
cktStack = []

# Queue to hold analyses to be processed
class Queue:
    def __init__(self):
        self.q = []

    def reset(self):
        self.q = []

analysisQueue = Queue()

# Characters allowed in parameter strings
allowedChars = '._+-*:'

# Grammar definition for numeric fields, using spice syntax mostly
#point = pp.Literal('.')
#e = pp.CaselessLiteral('e')
#plusorminus = pp.Literal('+') | pp.Literal('-')
#number = pp.Word(pp.nums) 
#integer = pp.Combine( pp.Optional(plusorminus) + number )
#floatnumber = pp.Combine(
#    ( (integer + pp.Optional(point + pp.Optional(number)))
#      | (pp.Optional(plusorminus) + point + number))
#    + pp.Optional( e + integer ) ) 

floatnumber = pp.Regex(r'[+-]?\d*(\.\d*)?([Ee][+-]?\d+)?')

# Multipliers
multT = pp.CaselessLiteral('t').setParseAction(lambda tok: 1e12)
multG = pp.CaselessLiteral('g').setParseAction(lambda tok: 1e9)
multMeg = pp.CaselessLiteral('meg').setParseAction(lambda tok: 1e6)
multK = pp.CaselessLiteral('k').setParseAction(lambda tok: 1e3)
multMil = pp.CaselessLiteral('mil').setParseAction(lambda tok: 25.4e-6)
multm = pp.CaselessLiteral('m').setParseAction(lambda tok: 1e-3)
multu = pp.CaselessLiteral('u').setParseAction(lambda tok: 1e-6)
multn = pp.CaselessLiteral('n').setParseAction(lambda tok: 1e-9)
multp = pp.CaselessLiteral('p').setParseAction(lambda tok: 1e-12)
multf = pp.CaselessLiteral('f').setParseAction(lambda tok: 1e-15)

multiplier = multT | multG | multMeg | multK \
    | multMil | multm | multu | multn | multp | multf

numvalue = floatnumber('value') \
    + pp.Optional(multiplier)('mult') \
    + pp.Optional(pp.Word(pp.alphas))

def string_to_number(valueString):
    """
    Quote from spice3 help file:
    A number field may be an integer  field  (12,  -44),  a
    floating  point field (3.14159), either an integer or float-
    ing point number followed by  an  integer  exponent  (1e-14,
    2.65e3),  or  either  an  integer or a floating point number
    followed by one of the following scale factors:
    
          12         9                    6         3               -6
    T = 10     G = 10             Meg = 10    K = 10      mil = 25.4
          -3                 -6         -9          -12         -15
    m = 10     u (or  M) = 10     n = 10      p = 10      f = 10
    
    Letters immediately following a number that  are  not  scale
    factors  are  ignored,  and  letters immediately following a
    scale factor are ignored.  Hence, 10, 10V, 10Volts, and 10Hz
    all  represent  the  same number, and M, MA, MSec, and MMhos
    all represent  the  same  scale  factor.   Note  that  1000,
    1000.0,  1000Hz,  1e3, 1.0e3, 1KHz, and 1K all represent the
    same number.
    """
    result = numvalue.parseString(valueString)
    value = float(result.value)
    if result.mult:
        value *= result.mult
    return value

def match_parameter(paramset, par):
    """
    Given a parameter set a tuple par=(param_name, value), it tries to
    match the type of value to the parameter type defined in paramset
    If this is not possible or the parameter does not exist it
    raises an exception
    """
    # Returns None if parameter does not exist
    try:
        ptype = paramset.get_type(par[0])
    except ParamError:
        raise ParseError('Unrecognized parameter: {0}\n'.format(par[0])
                         + 'Valid parameters: ' + paramset.list_parameters())
    else:
        try:
            convert_type(par, ptype)
        except NetVarException:
            # Netlist variable detected, let it go and let paramset to
            # take care
            pass


def convert_type(par, ptype):
    """
    To be used once the parameter type has been determined

    Raises NetVarException if a netlist variable is detected
    """
    if ptype == str:
        if par.vector:
            raise ParseError('"' + par[0] + 
                             '" must be a string: only the following '
                             'special characters allowed: ' + allowedChars)
    elif not par.vector:
        # Try to convert to ptype
        if ptype == float:
            try:
                par[1] = string_to_number(par[1])
            except (ValueError, pp.ParseException):
                raise NetVarException('"' + par[0] 
                                      + '" must be numeric: ' + par[1])
        elif ptype == int:
            try:
                par[1] = int(par[1])
            except ValueError:
                raise NetVarException('"' + par[0] 
                                      + '" must be an int: ' + par[1])
        elif ptype == bool:
            # Attempt to convert to a bool
            if (par[1].lower() == 'false') or (par[1] == '0'):
                par[1] = False
            elif (par[1].lower() == 'true') or (par[1] == '1'):
                par[1] = True
            else:
                raise NetVarException('Can not convert "' + par[0] 
                                      + '" value: ' + par[1]
                                      + ' to bool')
        elif ptype == list:
            raise NetVarException('Can not convert "' + par[0] + '" value: '
                                  + type(par[1]) + ' to list (vector)')
        else:
            assert False, par[0] + ': Parameter has unknown type'
    else:
        if ptype == list:
            # Convert to numeric values if possible
            v = []
            for item in par[1]:
                try:
                    val = string_to_number(item)   
                except (ValueError, pp.ParseException):
                    val = item
                v.append(val)
            par[1] = v
#            Old approach: force numeric values
#
#            # Assume the list contains all numeric values (otherwise
#            # it should have type==string in parameter definition)
#            try:
#                v = [string_to_number(val) for val in par[1]]
#                par[1] = v
#            except (ValueError, pp.ParseException):
#                raise ParseError('"', par[0] + 
#                                 '" must be a list of numeric values: ' 
#                                 + par[1])
        else:
            raise ParseError('Can not convert "' + par[0] + '" value: '
                             + type(par[1]))


# *******************************************************************
# Parse functions for individual lines
# *******************************************************************        

def parse_element(tok):
    # Try to create element
    try:
        dev = devClass[tok.devType](tok.instanceName)
    except KeyError:
        raise ParseError('Unknown device: ' + tok.devType)

    # Process paramenters first, as we need the model name (if any)
    modelname = None
    if tok.parameters:
        # first find if parameters are valid
        for par in tok.parameters:
            if par[0] == 'model':
                modelname = par[1]
            else:
                # try to match parameter type with expected type (if
                # parameter exists at all)
                match_parameter(dev, par)
                # if control reaches here we are OK
                dev.set_param(par[0], par[1])
                
    # Add to circuit and connect
    termList = tok.nodes
    cktStack[-1].add_elem(dev, modelname)
    cktStack[-1].connect(dev, termList)
    # Check that the number of terminals is correct, if possible
    dev.check_terms()


def parse_model(tok):
    model = cir.get_model(tok.modName)
    if not model:
        try:
            model = cir.Model(tok.modName,
                              tok.devType,
                              devClass[tok.devType].paramDict)
        except KeyError:
            raise ParseError('Could not create model "' + tok.modName 
                             + '" type: ' + tok.devType)
        else:
            cir.add_model(model)
    # Verify that model has the correct type
    if model.modelType != tok.devType:
         raise ParseError('Model "{0}" type {1} conflicts with another model of type {2}'.format(tok.modName, tok.devType, model.modelType))
    # Set parameters if any
    if tok.parameters:
        for par in tok.parameters:
            # try to match parameter type with expected type (if
            # parameter exists at all)
            match_parameter(model, par)
            # if control reaches here we are OK
            model.set_param(par[0], par[1])


def parse_options(tok):
    """
    Set attributes in glVars
    """
    # Set parameters 
    for par in tok.options:
        # try to match parameter type with expected type (if
        # parameter exists at all)
        match_parameter(glVar, par)
        # if control reaches here we are OK
        glVar.set_param(par[0], par[1])
    # Now make sure attributes are updated
    glVar.set_attributes(useDefaults = False)


def parse_vars(tok):
    """
    Set netlist variables in ParamSet.netVar
    """
    # Set parameters 
    for par in tok.vars:
        if par.vector:
            convert_type(par, list)
        elif (par[1] == 'True') or (par[1] == 'False'):
            convert_type(par, bool)
        else:
            # Convert type to float if possible
            convert_type(par, float)

        # if control reaches here we are OK
        ParamSet.netVar[par[0]] = par[1]


def parse_subcktDef(tok):
    # Try to create subcircuit
    subckt = cir.SubCircuit(tok.subName, tok.nodes)
    # Push subcircuit to stack (to be removed with .ends)
    cktStack.append(subckt)

def parse_subcktInst(tok):
    # Here we have to rerieve the subcircuit definition name from the
    # last token in 'nodes'
    defName = tok.nodes[-1]
    nodelist = tok.nodes[:-1]
    # import pdb; pdb.set_trace()
    # 1 connection node is OK if gnd is global: the code below not needed
    #    if len(nodelist) < 2:
    #        raise ParseError('Insufficient number of connections in '
    #                        + tok.instanceName + '. (Minimum is 2)')
    xckt = cir.Xsubckt(tok.instanceName, defName)
    cktStack[-1].add_subckt(xckt)
    cktStack[-1].connect(xckt, nodelist)

def parse_include(tok):
    parse_file(tok.filename, cktStack[-1])


def parse_analysis(tok):
    # Try to create element
    try:
        an = anClass[tok.anType]()
    except KeyError:
        raise ParseError('Unknown analysis: ' + tok.anType)
    # Process paramenters 
    if tok.parameters:
        # first find if parameters are valid
        for par in tok.parameters:
            # try to match parameter type with expected type (if
            # parameter exists at all)
            match_parameter(an, par)
            # if control reaches here we are OK
            an.set_param(par[0], par[1])

    # set attributes with parameter values plus defaults
    an.set_attributes()
    # Add analysis to queue
    analysisQueue.q.append(an)


def parse_plot(tok):
    # Create output request and add to circuit
    outreq = cir.OutRequest(tok.Type, tok.Vars)
    # Add requests always to main circuit
    cktStack[0].add_plot_request(outreq)

def parse_save(tok):
    # Create output request and add to circuit
    outreq = cir.OutRequest(tok.Type, tok.Vars)
    # Add requests always to main circuit
    cktStack[0].add_save_request(outreq)
        
def parse_ends(tok):
    # First make sure that we were indeed processing a subcircuit
    if (len(cktStack) > 1) and (hasattr(cktStack[-1], 'extConnectionList')):
        # Ok, the thing resembles a subcircuit
        cktStack.pop()
    else:
        raise ParseError('.ends found but not inside a subcircuit')


def parse_end(tok):
    # Empty stack. Stop processing everything and return
    while len(cktStack):
        cktStack.pop()


# *******************************************************************

def parse_file(filename, ckt):
    """
    Parse statements contained in the file pointed by filename and add
    to ckt. For now the parser is case-sensitive except for the dotted
    keywords (.model, etc.).
    """
    # *******************************************************************
    # Define some grammar
    # *******************************************************************
    # Defined here just in case (because function must be re-entrant)
    
    # Symbol definitions
    LPAR,RPAR,LBRACK,RBRACK,EQ,COLON,QUOTE = map(pp.Suppress,"()[]=:'")

    # ParString used for string parameter values
    parString = pp.Word(pp.alphanums + allowedChars) \
        | QUOTE + pp.Word(pp.alphanums + allowedChars + ' ()') + QUOTE

    vector = LBRACK + pp.delimitedList(parString) + RBRACK

    parName = pp.Word(pp.alphas, pp.alphanums + '_')
    parameters = pp.OneOrMore(
        pp.Group(parName + EQ + 
                 (parString('single') | pp.Group(vector)('vector'))))

    # Used for names of: devices (type and instances), nodes, subcircuits
    identifier = pp.Word(pp.alphanums + '_-')

    elemname = identifier('devType') + pp.Suppress(':') \
        + identifier('instanceName')
    
    nodes = pp.OneOrMore(identifier + ~pp.FollowedBy("="))

    # Output variables (for .plot, .save)
    intTerm = pp.Group(pp.Combine(identifier + pp.Literal(':') + identifier)
                       + COLON + identifier)
    oneVar = intTerm | identifier
    outVars = pp.OneOrMore(oneVar)
    # Comment line: any line that starts with # , * or //
    commentline = pp.Suppress(((pp.Literal('*') ^ pp.Literal('#')) 
                               + pp.Regex('.*')) ^ pp.dblSlashComment)
    
    # example: diode:d1 1 gnd model=mydiode isat= 2e-15 area = 2.
    elemline = elemname + nodes('nodes') + pp.Optional(parameters('parameters'))
    
    # example: .model mydiode diode ( isat=1e-15 cj0=5e-12 )
    modelline = pp.Suppress(pp.Keyword('.model', caseless=True)) \
        + identifier('modName') \
        + identifier('devType') \
        + pp.Optional(LPAR + parameters('parameters') + RPAR )

    # example: .options abstol=1e-8
    optionsline = pp.Suppress(pp.Keyword('.options', caseless=True)) \
        + parameters('options') 

    # example: .vars freq=1GHz
    varsline = pp.Suppress(pp.Keyword('.vars', caseless=True)) \
        + parameters('vars') 

    # example: .subckt LM741 in out vdd vee
    subcktDefLine = pp.Suppress(pp.Keyword('.subckt', caseless=True)) \
        + identifier('subName') + nodes('nodes')

    # example: xamp1 2 5 1 3 LM741
    # Treat the last node as the subcircuit definition name. Sorry
    subcktInstLine = pp.Word('xX', pp.alphanums + '_')('instanceName') \
        + nodes('nodes')

    # example: .include model.net
    fileString = pp.Word(pp.alphanums + '._-+')
    includeline = pp.Suppress(pp.Keyword('.include', caseless=True)) \
        + fileString('filename')

    # example: .analysis op
    analysisline = pp.Suppress(pp.Keyword('.analysis', caseless=True)) \
        + identifier('anType') + pp.Optional(parameters('parameters'))

    # example: .plot dc 10 out1
    plotline = pp.Suppress(pp.Keyword('.plot', caseless=True)) \
        + identifier('Type') + outVars('Vars') 

    # example: .save dc 10 out1
    saveline = pp.Suppress(pp.Keyword('.save', caseless=True)) \
        + identifier('Type') + outVars('Vars') 

    endsline = pp.Keyword('.ends', caseless=True)
    
    endline = pp.Keyword('.end', caseless=True)

    netlistLine = commentline \
        | elemline.setParseAction(parse_element) \
        | modelline.setParseAction(parse_model) \
        | optionsline.setParseAction(parse_options) \
        | varsline.setParseAction(parse_vars) \
        | subcktDefLine.setParseAction(parse_subcktDef) \
        | subcktInstLine.setParseAction(parse_subcktInst) \
        | includeline.setParseAction(parse_include) \
        | analysisline.setParseAction(parse_analysis) \
        | plotline.setParseAction(parse_plot) \
        | saveline.setParseAction(parse_save) \
        | endsline.setParseAction(parse_ends) \
        | endline.setParseAction(parse_end)

    # Each time this function is called it puts the working ckt in the
    # stack and takes it out when finished. So at any time we can
    # access the current circuit as cktStack[-1]
    cktStack.append(ckt)
    # Clean queue
    analysisQueue.reset()
    # Save file name in circuit
    ckt.filename = filename

    try:
        with open(filename, 'r') as f:
            lineNumber = 0
            lineAcc = ''
        
            for line in f:
                lineNumber += 1
                # Remove unneeded spaces from line
                line = line.strip()
                # Prepend lineAcc value (in case previous line ended with '\')
                if lineAcc:
                    line = lineAcc + ' ' + line
        
                if not line:
                    # nothing left to parse in this line, go to the next
                    continue
        
                # First line is the main circuit title (to follow
                # spice tradition)
                if (ckt.name == 'main') and (lineNumber == 1) \
                        and (len(cktStack) == 1):
                    ckt.title = line
                    continue
                
                # Join consecutive lines if line ends with '\'
                if line[-1] == '\\':
                    # remove backslash
                    line = line[:-1]
                    # must read next line before continuing
                    lineAcc = line
                    continue
                else:
                    # Reset the accumulator
                    lineAcc = ''
        
                # Most of the work is made here
                try: 
                    result = netlistLine.parseString(line, parseAll = True)
                except (ParseError, NetVarException) as pe:
                    raise ParseError('Parse error in file ' + filename
                                     + ', line ' + str(lineNumber) 
                                     + ':\n' + str(pe) + '\n"' + line + '"') 
                except cir.CircuitError as ce:
                    raise ParseError('Circuit error in file ' + filename
                                     + ', line ' + str(lineNumber) 
                                     + ':\n' + str(ce) + '\n"' + line + '"') 
                except pp.ParseException as pe:
                    mesg = 'Syntax error in file ' + filename
                    mesg += ', line ' + str(lineNumber) + '\n' 
                    mesg += '"' + line + '"'
                    if len(line) < 80:
                        mesg += '\n' + pe.col * '-' + '^'
                    raise ParseError(mesg) 
                
                if not cktStack:
                    # This means that .end was found and processing must end
                    # immediatly
                    break
                    
            # Take circuit out (if stack not empty due to .end)
            if cktStack:
                cktStack.pop()

    except IOError as ioe:
        raise ParseError('Parse error -> ' + str(ioe))

    return analysisQueue.q

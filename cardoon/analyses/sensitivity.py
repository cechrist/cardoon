"""
:mod:`sensitivity` -- Nodal sensitivity analysis
------------------------------------------------

.. module:: sensitivity
.. moduleauthor:: Carlos Christoffersen

Functions to calculate derivatives for sensitivity analysis. Elements
are assumed to have nodal arguments (``nD_*``).
"""

import numpy as np
import pycppad as ad
import nodal as nd

def get_i_derivatives(device, parList, xVec, force=False):
    """
    Returns DC current derivatives with respect to given parameters

    Inputs: 

       device (linear or nonlinear)
       parList: list of parameters as returned by param.get_float_attributes()
       xVec: nodal vector
       force: if True, force tape re-generation

    Outputs: 

       output currents Jacobian

    Possible side effect: operating point attributes lost in internal
    terminals if tape is created
    """
    # Create input vector from parameter values
    pVec = np.array([val for name, val in parList])
    inVec = pVec

    if device.isNonlinear:
        # Add controlling voltages to input: AD tape can be reused for
        # different operating points
        xin = np.zeros(device.nD_nxin)
        nd.set_xin(xin, device.nD_vpos, device.nD_vneg, xVec)
        #assert(np.max(xin - device.nD_xOP) == 0)
        inVec = np.concatenate((pVec, xin), axis=0)

    # Use taped function if it exists
    if (not hasattr(device, '_sensF')) or force:
        device._sensF = create_sensitivity_tape(device, parList, inVec)

    f = device._sensF

    # Generate current derivatives: must use forward call only because
    # reverse creates problems with calculations using condassign in
    # the wrong way (no longer needed, but forward still OK because no
    # need to calculate derivatives respect to nonlinear control
    # voltages.
    f.forward(0, inVec)
    w = np.zeros_like(inVec)
    Ji = np.zeros((len(xVec), len(pVec)), dtype = float)
    for i in range(len(pVec)):
        # In the future w[i] can be set equal to delta_p then forward
        # returns the delta_i
        w[i] = 1.
        doutdp = f.forward(1, w)
        #print(f.forward(1,w))
        # Combine derivatives in output Jacobian
        if device.linearVCCS:
            for k,vccs in enumerate(device.nD_linVCCS):
                # Controlling voltage: xVec[col1] - xVec[col2]
                vc = 0.
                col = vccs[1]
                if col >= 0:
                    vc = xVec[col]
                col = vccs[3]
                if col >= 0:
                    vc -= xVec[col]
                # Output current derivative = dg/dp * vc
                di = vc * doutdp[k]
                row = vccs[0]
                if row != -1:
                    Ji[row, i] += di
                row = vccs[2]
                if row != -1:
                    Ji[row, i] -= di
        if device.isDCSource:
            row = device.sourceOutput
            k = len(device.nD_linVCCS)
            if row[0] >= 0:
                Ji[row[0], i] -= doutdp[k]
            if row[1] >= 0:
                Ji[row[1], i] += doutdp[k]
        if device.isNonlinear:
            diNL = doutdp[-len(device.csOutPorts):]
            nd.set_i(Ji[:,i], device.nD_cpos, device.nD_cneg, diNL)
        # Unselect column
        w[i] = 0.

    #print(f.jacobian(inVec))
    return Ji

def create_sensitivity_tape(device, parList, inVec):
    """
    Create sensitivity AD tape

    Inputs:

        device: device with nodal attributes 
        parList: list of parameters to calculate sensitivities

        inVec: input vector including parameters and perhaps nonlinear
        control voltages at the end of the vector

    Side effects: operating point attributes lost in internal terminals
    """
    # Should store some parameter information for safety
    #
    # Create AD tape -----------------------------------------------------
    # 
    a_inVec = ad.independent(inVec)
    # set adouble attributes: a_inVec may be longer than parList but
    # zip() truncates to the shortest list
    for item, a_val in zip(parList, a_inVec):
        setattr(device, item[0], a_val)
    # Re-calculate coductances / parameters
    linVec = np.empty(shape=0)
    device.process_params()
    if device.linearVCCS:
        gList = []
        for vccs in device.linearVCCS:
            # Linear output current
            g = vccs[2]
            # Kludge: Make sure all elements are of the correct type
            if type(g) != ad.cppad_.a_float:
                g = ad.cppad_.a_float(g)
            gList.append(g)
        # Overwrite output vector
        linVec = np.array(gList)
    #import pdb; pdb.set_trace()

    # Later add for time- and frequency- domain
    if device.isDCSource:
        val = device.get_DCsource()
        # Kludge: Make sure all elements are of the correct type
        if type(val) != ad.cppad_.a_float:
            val = ad.cppad_.a_float(val)
        valVec = np.array([val])
        sourceVec = np.concatenate((linVec, valVec), axis=0)
    else:
        sourceVec = linVec
    
    if device.isNonlinear:
        # calculate nonlinear currents and concatenate to output
        (iVec, qVec) = device.eval_cqs(a_inVec[-device.nD_nxin:])
#        # Kludge: Make sure all elements are of the correct type (not needed?)
#        for k,i in enumerate(iVec):
#            if type(i) != ad.cppad_.a_float:
#                iVec[k] = ad.cppad_.a_float(i)
        a_outVec = np.concatenate((sourceVec, iVec), axis=0)
    else:
        # Just copy whatever we have
        a_outVec = sourceVec

    # tape stored in f
    f = ad.adfun(a_inVec, a_outVec)
    # Restore element to original state
    device.clean_attributes()
    device.set_attributes()
    device.process_params()
    nd.restore_RCnumbers(device)
    # End of crate tape ------------------------------------------------
    return f


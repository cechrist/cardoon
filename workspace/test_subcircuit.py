circuit = cir.get_mainckt()
analysisQueue = parse_net('lma411_cascade.net', circuit)
import analyses.nodal
circuit.init()
analyses.nodal.make_nodal_circuit(circuit)
print circuit.termDict
lma=circuit.cktDict['lma411']
analyses.nodal.make_nodal_circuit(lma)


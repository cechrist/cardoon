
Cardoon Electronic Circuit Simulator
====================================

Introduction
------------

Cardoon is a circuit simulator program being developed using
Python. The design is inspired from the ideas that worked well in
fREEDA <http://www.freeda.org/> and carrot
<http://vision.lakeheadu.ca/research.html>, and some improvements that
take advantage of the flexibility in Python. 

One of the main goals of this project is to obtain a program which is
easy to modify and test to experiment new simulation
algorithms/models. A secondary objective is to add features and make
the program efficient enough to be useful for general use, if
possible.

At this time only a few models are implemented, but some of them are
already more complete than the equivalents in freeda and carrot (such
as MOS EKV 2.6 and the BJT model). State- variable modelling (as in
freeda) is supported. Currently a state-variable diode model is
included. Electro-thermal models are automatically generated from
electrical models.  Nonlinear models use automatic differentiation for
jacobian calculation. In addition, automatic differentiation can also
be used to calculate sensitivities.

So far there is no real analysis implemented, except a simple
nonlinear device test routine. The immediate plan is to document and
consolidate the main design and after that implement simple
DC/AC/Transient(/noise?) to verify the Element interface works as
intended.

Cardoon is being developed with Python 2.6 on Debian Squeeze, but it
should work OK with Python 2.7 under other operating systems.

Comments are welcome!

Carlos Christoffersen <c.christoffersen@ieee.org>


License
-------

Cardoon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 or later:

<http://www.gnu.org/licenses/gpl.html>

(See also LICENSE file included with the source)



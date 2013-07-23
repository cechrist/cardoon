#!/usr/bin/env python
from __future__ import print_function
from distutils.core import setup
# Perhaps this will fail if cardoon already installed
from cardoon.simulator import release

setup(name='cardoon',
      version=release,
      description='Cardoon electronic circuit simulation library',
      author='Carlos Christoffersen',
      author_email='cechrist@vision.lakeheadu.ca',
      url='http://vision.lakeheadu.ca/cardoon/',
      packages=['cardoon', 'cardoon.devices', 'cardoon.analyses'],
     )

print("""

===================================================================

IMPORTANT: Before this software can be used, the required libraries
must be installed. Please read doc/install.rst

===================================================================
""")



Installation
============

Installation on GNU/Linux
-------------------------

Cardoon is being currently being developed in Debian GNU/Linux 7.1
(Wheezy) <http://www.debian.org>, but it should work in other
environments where the Python interpreter and the libraries are
available.  In addition to a Python interpreter (version 2.6 or 2.7),
several free libraries are required to run cardoon. In Debian Linux or
Debian-based distributions (Ubuntu, Linux Mint and others) all of
these libraries except pycppad can be installed as packages::

  apt-get install libsuperlu3 python-numpy python-scipy python-matplotlib \
  python-pyparsing ipython python-sphinx pylint libboost-python-dev

To install pycppad please follow the instructions on its website
(listed below).

The corresponding websites for each library are listed below:

* pycppad: Automatic differentiation

  - cppad:  http://www.coin-or.org/CppAD/Doc/cppad.xml

  - pycppad:  http://www.seanet.com/~bradbell/pycppad/pycppad.xml 

* numpy:  http://numpy.scipy.org/ (matrix and vector support)

* scipy:  http://www.scipy.org/ (sparse matrix support, FFT, etc.)

* matplotlib:  http://matplotlib.sourceforge.net/ (plotting)

* pyparsing:  http://pyparsing.wikispaces.com/ (parser, cardoon tested
  with up to version 1.5.2)

* ipython:  http://ipython.org/ (iterative shells)

* Sphinx: http://sphinx.pocoo.org/ (documentation)

* pyreverse from the pylint package: http://www.logilab.org/2560 (to
  generate UML diagrams)

Most of these libraries in turn depend on other libraries They should
be automatically installed. Some must be installed manually such as
libsuperlu3 (http://crd-legacy.lbl.gov/~xiaoye/SuperLU/) because it is
not in the python-scipy dependencies.

You may also install *git* to fetch the source code from the github
repository::

    git clone git://github.com/cechrist/cardoon.git

Alternatively the source code can be downloaded as a zip file from:

https://github.com/cechrist/cardoon

The simulator can be run directly from the source code tree by putting
the following script somewhere in your path (assume that the source
code is in ``/home/user/src/cardoon``)::

  # Sample script to run cardoon from the command line
  PYTHONPATH=$PYTHONPATH:/home/user/src/cardoon
  python -m cardoon $*

Alternatively, a ``setup.py`` script is provided. If you would want to
install the cardoon library system-wide, from the cardoon directory,
install as follows::

    python ./setup.py install --prefix=/usr/local

You can replace ``/usr/local`` by another directory. For example to
install in your home directory::

    python ./setup.py install --prefix=$HOME

If you install in a custom directory, be sure to set the
``PYTHONPATH`` variable to tell the python interpreter where to look
for the libraries. For example, if you installed with
``--prefix=/home/user`` then::

    $ PYTHONPATH=$PYTHONPATH:/home/user/lib/site-packages

To try if the installation was successful, run the interpreter and try
importing the library::

    cechrist@venus:~$ python
    Python 2.7.3 (default, Jan  2 2013, 13:56:14) 
    [GCC 4.7.2] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import cardoon.simulator 
    >>> 

You should get no error messages.


Installation on MS Windows
--------------------------

There are many possible setups, one of them is discussed here. If you
are not building the pycppad library, download and install
WinPython-32-bit-2.7.5.1 from http://code.google.com/p/winpython/ .
WinPython comes with most of the required libraries ready to use. Only
two additional libraries are needed:

* pycppad: Automatic differentiation

  - cppad:  http://www.coin-or.org/CppAD/Doc/cppad.xml

  - pycppad:  http://www.seanet.com/~bradbell/pycppad/pycppad.xml 

  - pycppad uses the boost.python library http://www.boost.org/libs/python/

  A pre-compiled version to be used with WinPython-32-bit-2.7.5.1 is
  available at: 
  http://vision.lakeheadu.ca/cardoon/pycppad-WinPython-32bit-2.7.5.1.zip

  Open the zip file and put all files in the ``python-2.7.5``
  subdirectory of the WinPython installation::

      cd WinPython-32bit-2.7.5.1\python-2.7.5
      unzip c:\Users\user1\Downloads\pycppad-WinPython-32bit-2.7.5.1.zip

  pycppad should be ready to use after this step. If you would rather
  compile the libraries yourself, some instructions are provided at:
  http://list.coin-or.org/pipermail/cppad/2013q2/000309.html

* pyparsing:  http://pyparsing.wikispaces.com/ (parser)

  Download the source from the official site (tested up to version
  1.5.2), open in some directory and run the following in a "WinPython
  Command Prompt" window::

    python setup.py install

  pyparsing should be ready to use after this step.

The source code for the cardoon simulator can be downloaded as a zip
file from:

https://github.com/cechrist/cardoon

Unpack the zip file in some directory. The simulator can be run
directly from the source code tree by putting the following batch file
somewhere in your path (assume that the source code is in
``d:\cardoon``)::

  @echo off
  rem Sample script to run cardoon from the command line (cardoon.bat)
  rem Replace d:\cardoon by the actual directory
  PYTHONPATH=%PYTHONPATH%:d:\cardoon
  python -m cardoon %*

Alternatively, a ``setup.py`` script is provided. If you would want to
install the cardoon library system-wide, from the cardoon directory,
install as follows::

    python setup.py install 

If you install in a custom directory, be sure to set the
``PYTHONPATH`` variable to tell the python interpreter where to look
for the libraries. 

To try if the installation was successful, run the interpreter and try
importing the library::


    D:\>python
    Python 2.7.5 (default, May 15 2013, 22:43:36) [MSC v.1500 32 bit (Intel)] on win
    32
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import cardoon.simulator
    >>> 

You should get no error messages.


Generating this documentation
-----------------------------

The main documentation files are kept in the ``doc``
directory. Documentation can be generated in html or LaTeX formats
(other formats are possible but not tested).  The documentation can be
generated as follows::

    cd doc
    make html


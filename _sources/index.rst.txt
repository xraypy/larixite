.. larixite  documentation master file

Larixite: Crystal structures for X-ray absorption spectroscopy
===============================================================

.. _xraylarch: https://xraypy.github.io/xraylarch/
.. _Larixite Webapp: https://millenia.cars.aps.anl.gov/larixite
.. _American Mineralogical Crystal Structure Database: https://rruff.geo.arizona.edu/AMS/amcsd.php

Larixite is a Python package to help use crystallographic data to build
clusters of atoms to generate inputs for X-ray absorption spectroscopy and
other scientific disciplines that use non-crystalline clusters of atoms.

This project includes:

1. an sqlite3 database of structures from the `American Mineralogical Crystal Structure Database`_ (AMCSD).
2. Python code to convert structures from the AMCSD database, other CIF files,
   or XYZ coordinates into atomic clusters for XAS calculations with FEFF,
   FDMNES, and other XAS calculation tools.
3. A basic web application to guide those conversions. See `Larixite WebApp`_.


Installation
------------------

Either install from PyPI with::

    > pip install larixite


or download and unpack this code and install with::


    > pip install .


Status
--------

Larixite has been in rapid development, but is also a spin-off from code that
has been in `xraylarch`_ for many years.  That is, while many parts of the code
are moving rapidly, much of the code is reasonably stable.


Web App
----------


The `Larixite WebApp`_ can be run locally for debugging or for local
deployment.  To do this, install the extra wed dependencies (essentially only
Flask is needed) with::

    > pip install ".[web]"


and run the script "run_local.py" with::

      > python run_local.py

will launch a local web server with the app running at http://127.0.0.1:11564/

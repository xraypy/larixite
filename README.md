# Larixite

Crystal structures and clusters of atoms for X-ray absorption spectroscopy.

The main purpose of larixite is to provide a Python package for using
crystallographic data or calculated clusters of atoms to generate inputs for
X-ray absorption spectroscopy and other scientific disciplines that use
non-crystalline clusters of atoms.

This project includes:

1. an sqlite3 database of structures from the [American Mineralogical
   Crystal Structure Database](https://rruff.geo.arizona.edu/AMS/amcsd.php) (AMCSD)
2. Python code to convert structures from the AMCSD database, other CIF files,
   or XYZ coordinates into atomic clusters for XAS calculations with FEFF,
   FDMNES, and other XAS calculation tools.
3. A basic web application to guide those conversions. See [Larixite Web App](https://millenia.cars.aps.anl.gov/larixite).


## install

Either install from PyPI with

    > pip install larixite


Download and unpack this code and install with


    > pip install .


## Status

Larixite has been in rapid development, but is also a spin-off from code that
has been in Xraylarch for many years.  That is, while many parts of the code
are moving rapidly, much of the code is reasonably stable.


## Web App

The [Larixite Web App](https://millenia.cars.aps.anl.gov/larixite) can be run
locally for debugging or for local deployment.  To do this, install the extra
wed dependencies (essentially only Flask is needed) with

    > pip install ".[web]"


and run the script "run_local.py" with

      > python run_local.py

will launch a local web server with the app running at http://127.0.0.1:11564/

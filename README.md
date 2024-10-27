# Xstructures

Crystal structures and clusters of atoms.

The main purpose of xstructures is to provide a Python package for using
crystallographic data or calculated clusters of atoms to generate inputs
for  X-ray absorption spectroscopy and other scientific disciplines that use
non-crystalline clusters of atoms.


This project includes

1. an sqlite3 database of crystal structure from the American Mineralogical
   Crystal Structure Database (AMCSD)
2. Python code to convert structures from the AMCSD database, other CIF files,
    or XYZ coordinates into atomic clusters for XAS calculations with FEFF,
    FDMNES, and other XAS calculation tools.
3. A basic web application to guide those conversions.


## install

Download and unpack this code and install with


    > pip install .


## Status

Xstructures is the alpha-level of development.  But, the project is also a
spin-off from code that has been in Xraylarch for many years, so many parts of
the code are reasonably stable.


## Web App

To run the web app locally,  install the extra we dependencies with

    > pip install ".[web]"


and run the script "run_local.py" with

      > python run_local.py

will launch a local web server with the app running at http://127.0.0.1:11564/

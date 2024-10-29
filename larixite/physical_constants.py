# Useful physical constants
# most of these are put into common X-ray units (Angstroms, ev)

import scipy.constants as consts
from numpy import pi

RAD2DEG = 180.0/pi
DEG2RAD = pi/180.0
PI = pi
TAU = 2*pi

# electron rest mass in eV
E_MASS = consts.electron_mass * consts.c**2 / consts.e

# Planck's Constant
#   h*c    ~= 12398.42 eV*Ang
#   hbar*c ~=  1973.27 eV*Ang
PLANCK_HC    = 1.e10 * consts.Planck * consts.c / consts.e
PLANCK_HBARC = PLANCK_HC / TAU

# will be able to import these from xraydb when v 4.5.1 is required
ATOM_SYMS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
           'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
           'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
           'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
           'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La',
           'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
           'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
           'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
           'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md',
           'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
           'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

ATOM_NAMES = ['hydrogen', 'helium', 'lithium', 'beryllium', 'boron', 'carbon',
            'nitrogen', 'oxygen', 'fluorine', 'neon', 'sodium', 'magnesium',
            'aluminum', 'silicon', 'phosphorus', 'sulfur', 'chlorine', 'argon',
            'potassium', 'calcium', 'scandium', 'titanium', 'vanadium',
            'chromium', 'manganese', 'iron', 'cobalt', 'nickel', 'copper',
            'zinc', 'gallium', 'germanium', 'arsenic', 'selenium', 'bromine',
            'krypton', 'rubidium', 'strontium', 'yttrium', 'zirconium',
            'niobium', 'molybdenum', 'technetium', 'ruthenium', 'rhodium',
            'palladium', 'silver', 'cadmium', 'indium', 'tin', 'antimony',
            'tellurium', 'iodine', 'xenon', 'cesium', 'barium', 'lanthanum',
            'cerium', 'praseodymium', 'neodymium', 'promethium', 'samarium',
            'europium', 'gadolinium', 'terbium', 'dysprosium', 'holmium',
            'erbium', 'thulium', 'ytterbium', 'lutetium', 'hafnium',
            'tantalum', 'tungsten', 'rhenium', 'osmium', 'iridium', 'platinum',
            'gold', 'mercury', 'thallium', 'lead', 'bismuth', 'polonium',
            'astatine', 'radon', 'francium', 'radium', 'actinium', 'thorium',
            'protactinium', 'uranium', 'neptunium', 'plutonium', 'americium',
            'curium', 'berkelium', 'californium', 'einsteinium', 'fermium',
            'mendelevium', 'nobelium', 'lawrencium', 'rutherfordium',
            'dubnium', 'seaborgium', 'bohrium', 'hassium', 'meitnerium',
            'darmstadtium', 'roentgenium', 'copernicium', 'nihonium',
            'flerovium', 'moscovium', 'livermorium', 'tennessine', 'oganesson']

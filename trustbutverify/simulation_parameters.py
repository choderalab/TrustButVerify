"""
Default parameter values for MD simulations.
"""
from simtk import unit as u

PADDING = 1.0 * u.nanometers
CUTOFF = 0.95 * u.nanometers
PRESSURE = 1.0 * u.atmospheres
FRICTION = 1.0 / u.picoseconds
EQUILIBRATION_FRICTION = 1.0 / u.picoseconds
BAROSTAT_FREQUENCY = 25

EQUILIBRATION_TIMESTEP = 1.0 * u.femtoseconds
TIMESTEP = 2.0 * u.femtoseconds

base_path = "/home/kyleb/dat/TrustButVerify/"

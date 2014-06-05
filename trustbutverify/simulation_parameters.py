"""
Default parameter values for MD simulations.
"""
import os
from simtk import unit as u

PADDING = 1.0 * u.nanometers
CUTOFF = 0.95 * u.nanometers
PRESSURE = 1.0 * u.atmospheres
FRICTION = 1.0 / u.picoseconds
EQUIL_FRICTION = 1.0 / u.picoseconds
BAROSTAT_FREQUENCY = 25

EQUIL_TIMESTEP = 1.0 * u.femtoseconds
TIMESTEP = 2.0 * u.femtoseconds

N_STEPS = 10000
N_EQUIL_STEPS = 10000

EQUIL_OUTPUT_FREQUENCY = 5000
PROTEIN_OUTPUT_FREQUENCY = 1000
OUTPUT_FREQUENCY = PROTEIN_OUTPUT_FREQUENCY * 10

base_path = os.path.join(os.environ["HOME"], "/dat/TrustButVerify/")

water_mapping = {"tip3p":"tip3p", "tip4pew":"tip4pew", "tip3p-fb":"tip3p", "tip4p-fb":"tip4pew"}  # This is a hack for OpenMM to map water models to their "base" names.

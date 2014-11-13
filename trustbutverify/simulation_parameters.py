"""
Default parameter values for MD simulations.
"""
import os
from simtk import unit as u

pH = 7.0  # Default

PADDING = 1.0 * u.nanometers
CUTOFF = 0.95 * u.nanometers
PRESSURE = 1.0 * u.atmospheres
FRICTION = 1.0 / u.picoseconds
EQUIL_FRICTION = 10.0 / u.picoseconds
BAROSTAT_FREQUENCY = 25

EQUIL_TIMESTEP = 1.0 * u.femtoseconds
TIMESTEP = 2.0 * u.femtoseconds

N_STEPS = 25000000
N_STEPS_PEPTIDES = 25000000
N_EQUIL_STEPS = 25000

EQUIL_OUTPUT_FREQUENCY = 5000

PROTEIN_OUTPUT_FREQUENCY = 50000  # This is for just protein atoms
OUTPUT_FREQUENCY = PROTEIN_OUTPUT_FREQUENCY * 10  # This is for all atoms, including solvent

PROTEIN_OUTPUT_FREQUENCY_PEPTIDES = 5000
OUTPUT_FREQUENCY_PEPTIDES = PROTEIN_OUTPUT_FREQUENCY_PEPTIDES * 10



base_path = os.path.join(os.environ["HOME"], "dat/TrustButVerify50ns/")
CS_CACHE_PATH = os.path.join(base_path, "cached_chemical_shifts/")


water_mapping = {"tip3p":"tip3p", "tip4pew":"tip4pew", "tip3pfb":"tip3p", "tip4pfb":"tip4pew"}  # This is a hack for OpenMM to map water models to their "base" names.

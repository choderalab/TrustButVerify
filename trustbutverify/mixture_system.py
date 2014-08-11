import tempfile
import numpy as np
import os
import mdtraj as md

import pdbfixer
import pdbbuilder
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

from .target import Target
from .simulation_parameters import *
import utils

from protein_system import System
import gaff2xml

N_STEPS_MIXTURES = 100000000
N_EQUIL_STEPS_MIXTURES = 5000000
OUTPUT_FREQUENCY_MIXTURES = 500

class MixtureSystem(System):
    def __init__(self, smiles_strings, n_monomers, temperature, pressure, output_frequency=OUTPUT_FREQUENCY_MIXTURES, n_steps=N_STEPS_MIXTURES, equil_output_frequency=OUTPUT_FREQUENCY_MIXTURES, **kwargs):
        super(MixtureSystem, self).__init__(temperature=temperature, pressure=pressure, output_frequency=output_frequency, n_steps=n_steps, equil_output_frequency=equil_output_frequency, **kwargs)
        
        self.smiles_strings = smiles_strings
        self._target_name = "hey"
        self.n_monomers = n_monomers

    def build(self):
        self.box_pdb_filename = self.identifier + ".pdb"
        utils.make_path(self.box_pdb_filename)

        if os.path.exists(self.box_pdb_filename):
            return
        
        self.monomer_trajectories = []
        self.pdb_filenames = []
        ligand_trajectories, self.ffxml = gaff2xml.utils.smiles_to_mdtraj_ffxml(self.smiles_strings)    
        for k, ligand_traj in enumerate(ligand_trajectories):
            pdb_filename = tempfile.mktemp(suffix=".pdb")
            ligand_traj.save(pdb_filename)
            self.pdb_filenames.append(pdb_filename)

        self.packed_trj = gaff2xml.packmol.pack_box(self.pdb_filenames, self.n_monomers)
        self.packed_trj.save(self.box_pdb_filename)
        
    def equilibrate(self):
        self.equil_dcd_filename = "equil.dcd"
        self.equil_pdb_filename = "equil.pdb"
        utils.make_path(self.equil_pdb_filename)
        
        if os.path.exists(self.equil_pdb_filename):
            return

        positions = self.packed_trj.openmm_positions(0)
        topology = self.packed_trj.top.to_openmm()
        topology.setUnitCellDimensions(mm.Vec3(*self.packed_trj.unitcell_lengths[0]) * u.nanometer)
        
        self.ffxml.seek(0)
        ff = app.ForceField(self.ffxml)

        system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=self.cutoff, constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(self.temperature, self.equil_friction, self.equil_timestep)
        system.addForce(mm.MonteCarloBarostat(self.pressure, self.temperature, self.barostat_frequency))

        simulation = app.Simulation(topology, system, integrator, platform=self.platform)
        simulation.context.setPositions(positions)
        
        print('Minimizing.')
        simulation.minimizeEnergy()

        simulation.context.setVelocitiesToTemperature(self.temperature)
        print('Equilibrating.')

        simulation.reporters.append(app.DCDReporter(self.equil_dcd_filename, self.equil_output_frequency))
        simulation.step(self.n_equil_steps)
        
        # Re-write a better PDB with correct box sizes.
        traj = md.load(self.equil_dcd_filename, top=self.box_pdb_filename)[-1]
        traj.save(self.equil_pdb_filename)


    def production(self):
        self.production_dcd_filename = "production.dcd"        

        utils.make_path(self.production_dcd_filename)

        if os.path.exists(self.production_dcd_filename):
            return
        
        self.ffxml.seek(0)
        ff = app.ForceField(self.ffxml)
        
        pdb = app.PDBFile(self.equil_pdb_filename)
        
        system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=self.cutoff, constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(self.temperature, self.friction, self.timestep)
        system.addForce(mm.MonteCarloBarostat(self.pressure, self.temperature, self.barostat_frequency))

        simulation = app.Simulation(pdb.topology, system, integrator, platform=self.platform)
        simulation.context.setPositions(pdb.positions)

        simulation.context.setVelocitiesToTemperature(self.temperature)
        print('Production.')
        simulation.reporters.append(app.DCDReporter(self.production_dcd_filename, self.output_frequency))
        simulation.step(self.n_steps)

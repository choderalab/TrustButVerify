import numpy as np
import os
import mdtraj as md

import pdbfixer
from simtk.openmm import app
import simtk.openmm as mm

from .target import Target
from .simulation_parameters import *
import utils


class ProteinTarget(Target):
    def __init__(self, pdb_id, temperature, padding=PADDING, cutoff=CUTOFF, keep_chains=None, pressure=PRESSURE, timestep=TIMESTEP, equilibration_timestep=EQUILIBRATION_TIMESTEP):
        self.pdb_id = pdb_id
        self._target_name = pdb_id
        
        if keep_chains is None:
            keep_chains = np.arange(1)
        
        self.keep_chains = keep_chains


    def build(self, ff_name, water_name):
        out_filename = self.get_initial_pdb_filename(ff_name, water_name)
        utils.make_path(out_filename)

        if os.path.exists(out_filename):
            return
        
        fixer = pdbfixer.PDBFixer(pdbid=self.pdb_id)

        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.removeHeterogens(True)
        fixer.addMissingHydrogens()

        n_chains = len(list(fixer.topology.chains()))
        chains_to_remove = np.setdiff1d(np.arange(n_chains), self.keep_chains)
        fixer.removeChains(chains_to_remove)
                
        app.PDBFile.writeFile(fixer.topology, fixer.positions, open(out_filename, 'w'))


    def equilibrate(self, ff_name, water_name):
        
        input_pdb_filename = self.get_initial_pdb_filename(ff_name, water_name)
        equilibration_pdb_filename = self.get_equilibration_pdb_filename(ff_name, water_name)
        equilibration_dcd_filename = self.get_equilibration_dcd_filename(ff_name, water_name)
        equilibration_protein_pdb_filename = self.get_equilibration_protein_pdb_filename(ff_name, water_name)
        
        utils.make_path(equilibration_pdb_filename)
        
        if os.path.exists(equilibration_pdb_filename):
            return
        
        ff = app.ForceField('%s.xml' % ff_name, '%s.xml' % water_name)
        pdb = app.PDBFile(get_initial_pdb_filename(ff_name, water_name))
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(ff, model=base_waters[water_name], padding=self.padding)
        topology = modeller.getTopology()
        positions = modeller.getPositions()

        system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=self.cutoff, constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(self.temperature, self.equilibration_friction, self.equilibration_timestep)
        system.addForce(mm.MonteCarloBarostat(self.pressure, self.temperature, self.barostat_frequency))

        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        
        print('Minimizing...')
        simulation.minimizeEnergy()

        simulation.context.setVelocitiesToTemperature(self.temperature)
        print('Running.')
        
        simulation.reporters.append(app.PDBReporter(equilibration_pdb_filename, self.equilibrate_output_frequency))
        simulation.reporters.append(app.DCDReporter(equilibration_dcd_filename, self.equilibrate_output_frequency))
        simulation.step(n_equil_steps)
        del simulation
        del system
        traj = md.load(equilibration_dcd_filename, top=equilibration_pdb_filename)
        traj.save(equilibration_pdb_filename)
        
        top, bonds = traj.top.to_dataframe()
        atom_indices = top.index[top.chainID == 0].values
        traj.restrict_atoms(atom_indices)
        traj.save(equilibration_protein_pdb_filename)
        
        
    def production(self, ff_name, water_name):

        equilibration_pdb_filename = self.get_equilibration_pdb_filename(ff_name, water_name)
        production_dcd_filename = self.get_production_dcd_filename(ff_name, water_name)
        production_protein_dcd_filename = self.get_production_protein_dcd_filename(ff_name, water_name)
        
        utils.make_path(production_pdb_filename)

        if os.path.exists(production_protein_pdb_filename):
            return

        ff = app.ForceField('%s.xml' % ff_name, '%s.xml' % water_name)
        
        traj = md.load(pdb_filename)
        top, bonds = traj.top.to_dataframe()
        atom_indices = top.index[top.chainID == 0].values

        pdb = app.PDBFile(pdb_filename)
        
        system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=self.cutoff, constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(self.temperature, self.friction, self.timestep)
        system.addForce(mm.MonteCarloBarostat(self.pressure, self.temperature, self.barostat_frequency))

        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)

        simulation.context.setVelocitiesToTemperature(self.temperature)
        print('Running.')
        simulation.reporters.append(md.reporters.DCDReporter(dcd_filename, self.output_frequency, atomSubset=atom_indices))
        simulation.reporters.append(app.DCDReporter(dcd_filename_allatoms, self.output_frequency_allatoms))
        simulation.step(self.n_steps)        
        

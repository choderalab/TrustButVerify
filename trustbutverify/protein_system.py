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


class System(Target):
    def __init__(self, temperature, padding=PADDING, cutoff=CUTOFF, keep_chains=None, pressure=PRESSURE, timestep=TIMESTEP, equil_timestep=EQUIL_TIMESTEP, 
        barostat_frequency=BAROSTAT_FREQUENCY, friction=FRICTION, equil_friction=EQUIL_FRICTION, n_steps=N_STEPS, n_equil_steps=N_EQUIL_STEPS, equil_output_frequency=EQUIL_OUTPUT_FREQUENCY, 
        output_frequency=OUTPUT_FREQUENCY, protein_output_frequency=PROTEIN_OUTPUT_FREQUENCY, ionic_strength=0.0 * u.molar, pH=pH):

        self.padding = padding
        self.cutoff = cutoff
        self.pressure = pressure
        self.timestep = timestep
        self.equil_timestep = equil_timestep
        self.friction = friction
        self.equil_friction = equil_friction
        self.barostat_frequency = barostat_frequency
        self.ionic_strength = ionic_strength
        self.pH = pH
        
        self.n_equil_steps = n_equil_steps
        self.n_steps = n_steps
        
        self.equil_output_frequency = equil_output_frequency
        self.output_frequency = output_frequency
        self.protein_output_frequency = protein_output_frequency
        
        self.temperature = temperature


    def equilibrate(self, ff_name, water_name):
        
        input_pdb_filename = self.get_initial_pdb_filename(ff_name, water_name)
        equil_pdb_filename = self.get_equil_pdb_filename(ff_name, water_name)
        equil_dcd_filename = self.get_equil_dcd_filename(ff_name, water_name)
        equil_protein_pdb_filename = self.get_equil_protein_pdb_filename(ff_name, water_name)
        
        utils.make_path(equil_pdb_filename)
        
        if os.path.exists(equil_pdb_filename):
            return
        
        ff = app.ForceField('%s.xml' % ff_name, '%s.xml' % water_name)
        pdb = app.PDBFile(input_pdb_filename)
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(ff, model=water_mapping[water_name], padding=self.padding, ionicStrength=self.ionic_strength)
        topology = modeller.getTopology()
        positions = modeller.getPositions()

        system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=self.cutoff, constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(self.temperature, self.equil_friction, self.equil_timestep)
        system.addForce(mm.MonteCarloBarostat(self.pressure, self.temperature, self.barostat_frequency))

        simulation = app.Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        
        print('Minimizing.')
        simulation.minimizeEnergy()

        simulation.context.setVelocitiesToTemperature(self.temperature)
        print('Equilibrating.')
        
        simulation.reporters.append(app.PDBReporter(equil_pdb_filename, self.n_equil_steps - 1))
        simulation.reporters.append(app.DCDReporter(equil_dcd_filename, self.equil_output_frequency))
        simulation.step(self.n_equil_steps)
        del simulation
        del system
        traj = md.load(equil_dcd_filename, top=equil_pdb_filename)[-1]
        traj.save(equil_pdb_filename)
        
        top, bonds = traj.top.to_dataframe()
        atom_indices = top.index[top.chainID == 0].values
        traj.restrict_atoms(atom_indices)
        traj.save(equil_protein_pdb_filename)
        
        
    def production(self, ff_name, water_name):

        equil_pdb_filename = self.get_equil_pdb_filename(ff_name, water_name)
        production_dcd_filename = self.get_production_dcd_filename(ff_name, water_name)
        production_protein_dcd_filename = self.get_production_protein_dcd_filename(ff_name, water_name)
        
        utils.make_path(production_dcd_filename)

        if os.path.exists(production_protein_dcd_filename):
            return

        ff = app.ForceField('%s.xml' % ff_name, '%s.xml' % water_name)
        
        traj = md.load(equil_pdb_filename)
        top, bonds = traj.top.to_dataframe()
        atom_indices = top.index[top.chainID == 0].values

        pdb = app.PDBFile(equil_pdb_filename)
        
        system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=self.cutoff, constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(self.temperature, self.friction, self.timestep)
        system.addForce(mm.MonteCarloBarostat(self.pressure, self.temperature, self.barostat_frequency))

        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)

        simulation.context.setVelocitiesToTemperature(self.temperature)
        print('Production.')
        simulation.reporters.append(md.reporters.DCDReporter(production_protein_dcd_filename, self.protein_output_frequency, atomSubset=atom_indices))
        simulation.reporters.append(app.DCDReporter(production_dcd_filename, self.output_frequency))
        simulation.step(self.n_steps)
    
    def load(self, ff_name, water_name):
        equil_protein_pdb_filename = self.get_equil_protein_pdb_filename(ff_name, water_name)
        production_protein_dcd_filename = self.get_production_protein_dcd_filename(ff_name, water_name)
        return md.load(production_protein_dcd_filename, top=equil_protein_pdb_filename)
    

class ProteinSystem(System):
    def __init__(self, pdb_id, temperature, pdb_filename=None, keep_chains=None, **kwargs):
        
        super(ProteinSystem, self).__init__(temperature, **kwargs)
        
        self.pdb_id = pdb_id
        self._target_name = pdb_id
        self.pdb_filename = pdb_filename
        
        
        if keep_chains is None:
            keep_chains = np.arange(1)
        
        self.keep_chains = keep_chains


    def build(self, ff_name, water_name):
        out_filename = self.get_initial_pdb_filename(ff_name, water_name)
        utils.make_path(out_filename)

        if os.path.exists(out_filename):
            return
        
        if self.pdb_filename is not None:
            fixer = pdbfixer.PDBFixer(filename=self.pdb_filename)
        else:
            fixer = pdbfixer.PDBFixer(pdbid=self.pdb_id)

        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.removeHeterogens(True)
        fixer.addMissingHydrogens(pH=self.pH)

        n_chains = len(list(fixer.topology.chains()))
        chains_to_remove = np.setdiff1d(np.arange(n_chains), self.keep_chains)
        fixer.removeChains(chains_to_remove)
                
        app.PDBFile.writeFile(fixer.topology, fixer.positions, open(out_filename, 'w'))
    


class PeptideSystem(System):
    def __init__(self, sequence, temperature, N_cap=None, C_cap=None, output_frequency=OUTPUT_FREQUENCY_PEPTIDES, protein_output_frequency=PROTEIN_OUTPUT_FREQUENCY_PEPTIDES, n_steps=N_STEPS_PEPTIDES, **kwargs):
        
        super(PeptideSystem, self).__init__(temperature, output_frequency=output_frequency, protein_output_frequency=protein_output_frequency, n_steps=n_steps, **kwargs)

        self._target_name = "%s_%s_%s" % (N_cap, sequence, C_cap)
        self.sequence = sequence
        self.N_cap = N_cap
        self.C_cap = C_cap


    def build(self, ff_name, water_name):
        out_filename = self.get_initial_pdb_filename(ff_name, water_name)
        utils.make_path(out_filename)

        if os.path.exists(out_filename):
            return

        pdbbuilder.build_pdb(self.sequence, out_filename, self.N_cap, self.C_cap, pH=pH)
        
    

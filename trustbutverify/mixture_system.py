import tempfile
import numpy as np
import os
import mdtraj as md

import pdbfixer
#import pdbbuilder
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

from .target import Target
from .simulation_parameters import *
from .cirpy import resolve 
import utils

from protein_system import System
import gaff2xml
import cStringIO
import itertools

N_STEPS_MIXTURES = 100000 # 0.2ns
N_EQUIL_STEPS_MIXTURES = 100000 # 0.1ns
OUTPUT_FREQUENCY_MIXTURES = 500

class MixtureSystem(System):
    def __init__(self, cas_strings, n_monomers, temperature, pressure=PRESSURE, output_frequency=OUTPUT_FREQUENCY_MIXTURES, n_steps=N_STEPS_MIXTURES, equil_output_frequency=OUTPUT_FREQUENCY_MIXTURES, **kwargs):
        super(MixtureSystem, self).__init__(temperature=temperature, pressure=pressure, output_frequency=output_frequency, n_steps=n_steps, equil_output_frequency=equil_output_frequency, **kwargs)

        self._main_dir = os.getcwd()
        
        self.cas_strings = cas_strings
        self.smiles_strings = []
        for mlc in cas_strings:
            self.smiles_strings.append(resolve(mlc, 'smiles'))

        self.n_monomers = n_monomers
        identifier = list(itertools.chain(cas_strings)) # this is formatted completely incorrectly but will prevent overwrite for now
        self._target_name = '_'.join(identifier)

    def build(self):
        utils.make_path('monomers/')
        utils.make_path('boxes/')
        utils.make_path('ffxml/')
        self.monomer_pdb_filenames = ["monomers/"+string+".pdb" for string in self.cas_strings]
        self.box_pdb_filename = "boxes/" + self.identifier + ".pdb"
        self.ffxml_filename = "ffxml/" + '_'.join(self.cas_strings) + ".xml"

        utils.make_path(self.box_pdb_filename)

        rungaff = False
        if not os.path.exists(self.ffxml_filename):
            rungaff = True
        for filename in self.monomer_pdb_filenames:
            if not os.path.exists(filename):
                rungaff = True

        if rungaff:
            with gaff2xml.utils.enter_temp_directory():  # Avoid dumping 50 antechamber files in local directory.
                ligand_trajectories, self.ffxml = gaff2xml.utils.smiles_to_mdtraj_ffxml(self.smiles_strings)    
            if not os.path.exists(self.ffxml_filename):
                outfile = open(self.ffxml_filename, 'w')
                outfile.write(self.ffxml.read())
                outfile.close()
                self.ffxml.seek(0)

            for k, ligand_traj in enumerate(ligand_trajectories): # will the ligand trajectories always be in the same order as the smiles_strings? yes, right?
                pdb_filename = self.monomer_pdb_filenames[k] # so I can do this?
                ligand_traj.save(pdb_filename)
        else:
            fid = open(self.ffxml_filename,'r')
            self.ffxml = cStringIO.StringIO()
            self.ffxml.write(fid.read())
            fid.close()

        self.packed_trj = gaff2xml.packmol.pack_box(self.monomer_pdb_filenames, self.n_monomers)
        self.packed_trj.save(self.box_pdb_filename)
        
        
    def equilibrate(self):
        utils.make_path('equil/')
        self.equil_dcd_filename = "equil/"+self.identifier +"_equil.dcd"
        self.equil_pdb_filename = "equil/"+self.identifier +"_equil.pdb"
        utils.make_path(self.equil_pdb_filename)
        
        if os.path.exists(self.equil_pdb_filename):
            pass

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
        utils.make_path('production/')
        self.production_dcd_filename = "production/"+self.identifier +"_production.dcd"        

        utils.make_path(self.production_dcd_filename)

        if os.path.exists(self.production_dcd_filename):
            pass
        
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

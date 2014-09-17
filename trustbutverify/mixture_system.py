import tempfile
import numpy as np
import os
import mdtraj as md

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

from .target import Target
from .simulation_parameters import *
from .cirpy import resolve 
import utils
from .nonbondedstatedatareporter import NonbondedStateDataReporter

from protein_system import System
import gaff2xml
import itertools

from pymbar import timeseries as ts
import pandas as pd

N_STEPS_MIXTURES = 500000 # 1 ns (at a time)
N_EQUIL_STEPS_MIXTURES = 5000000 # 5ns
OUTPUT_FREQUENCY_MIXTURES = 5000 # 10ps
OUTPUT_DATA_FREQUENCY_MIXTURES = 125 # 0.25ps
STD_ERROR_TOLERANCE = 0.1 # g/mL

class MixtureSystem(System):
    def __init__(self, cas_strings, n_monomers, temperature, pressure=PRESSURE, output_frequency=OUTPUT_FREQUENCY_MIXTURES, output_data_frequency=OUTPUT_DATA_FREQUENCY_MIXTURES, n_steps=N_STEPS_MIXTURES, equil_output_frequency=OUTPUT_FREQUENCY_MIXTURES, stderr_tolerance = STD_ERROR_TOLERANCE, **kwargs):
        super(MixtureSystem, self).__init__(temperature=temperature, pressure=pressure, output_frequency=output_frequency, n_steps=n_steps, equil_output_frequency=equil_output_frequency, **kwargs)

        self._main_dir = os.getcwd()
        self.cas_strings = cas_strings
        self.n_monomers = n_monomers
        identifier = list(itertools.chain(cas_strings, [str(n) for n in n_monomers], [str(temperature).split(' ')[0]]))
        self._target_name = '_'.join(identifier)
        self.output_data_frequency = output_data_frequency
        self.stderr_tolerance = stderr_tolerance
        self.ran_equilibrate = False

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
            self.smiles_strings = []
            for mlc in self.cas_strings:
                self.smiles_strings.append(resolve(mlc, 'smiles'))
            with gaff2xml.utils.enter_temp_directory():  # Avoid dumping 50 antechamber files in local directory.
                ligand_trajectories, ffxml = gaff2xml.utils.smiles_to_mdtraj_ffxml(self.smiles_strings)    
            if not os.path.exists(self.ffxml_filename):
                outfile = open(self.ffxml_filename, 'w')
                outfile.write(ffxml.read())
                outfile.close()
                ffxml.seek(0)
            for k, ligand_traj in enumerate(ligand_trajectories):
                pdb_filename = self.monomer_pdb_filenames[k]
                if not os.path.exists(pdb_filename):
                    ligand_traj.save(pdb_filename)

        self.ffxml = app.ForceField(self.ffxml_filename)

        if "7732-18-5" in self.cas_strings:
            self.ffxml.loadFile("tip3p.xml")

        if not os.path.exists(self.box_pdb_filename):
            self.packed_trj = gaff2xml.packmol.pack_box(self.monomer_pdb_filenames, self.n_monomers)
            self.packed_trj.save(self.box_pdb_filename)
        else:
            self.packed_trj = md.load(self.box_pdb_filename)

    def equilibrate(self):
        self.ran_equilibrate = True
        utils.make_path('equil/')
        self.equil_dcd_filename = "equil/"+self.identifier +"_equil.dcd"
        self.equil_pdb_filename = "equil/"+self.identifier +"_equil.pdb"
        utils.make_path(self.equil_pdb_filename)
        
        if os.path.exists(self.equil_pdb_filename):
            return

        positions = self.packed_trj.openmm_positions(0)
        topology = self.packed_trj.top.to_openmm()
        topology.setUnitCellDimensions(mm.Vec3(*self.packed_trj.unitcell_lengths[0]) * u.nanometer)
        
        ff = self.ffxml

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
        self.production_pdb_filename = "production/"+self.identifier +"_production.pdb"
        self.production_data_filename = "production/"+self.identifier +"_production.dat"

        utils.make_path(self.production_dcd_filename)

        if os.path.exists(self.production_pdb_filename):
            return        

        if self.ran_equilibrate:
            pdb = app.PDBFile(self.equil_pdb_filename)
            topology = pdb.topology
            positions = pdb.positions
        else:
            positions = self.packed_trj.openmm_positions(0)
            topology = self.packed_trj.top.to_openmm()
            topology.setUnitCellDimensions(mm.Vec3(*self.packed_trj.unitcell_lengths[0]) * u.nanometer)
        
        ff = self.ffxml

        system = ff.createSystem(topology, nonbondedMethod=app.PME, nonbondedCutoff=self.cutoff, constraints=app.HBonds)
        integrator = mm.LangevinIntegrator(self.temperature, self.friction, self.timestep)
        system.addForce(mm.MonteCarloBarostat(self.pressure, self.temperature, self.barostat_frequency))

        simulation = app.Simulation(topology, system, integrator, platform=self.platform)
        simulation.context.setPositions(positions)

        if not self.ran_equilibrate:
            print('Minimizing.')
            simulation.minimizeEnergy()

        simulation.context.setVelocitiesToTemperature(self.temperature)
        print('Production.')
        simulation.reporters.append(app.DCDReporter(self.production_dcd_filename, self.output_frequency))
        simulation.reporters.append(NonbondedStateDataReporter(self.production_data_filename, self.output_data_frequency, step=True, potentialEnergy=True, temperature=True, density=True))

        converged = False
        while not converged:
            simulation.step(self.n_steps)
            d = pd.read_csv(self.production_data_filename, names=["step", "U_NB", "U", "Temperature", "Density"], skiprows=1)
            density_ts = np.array(d.Density)
            [t0, g, Neff] = ts.detectEquilibration(density_ts, nskip=1000)
            density_ts = density_ts[t0:]
            density_mean_stderr = density_ts.std() / np.sqrt(Neff)
            if density_mean_stderr < self.stderr_tolerance:
                converged = True

        del(simulation)
        if self.ran_equilibrate:
            traj = md.load(self.production_dcd_filename, top=self.equil_pdb_filename)[-1]
        else:
            traj = md.load(self.production_dcd_filename, top=self.box_pdb_filename)[-1]
        traj.save(self.production_pdb_filename)

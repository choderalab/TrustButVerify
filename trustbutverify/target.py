import mdtraj as md
import numpy as np
import os
from .simulation_parameters import *


class Target(object):
    
    @property
    def base_path(self):
        return base_path

    @property
    def pdb_path(self):
        return os.path.join(self.base_path, "pdb")
        
    @property
    def equil_path(self):
        return os.path.join(self.base_path, "equil")

    @property
    def production_path(self):
        return os.path.join(self.base_path, "production")
    
    def get_initial_pdb_filename(self, ff_name, water_name):
        return os.path.join(self.pdb_path, self.get_base_filename(ff_name, water_name) + ".pdb")  # Actually, the starting model doesn't have an FF dependence but might as well make it consistent here...
    
    def get_equil_pdb_filename(self, ff_name, water_name):
        return os.path.join(self.equil_path, self.get_base_filename(ff_name, water_name) + ".pdb")

    def get_equil_protein_pdb_filename(self, ff_name, water_name):
        return os.path.join(self.equil_path, self.get_base_filename(ff_name, water_name) + "_protein" + ".pdb")

    def get_production_pdb_filename(self, ff_name, water_name):
        return os.path.join(self.production_path, self.get_base_filename(ff_name, water_name) + ".pdb")

    def get_equil_dcd_filename(self, ff_name, water_name):
        return os.path.join(self.equil_path, self.get_base_filename(ff_name, water_name) + ".dcd")                

    def get_production_dcd_filename(self, ff_name, water_name):
        return os.path.join(self.production_path, self.get_base_filename(ff_name, water_name) + ".dcd")                

    def get_production_protein_dcd_filename(self, ff_name, water_name):
        return os.path.join(self.production_path, self.get_base_filename(ff_name, water_name) + "_protein" + ".dcd")                

    def get_base_filename(self, ff_name=None, water_name=None):
        return "%s_%s_%s" % (ff_name, water_name, self._target_name)
    
    def build(self):
        pass
    
    def production(self):
        pass

    def equilibrate(self):
        pass
        
    def predict(self):
        pass
    
    @property
    def experimental(self):
        pass


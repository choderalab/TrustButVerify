import numpy as np
import pandas as pd
import mdtraj as md
import nmrpystar

class Analyzer(object):
    def __init__(self, identifier, bmrb_filename):
        self.identifier = identifier
        self.bmrb_filename = bmrb_filename    
    
class ChemicalShiftAnalyzer(Analyzer):


    def analyze(self, traj):
        pass
    
    def load_expt(self):
        parsed = nmrpystar.parse(open(self.bmrb_filename).read())
        print(parsed.status)

        q = parsed.value.saves["assigned_chem_shift_list_1"].loops[1]
        x = pd.DataFrame(q.rows, columns=q.keys)
        x = x[["Atom_chem_shift.Seq_ID", "Atom_chem_shift.Atom_ID", "Atom_chem_shift.Val"]]
        x.rename(columns={"Atom_chem_shift.Seq_ID":"resSeq", "Atom_chem_shift.Atom_ID":"name", "Atom_chem_shift.Val":"value"}, inplace=True)

        # Need to make dtypes match to do eventual comparison.
        x["resSeq"] = x["resSeq"].astype('int')
        x["value"] = x["value"].astype('float')
        x["expt"] = "CS"
        x["system"] = self.identifier

        expt = x.set_index(["system", "expt", "resSeq", "name"]).value
        return expt
        
    
class ScalarCouplingAnalyzer(Analyzer):

    def load_expt(self):
        parsed = nmrpystar.parse(open(self.bmrb_filename).read())
        print(parsed.status)

        q = parsed.value.saves["coupling_constant_list_1"].loops[1]
        x = pd.DataFrame(q.rows, columns=q.keys)
        x = x[["Coupling_constant.Seq_ID_1", "Coupling_constant.Val", "Coupling_constant.Val_err"]]
        x.rename(columns={"Coupling_constant.Seq_ID_1":"resSeq", "Coupling_constant.Val":"value", "Coupling_constant.Val_err":"err"}, inplace=True)

        # Need to make dtypes match to do eventual comparison.
        x["resSeq"] = x["resSeq"].astype('int')
        x["value"] = x["value"].astype('float')
        x["expt"] = "3JHNHA"
        x["system"] = self.identifier

        expt = x.set_index(["system", "expt", "resSeq"]).value
        
        return expt

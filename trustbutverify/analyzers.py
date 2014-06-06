import numpy as np
import pandas as pd
import mdtraj as md
import nmrpystar

def multi_index_to_str(multi_index):
    return ["_".join([str(a) for a in ind]) for ind in multi_index]

class Analyzer(object):
    def __init__(self, identifier, bmrb_filename):
        self.identifier = identifier
        self.bmrb_filename = bmrb_filename    
    
class ChemicalShiftAnalyzer(Analyzer):


    def analyze(self, traj):

        top, bonds = traj.top.to_dataframe()
        prediction = md.nmr.chemical_shifts_shiftx2(traj).mean(1).reset_index()  # Average over time dimensions and turn into dataframe
        
        prediction.rename(columns={0:"value"}, inplace=True)  # Give a name to the colum with the actual values.
        prediction["expt"] = "CS"
        prediction["system"] = self.identifier
        
        prediction = prediction.set_index(["system", "expt", "resSeq", "name"]).value
        prediction = pd.Series(prediction.values, multi_index_to_str(prediction.index))
        
        return prediction
    
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
        
        expt = pd.Series(expt.values, multi_index_to_str(expt.index))
        
        return expt
        
    
class ScalarCouplingAnalyzer(Analyzer):
    def analyze(self, traj):
        top, bonds = traj.top.to_dataframe()
        ind, values = md.compute_J3_HN_HA(traj)
        prediction = pd.DataFrame({"value":values.mean(0)})
        prediction["resSeq"] = top.ix[ind[:,-1]].resSeq.values  # Set the residue numbers to the last (fourth) atom in the dihedral
        prediction["expt"] = "3JHNHA"
        prediction["system"] = self.identifier
        
        prediction = prediction.set_index(["system", "expt", "resSeq"]).value
        prediction = pd.Series(prediction.values, multi_index_to_str(prediction.index))
        
        return prediction        
    

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



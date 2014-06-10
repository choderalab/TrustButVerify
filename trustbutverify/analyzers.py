import numpy as np
import pandas as pd
import mdtraj as md
import nmrpystar
from sklearn.externals.joblib import Memory

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]

CACHEDIR = "/home/kyleb/.cached_chemical_shifts/"
memory = Memory(cachedir=CACHEDIR, verbose=0)


def multi_index_to_str(multi_index):
    return ["_".join([str(a) for a in ind]) for ind in multi_index]

class Analyzer(object):
    def __init__(self, identifier, data_filename):
        self.identifier = identifier
        self.data_filename = data_filename    

@memory.cache
def chemical_shift_function(traj, identifier):
        top, bonds = traj.top.to_dataframe()
        prediction = md.nmr.chemical_shifts_shiftx2(traj).mean(1).reset_index()  # Average over time dimensions and turn into dataframe
        
        prediction.rename(columns={0:"value"}, inplace=True)  # Give a name to the colum with the actual values.
        prediction["expt"] = "CS"
        prediction["system"] = identifier
        
        prediction = prediction.set_index(["system", "expt", "resSeq", "name"]).value
        prediction = pd.Series(prediction.values, multi_index_to_str(prediction.index))
        
        return prediction


class ChemicalShiftAnalyzer(Analyzer):

    
    def analyze(self, traj):
        return chemical_shift_function(traj, self.identifier)
    
    def load_expt(self):
        parsed = nmrpystar.parse(open(self.data_filename).read())
        print(parsed.status)
        
        if "assigned_chemical_shifts" in parsed.value.saves:
            q = parsed.value.saves["assigned_chemical_shifts"].loops[1]
        elif "assigned_chem_shift_list_1" in parsed.value.saves:
            q = parsed.value.saves["assigned_chem_shift_list_1"].loops[1]
        else:
            raise(KeyError("Can't find chemical shift assignments in BMRB file %s" % self.data_filename))
        
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
        parsed = nmrpystar.parse(open(self.data_filename).read())
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
        
        expt = pd.Series(expt.values, multi_index_to_str(expt.index))
        
        return expt

class BuzzScalarCouplingAnalyzer(ScalarCouplingAnalyzer):
    def load_expt(self):
        expt = pd.read_csv(self.data_filename, index_col=0)

        aa = self.identifier.split("_")[1]
        expt = expt.ix[[aa]]
        
        expt["value"] = expt["coupling"]
        expt["resSeq"] = 1
        expt["system"] = self.identifier
        expt["expt"] = "3JHNHA"
        
        expt.ix["H"].value = 7.76  # We're using the pH 2.9 result for HIS because that will allow us to simulate fully protonated HIS
        # Rather than need to do a constant pH simulation near the midpoint of the HIS titration curve.
        
        
        expt = expt.set_index(["system", "expt", "resSeq"]).value
        
        expt = pd.Series(expt.values, multi_index_to_str(expt.index))
        return expt


class OhScalarCouplingAnalyzer(ScalarCouplingAnalyzer):
    def load_expt(self):
        #  To DO: FIX HARDCODED PATH!!!
        larger = pd.read_csv("/home/kyleb/src/choderalab/ForcefieldData/nmr/ace_x_y_nh2/data/larger_couplings.csv")
        smaller = pd.read_csv("/home/kyleb/src/choderalab/ForcefieldData/nmr/ace_x_y_nh2/data/smaller_couplings.csv")

        expt = []
        for aa in amino_acids:
            
            value = smaller.ix["G"][aa]
            xyz = ["G%s" % aa, 0, value]
            expt.append(xyz)
            
            value = larger.ix["G"][aa]
            xyz = ["G%s" % aa, 1, value]
            expt.append(xyz)
            
            value = larger.ix[aa]["G"]
            xyz = ["%sG" % aa, 0, value]
            expt.append(xyz)
            
            value = smaller.ix[aa]["G"]
            xyz = ["%sG" % aa, 1, value]
            expt.append(xyz)

        expt = pd.DataFrame(expt, columns=["seq", "resSeq", "value"])
        seq = self.identifier.split("_")[1]
        expt = expt[expt.seq == seq]
        expt["system"] = self.identifier
        expt["expt"] = "3JHNHA"
        expt = expt.set_index(["system", "expt", "resSeq"]).value
        expt = expt.drop_duplicates()
        expt = pd.Series(expt.values, multi_index_to_str(expt.index))
        
        return expt
        

def accumulate_experiments(analyzers_dict):
    data = []
    for key, analyzers in analyzers_dict.items():
        data.append(pd.concat([analyzer.load_expt() for analyzer in analyzers]))
    
    return pd.concat(data).drop_duplicates()

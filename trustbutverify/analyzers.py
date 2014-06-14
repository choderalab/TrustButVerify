import numpy as np
import pandas as pd
import mdtraj as md
import nmrpystar
from sklearn.externals.joblib import Memory
from .simulation_parameters import CS_CACHE_PATH


CHEMICAL_SHIFT_MODEL = "shiftx2"

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]

memory = Memory(cachedir=CS_CACHE_PATH, verbose=0)


def multi_index_to_str(multi_index):
    return ["_".join([str(a) for a in ind]) for ind in multi_index]

class Analyzer(object):
    def __init__(self, identifier, data_filename):
        self.identifier = identifier
        self.data_filename = data_filename    


@memory.cache
def chemical_shift_function(traj, identifier, model):
        prediction = md.compute_chemical_shifts(traj, model=model)
        return prediction


class ChemicalShiftAnalyzer(Analyzer):
    @staticmethod
    def find_assigned_shifts(parsed):
        for key, val in parsed.value.saves.items():
            if "Assigned_chem_shift_list.Sf_category" in val.datums:
                if val.datums["Assigned_chem_shift_list.Sf_category"] == "assigned_chemical_shifts":
                    return val.loops[1]

    @staticmethod
    def old_find_assigned_shifts(parsed):
        if "assigned_chemical_shifts" in parsed.value.saves:
            q = parsed.value.saves["assigned_chemical_shifts"].loops[1]
            print parsed.value.saves["assigned_chemical_shifts"].datums["Assigned_chem_shift_list.Sf_category"]
        elif "assigned_chem_shift_list_1" in parsed.value.saves:
            q = parsed.value.saves["assigned_chem_shift_list_1"].loops[1]
            print parsed.value.saves["assigned_chem_shift_list_1"].datums["Assigned_chem_shift_list.Sf_category"]
        else:
            raise(KeyError("Can't find chemical shift assignments in BMRB file %s" % self.data_filename))        
    
    def analyze(self, traj):
        prediction = chemical_shift_function(traj, self.identifier, CHEMICAL_SHIFT_MODEL).mean(1).reset_index()  # Average over time dimensions and turn into dataframe
        top, bonds = traj.top.to_dataframe()
        prediction.rename(columns={0:"value"}, inplace=True)  # Give a name to the colum with the actual values.
        prediction["expt"] = "CS"
        prediction["system"] = self.identifier


        multi_index = prediction.set_index(["system", "expt", "resSeq", "name"]).index
        prediction["identifier"] = multi_index_to_str(multi_index)
        prediction = prediction.set_index("identifier")
        
        sigma_dict = pd.Series({"N":2.0862, "CA":0.7743, "CB":0.8583, "C":0.8699, "H":0.3783, "HA":0.1967})  # From http://www.shiftx2.ca/performance.html        
        prediction["sigma"] = sigma_dict[prediction.name].values
        
        prediction.rename(columns={"name":"atom"}, inplace=True)  # Use a more descriptive name for the chemical shift atom name

        resSeq_to_AA = top.groupby("resSeq").first().resName
        prediction["AA"] = resSeq_to_AA[prediction.resSeq].values
        
        return prediction
    
    def load_expt(self):
        parsed = nmrpystar.parse(open(self.data_filename).read())
        print(parsed.status)

        q = ChemicalShiftAnalyzer.find_assigned_shifts(parsed)
        
        x = pd.DataFrame(q.rows, columns=q.keys)
        x = x[["Atom_chem_shift.Seq_ID", "Atom_chem_shift.Atom_ID", "Atom_chem_shift.Val"]]
        x.rename(columns={"Atom_chem_shift.Seq_ID":"resSeq", "Atom_chem_shift.Atom_ID":"name", "Atom_chem_shift.Val":"value"}, inplace=True)

        # Need to make dtypes match to do eventual comparison.
        x["resSeq"] = x["resSeq"].astype('int')
        x["value"] = x["value"].astype('float')
        x["expt"] = "CS"
        x["system"] = self.identifier

        expt = x.set_index(["system", "expt", "resSeq", "name"]).value
        
        expt = pd.Series(expt.values, multi_index_to_str(expt.index), name="value")
        
        return expt
        
    
class ScalarCouplingAnalyzer(Analyzer):
    def analyze(self, traj):
        top, bonds = traj.top.to_dataframe()
        ind, values = md.compute_J3_HN_HA(traj)
        prediction = pd.DataFrame({"value":values.mean(0)})

        prediction["resSeq"] = top.ix[ind[:, -1]].resSeq.values  # Set the residue numbers to the last (fourth) atom in the dihedral
        
        if top.ix[0].resName == "ACE":
            prediction["resSeq"] -= 1  # HARDCODED Hack to account for the ACE residue!!!!!!!!!!  Fix me later!
        
        prediction["AA"] = top.ix[ind[:, -1]].resName.values
        prediction["expt"] = "3JHNHA"
        prediction["system"] = self.identifier
        prediction["sigma"] = 0.36
        
        multi_index = prediction.set_index(["system", "expt", "resSeq"]).index
        prediction["identifier"] = multi_index_to_str(multi_index)
        prediction = prediction.set_index("identifier")

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

        expt.ix["H"] = 7.76  # We're using the pH 2.9 result for HIS because that will allow us to simulate fully protonated HIS
        # Rather than need to do a constant pH simulation near the midpoint of the HIS titration curve.        
        
        expt = expt.ix[[aa]]
        
        expt["value"] = expt["coupling"]
        expt["resSeq"] = 1
        expt["system"] = self.identifier
        expt["expt"] = "3JHNHA"

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
            xyz = ["G%s" % aa, 1, value]  # Using indices 1 and 2 here: {Ace:0, X:1, Y:2, NH2:3}
            expt.append(xyz)
            
            value = larger.ix["G"][aa]
            xyz = ["G%s" % aa, 2, value]
            expt.append(xyz)
            
            value = larger.ix[aa]["G"]
            xyz = ["%sG" % aa, 1, value]
            expt.append(xyz)
            
            value = smaller.ix[aa]["G"]
            xyz = ["%sG" % aa, 2, value]
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
    
    return pd.concat(data)

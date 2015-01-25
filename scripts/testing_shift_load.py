import pandas as pd
import nmrpystar

bmrb_filename = "/home/kyleb/src/choderalab/ForcefieldData/nmr/2EVN/6338.str"
parsed = nmrpystar.parse(open(bmrb_filename).read())
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

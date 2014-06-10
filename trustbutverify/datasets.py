from trustbutverify import protein_system, analyzers
from simtk import unit as u
import itertools

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]
amino_acids_noG = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "A", "I", "L", "M", "F", "W", "Y", "V"]

forcefields = ["amber99sbnmr", "amber99sbildn", "amber03", "amber96"]
water_models = ["tip3p", "tip3pfb", "tip4pew", "tip4pfb"]

targets = []
all_analyzers = {}

for aa in amino_acids:    
    all_analyzers["ACE_%s_NME" % aa] = [analyzers.BuzzScalarCouplingAnalyzer("ACE_%s_NME" % aa, "/home/kyleb/src/choderalab/ForcefieldData/nmr/ace_X_NME/experimental_data/baldwin_table1_2006_couplings.csv")]

for aa in amino_acids_noG:
    all_analyzers["ACE_%sG_NH2" % aa] = [analyzers.OhScalarCouplingAnalyzer("ACE_%sG_NH2" % aa, None)]
    all_analyzers["ACE_G%s_NH2" % aa] = [analyzers.OhScalarCouplingAnalyzer("ACE_G%s_NH2" % aa, None)]

all_analyzers["1am7"] = []
all_analyzers["1am7"].append(analyzers.ChemicalShiftAnalyzer("1am7", "/home/kyleb/src/choderalab/ForcefieldData/nmr/bacteriophage_lysozyme/16664.str"))
all_analyzers["1am7"].append(analyzers.ScalarCouplingAnalyzer("1am7", "/home/kyleb/src/choderalab/ForcefieldData/nmr/bacteriophage_lysozyme/19127.str"))

all_analyzers["1d3z"] = [analyzers.ChemicalShiftAnalyzer("1d3z", "/home/kyleb/src/choderalab/ForcefieldData/nmr/ubiquitin/bmrb17439_v3.str")]
all_analyzers["2evn"] = [analyzers.ChemicalShiftAnalyzer("2evn", "/home/kyleb/src/choderalab/ForcefieldData/nmr/2EVN/6338.str")]

for aa in amino_acids:    
    targets.append(protein_system.PeptideSystem(sequence="%s" % aa, N_cap="ACE", C_cap="NME", temperature=303 * u.kelvin, pH=5.0))

for aa in amino_acids_noG:    
    targets.append(protein_system.PeptideSystem(sequence="%s%s" % ("G", aa), N_cap="ACE", C_cap="NH2", temperature=298 * u.kelvin, pH=5.0))
    targets.append(protein_system.PeptideSystem(sequence="%s%s" % (aa, "G"), N_cap="ACE", C_cap="NH2", temperature=298 * u.kelvin, pH=5.0))

targets.append(protein_system.ProteinSystem("1am7", temperature=293 * u.kelvin, pdb_filename="./prepared_pdbs/1AM7_hplusplus.pdb"))  # pH is actually 5.45, but need to use None here to make OpenMM not re-protonate.
targets.append(protein_system.ProteinSystem("1d3z", temperature=298 * u.kelvin, ionic_strength = 0.185 * u.molar))
targets.append(protein_system.ProteinSystem("2evn", temperature=300 * u.kelvin, ionic_strength = 0.05 * u.molar))


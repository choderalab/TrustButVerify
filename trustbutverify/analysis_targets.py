from trustbutverify import analyzers
from simtk import unit as u
import itertools

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]
amino_acids_noG = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "A", "I", "L", "M", "F", "W", "Y", "V"]

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

all_analyzers["1bpi"] = [analyzers.ChemicalShiftAnalyzer("1bpi", "/home/kyleb/src/choderalab/ForcefieldData/nmr/1BPI/4968.str")]

all_analyzers["2lav"] = [analyzers.ChemicalShiftAnalyzer("2lav", "/home/kyleb/src/choderalab/ForcefieldData/nmr/vrk1/16715.str")]

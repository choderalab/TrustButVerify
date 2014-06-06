from trustbutverify import protein_system
from simtk import unit as u
import itertools

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]
amino_acids_noG = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "A", "I", "L", "M", "F", "W", "Y", "V"]

forcefields = ["amber99sbnmr", "amber99sbildn", "amber03", "amber96"]
water_models = ["tip3p", "tip3pfb", "tip4pew", "tip4pfb"]

targets = []

for aa in amino_acids:
    targets.append(protein_system.PeptideSystem(sequence="%s" % aa, N_cap="ACE", C_cap="NME", temperature=303 * u.kelvin, compute_J3HNHA=True))

for aa in amino_acids_noG:
    targets.append(protein_system.PeptideSystem(sequence="%s%s" % ("G", aa), N_cap="ACE", C_cap="NH2", temperature=298 * u.kelvin, compute_J3HNHA=True))
    targets.append(protein_system.PeptideSystem(sequence="%s%s" % (aa, "G"), N_cap="ACE", C_cap="NH2", temperature=298 * u.kelvin, compute_J3HNHA=True))

targets.append(protein_system.ProteinSystem("1am7", temperature=293 * u.kelvin, compute_shifts=True, compute_J3HNHA=True))
targets.append(protein_system.ProteinSystem("1d3z", temperature=298 * u.kelvin, ionic_strength = 0.185 * u.molar, compute_shifts=True))
targets.append(protein_system.ProteinSystem("2evn", temperature=300 * u.kelvin, ionic_strength = 0.05 * u.molar, compute_shifts=True))


from trustbutverify import protein_system
from simtk import unit as u
import itertools

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]

forcefields = ["amber10", "amber03", "amber96"]
water_models = ["tip3p", "tip3pfb", "tip4pew", "tip4pfb"]

temperature = 300 * u.kelvin

targets = []
#targets.append(protein_system.ProteinSystem("1am7", temperature=293 * u.kelvin))
#targets.append(protein_system.ProteinSystem("1d3z", temperature=298 * u.kelvin, ionic_strength = 0.185 * u.molar))
#targets.append(protein_system.ProteinSystem("2evn", temperature=300 * u.kelvin, ionic_strength = 0.05 * u.molar))

for aa in amino_acids:
    targets.append(protein_system.PeptideSystem(sequence="%s" % aa, N_cap="ACE", C_cap="NME", temperature=303 * u.kelvin))

for ff, water, target in itertools.product(forcefields, water_models, targets):
    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)

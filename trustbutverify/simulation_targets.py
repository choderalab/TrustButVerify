from trustbutverify import protein_system
from simtk import unit as u
import itertools

amino_acids = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "A", "I", "L", "M", "F", "W", "Y", "V"]
amino_acids_noG = ["R","H", "K", "D", "E", "S", "T", "N", "Q", "C", "A", "I", "L", "M", "F", "W", "Y", "V"]

forcefields = ["amber99sbnmr", "amber99sbildn", "amber03", "amber96", "a99sb-v2-r1"]
water_models = ["tip3p", "tip3pfb", "tip4pew", "tip4pfb"]

targets = []

for aa in amino_acids:    
    targets.append(protein_system.PeptideSystem(sequence="%s" % aa, N_cap="ACE", C_cap="NME", temperature=303 * u.kelvin, pH=5.0))

for aa in amino_acids_noG:    
    targets.append(protein_system.PeptideSystem(sequence="%s%s" % ("G", aa), N_cap="ACE", C_cap="NH2", temperature=298 * u.kelvin, pH=5.0))
    targets.append(protein_system.PeptideSystem(sequence="%s%s" % (aa, "G"), N_cap="ACE", C_cap="NH2", temperature=298 * u.kelvin, pH=5.0))

targets.append(protein_system.ProteinSystem("1am7", temperature=293 * u.kelvin, pdb_filename="./prepared_pdbs/1AM7_hplusplus.pdb"))  # pH is actually 5.45, but need to use None here to make OpenMM not re-protonate.
targets.append(protein_system.ProteinSystem("1d3z", temperature=298 * u.kelvin, ionic_strength = 0.185 * u.molar, pH=7.2))
targets.append(protein_system.ProteinSystem("2evn", temperature=300 * u.kelvin, ionic_strength = 0.05 * u.molar, pH=6.0))
targets.append(protein_system.ProteinSystem("1bpi", temperature=303 * u.kelvin, ionic_strength = 0.095 * u.molar, pH=5.8))
targets.append(protein_system.ProteinSystem("2lav", temperature=298 * u.kelvin, pdb_filename="./prepared_pdbs/0.15_80_10_pH6.8_2lav.result.pdb", pH=6.8))

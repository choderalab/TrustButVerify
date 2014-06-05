from trustbutverify import protein_system
from simtk import unit as u
import itertools

forcefields = ["amber10", "amber03"]
water_models = ["tip3p", "tip3pfb"]

temperature = 300 * u.kelvin

targets = []
targets.append(protein_system.ProteinSystem("1vii", temperature))
targets.append(protein_system.ProteinSystem("1am7", temperature))
targets.append(protein_system.ProteinSystem("1ubq", temperature))

for ff, water, target in itertools.product(forcefields, water_models, targets):
    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)

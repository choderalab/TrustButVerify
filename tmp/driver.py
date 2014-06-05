from trustbutverify import protein_system
from simtk import unit as u
import itertools

forcefields = ["amber10", "amber03"]
water_models = ["tip3p", "tip3pfb"]
targets = ["1vii", "1am7", "1ubq"]

temperature = 300 * u.kelvin

for ff, water, target in itertools.product(forcefields, water_models, targets):
    target = protein_system.ProteinSystem(target, temperature)
    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)

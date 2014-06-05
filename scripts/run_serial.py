from trustbutverify.datasets import forcefields, water_models, targets
from simtk import unit as u
import itertools

for ff, water, target in itertools.product(forcefields, water_models, targets):
    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)

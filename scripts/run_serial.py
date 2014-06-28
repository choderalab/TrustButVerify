from trustbutverify.simulation_targets import forcefields, water_models, targets
from simtk import unit as u
import itertools

for k, (ff, water, target) in enumerate(itertools.product(forcefields, water_models, targets)):
    print(k, ff, water, target.identifier)
    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)

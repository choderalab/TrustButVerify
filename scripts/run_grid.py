from trustbutverify.simulation_targets import forcefields, water_models, targets
from simtk import unit as u
import itertools
import sys

rank = int(sys.argv[1])

for k, (ff, water, target) in enumerate(itertools.product(forcefields, water_models, targets)):

    if k != rank:
        continue

    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)

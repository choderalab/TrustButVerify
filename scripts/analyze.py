from trustbutverify.datasets import forcefields, water_models, targets
from simtk import unit as u
import itertools

for k, (ff, water, target) in enumerate(itertools.product(forcefields, water_models, targets)):
    t = target.load(ff, water)

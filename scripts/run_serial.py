from trustbutverify.datasets import forcefields, water_models, targets
from trustbutverify.analyzers import accumulate_experiments
from simtk import unit as u
import itertools
import pandas as pd

for ff, water, target in itertools.product(forcefields, water_models, targets):
    print(ff, water, target.identifier)
    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)

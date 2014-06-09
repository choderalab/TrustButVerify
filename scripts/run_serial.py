from trustbutverify.datasets import forcefields, water_models, targets
from trustbutverify.analyzers import accumulate_experiments
from simtk import unit as u
import itertools
import pandas as pd

expt = accumulate_experiments(all_analyzers)

for ff, water, target in itertools.product(forcefields, water_models, targets):
    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)

from trustbutverify.datasets import forcefields, water_models, targets, all_analyzers
from simtk import unit as u
import itertools
import pandas as pd

predictions = []
for ff, water, target in itertools.product(forcefields, water_models, targets):
    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)
    predictions.append(target.analyze(ff, water, all_analyzers))

predictions = pd.concat(predictions)

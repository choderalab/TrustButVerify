from trustbutverify.datasets import forcefields, water_models, targets, all_analyzers
from trustbutverify.analyzers import accumulate_experiments
from simtk import unit as u
import itertools
import pandas as pd

expt = accumulate_experiments(all_analyzers)

predictions = []
for ff, water, target in itertools.product(forcefields, water_models, targets):
    target.build(ff, water)
    target.equilibrate(ff, water)
    target.production(ff, water)
    predictions.append(target.analyze(ff, water, all_analyzers))

predictions = pd.concat(predictions)

predictions["expt"] = expt.ix[predictions.identifier].values
predictions["delta"] = predictions["value"] - predictions["expt"]
predictions["delta2"] = predictions["delta"] ** 2

predictions = predictions.dropna()

rms = (predictions.groupby(["ff", "water"]).delta2).mean() ** 0.5

from trustbutverify.simulation_targets import forcefields, water_models, targets
from trustbutverify.analysis_targets import all_analyzers
from trustbutverify.analyzers import accumulate_experiments
from simtk import unit as u
import itertools
import pandas as pd

expt = accumulate_experiments(all_analyzers)

predictions = []
for ff, water, target in itertools.product(forcefields, water_models, targets):
    print(ff, water, target.identifier)
    predictions.append(target.analyze(ff, water, all_analyzers))

predictions = pd.concat(predictions)

predictions["measured"] = expt.ix[predictions.index].values
predictions["delta"] = predictions["measured"] - predictions["value"]
predictions["z_score"] = predictions["delta"] / predictions["sigma"]
predictions["chi2"] = predictions["z_score"] ** 2

predictions = predictions.dropna()

rms = (predictions.groupby(["ff", "water"]).chi2).mean() ** 0.5
rms_by_system = ((predictions.groupby(["ff", "water", "system"]).chi2).mean() ** 0.5).reset_index().pivot_table(cols=["ff", "water"], rows="system")
rms_by_system = ((predictions.groupby(["ff", "system"]).chi2).mean() ** 0.5).reset_index().pivot_table(cols=["ff"], rows="system")

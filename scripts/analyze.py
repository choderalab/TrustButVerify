from trustbutverify.simulation_targets import forcefields, water_models, targets
from trustbutverify.analysis_targets import all_analyzers
from trustbutverify.analyzers import accumulate_experiments
import trustbutverify.analyzers
from simtk import unit as u
import itertools
import pandas as pd

trustbutverify.analyzers.CHEMICAL_SHIFT_MODEL = "shiftx2"
trustbutverify.analyzers.FIRST_FRAMES = None

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

predictions = predictions.dropna(subset=['chi2'])  # Drop entries with no match between simulation and expt

rms = (predictions.groupby(["ff", "water"]).chi2).mean() ** 0.5
rms_by_ff_only = (predictions.groupby(["ff"]).chi2).mean() ** 0.5
rms_by_system = ((predictions.groupby(["ff", "water", "system"]).chi2).mean() ** 0.5).reset_index().pivot_table(cols=["ff", "water"], rows="system")
rms_by_system_water_averaged = ((predictions.groupby(["ff", "system"]).chi2).mean() ** 0.5).reset_index().pivot_table(cols=["ff"], rows="system")
rms_by_aa = ((predictions.groupby(["ff", "AA"]).chi2).mean() ** 0.5).reset_index().pivot_table(cols=["ff"], rows=["AA"])
rms_by_expt = ((predictions.groupby(["expt", "ff"]).chi2).mean() ** 0.5).reset_index().pivot_table(cols=["ff"], rows="expt")

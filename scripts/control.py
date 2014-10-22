from simtk import unit as u
from trustbutverify import mixture_system

import pandas as pd

import itertools
import sys
import os

data = pd.read_hdf("./data.h5","data")
data = data.set_index("CAS")

for k, (index) in enumerate(data.index):
    cas_strings = [str(index)]
    temperature = data.ix[index].Temperature
    model = mixture_system.MixtureSystem(cas_strings, [1000], temperature)
    model.build()
    model.equilibrate()
    model.production()



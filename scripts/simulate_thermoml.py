import pandas as pd
from simtk import unit as u
from trustbutverify import mixture_system
import sys

data = pd.read_csv("./densities.csv")
rank = int(sys.argv[1])

for k0, k1, components, smiles, cas, temperature, pressure, density in data.itertuples():
    print(k0, k1, components, smiles, cas, temperature, pressure, density)
    model = mixture_system.MixtureSystem([cas], [1000], temperature * u.kelvin, pressure = pressure * u.kilopascal)
    model.build()
    model.equilibrate()
    model.production()

from simtk import unit as u
from trustbutverify import mixture_system

model = mixture_system.MixtureSystem(["CCC", "CCCC"], [10, 20], 300 * u.kelvin, 1.0 * u.atmospheres)
model.build()
model.equilibrate()
model.production()

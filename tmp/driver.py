from trustbutverify import protein_target
from simtk import unit as u

target = protein_target.ProteinTarget("1am7", 300 * u.kelvin)
target.build("amber10", "tip3p")
target.equilibrate("amber10", "tip3p")

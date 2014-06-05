from trustbutverify import protein_target
from simtk import unit as u

ff = "amber10"
water = "tip3p"
pdbid = "1vii"
temperature = 300 * u.kelvin

target = protein_target.ProteinTarget(pdbid, temperature)
target.build(ff, water)
target.equilibrate(ff, water)
target.production(ff, water)

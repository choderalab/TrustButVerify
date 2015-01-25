from trustbutverify import analyzers
import mdtraj as md
import pandas as pd

analyzer = analyzers.ChemicalShiftAnalyzer("1am7", "/home/kyleb/src/choderalab/ForcefieldData/nmr/bacteriophage_lysozyme/16664.str")
expt = analyzer.load_expt()


analyzer = analyzers.ScalarCouplingAnalyzer("1am7", "/home/kyleb/src/choderalab/ForcefieldData/nmr/bacteriophage_lysozyme/19127.str")
expt = analyzer.load_expt()

traj = md.load("/home/kyleb/dat/TrustButVerify/production/amber03_tip3pfb_ACE_AG_NH2.dcd", top="/home/kyleb/dat/TrustButVerify/equil/amber03_tip3pfb_ACE_AG_NH2.pdb")
analyzer.analyze(traj)

from trustbutverify import analyzers

analyzer = analyzers.ChemicalShiftAnalyzer("1am7", "/home/kyleb/src/choderalab/ForcefieldData/nmr/bacteriophage_lysozyme/16664.str")
expt = analyzer.load_expt()


analyzer = analyzers.ScalarCouplingAnalyzer("1am7", "/home/kyleb/src/choderalab/ForcefieldData/nmr/bacteriophage_lysozyme/19127.str")
expt = analyzer.load_expt()

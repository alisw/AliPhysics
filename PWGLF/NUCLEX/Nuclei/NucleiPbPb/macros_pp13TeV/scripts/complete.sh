#!/bin/bash
root -l -b -q Spectra.cc+
root -l -b -q Systematics.cc+
root -l -b -q Systematics.cc+'(true)'
root -l -b -q JoinSystematics.cc+
root -l -b -q Final.cc+
root -l -b -q BWFits.cc+
root -l -b -q BWFits.cc+'(true)'
root -l -b -q Final.cc+
root -l -b -q Final_check.cc+
root -l -b -q Ratio.cc+
root -b -l << EOF
.x myB2.cc+
.L MakeB2FixedPlot.cc+
MakeB2FixedPlot(0,0)
MakeB2FixedPlot(1,0)
MakeB2FixedPlot(2,0)
MakeB2FixedPlot(0,1)
MakeB2FixedPlot(1,1)
MakeB2FixedPlot(2,1)
EOF
root -l -b -q YieldsPlot.cc+
root -l -b -q MakeDsuPplot.cc+
root -l -b -q MakeMeanPt.cc+
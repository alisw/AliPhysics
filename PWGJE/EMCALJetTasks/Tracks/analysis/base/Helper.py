#**************************************************************************
#* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
#*                                                                        *
#* Author: The ALICE Off-line Project.                                    *
#* Contributors are mentioned in the code where appropriate.              *
#*                                                                        *
#* Permission to use, copy, modify and distribute this software and its   *
#* documentation strictly for non-commercial purposes is hereby granted   *
#* without fee, provided that the above copyright notice appears in all   *
#* copies and that both the copyright notice and this permission notice   *
#* appear in the supporting documentation. The authors make no claims     *
#* about the suitability of this software for any purpose. It is          *
#* provided "as is" without express or implied warranty.                  *
#**************************************************************************
"""
Helper tools for spectrum plotting

@author Markus Fasel
"""

from ROOT import TFile, TGraphErrors, gDirectory
from copy import deepcopy

def NormaliseBinWidth(hist):
    """
    Normalise each bin by its width
    """
    for mybin in range(1,hist.GetXaxis().GetNbins()+1):
        bw = hist.GetXaxis().GetBinWidth(mybin)
        hist.SetBinContent(mybin, hist.GetBinContent(mybin)/bw)
        hist.SetBinError(mybin, hist.GetBinError(mybin)/bw)

def MakeRatio(num, den, isBinomial = False):
    """
    Calculate ratio between 2 histograms
    Option indicates whether we use binomial error calculation or gaussian error calculation
    """
    result = deepcopy(num)
    option = ""
    if isBinomial:
        option = "B"
    result.Divide(num, den, 1., 1., option)
    return result

def HistToGraph(hist, xmin = None, xmax = None):
    output = TGraphErrors()
    npoints = 0
    for mybin in range(1, hist.GetXaxis().GetNbins()+1):
        if xmin and hist.GetXaxis().GetBinLowEdge(mybin) < xmin:
            continue
        if xmax and hist.GetXaxis().GetBinLowEdge(mybin) > xmax:
            break
        output.SetPoint(npoints, hist.GetXaxis().GetBinCenter(mybin), hist.GetBinContent(mybin))
        output.SetPointError(npoints, hist.GetXaxis().GetBinWidth(mybin)/2., hist.GetBinError(mybin))
        npoints = npoints + 1
    return output

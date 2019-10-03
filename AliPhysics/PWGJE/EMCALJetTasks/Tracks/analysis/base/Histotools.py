#**************************************************************************
#* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
Module containing simple histogram tools. Currently implemented:
- Rebinner (rebinning only in a subrange)
- Reducing tool (creating histograms in a subrange)
- Combiner tool (combining histograms from different adjacent ranges)

@contact: Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
@organization: ALICE Collaboration
@since: April 14, 2015
"""

from ROOT import TH1D
from numpy import array
from copy import deepcopy

class HistoHandler(object):
  """
  Base class for histogram tools
  inhering classes:
  - Histogram rebinner, rebinning histograms only in a certain range
  - Histogram cutter, making a new histogram in a restricted range
  """
  
  def __init__(self, inputspectrum, xmin, xmax):
    """ 
    Constructor
    @param inputspectrum: input spectrum
    @param xmin: Min. limit of the range
    @param xmax: Max. limit of the range
    """
    self._inputspectrum = deepcopy(inputspectrum)
    self._xrange = {"min":xmin, "max":xmax}
    
  def _GuessBinning(self):
    """
    Obtain binning from the input histogram
    @return: the bin edges of the input histogram
    """
    edges = []
    for i in range(self._inputspectrum.GetXaxis().FindBin(self._xrange["min"]), self._inputspectrum.GetXaxis().FindBin(self._xrange["max"])+1):
      if not len(edges):
          edges.append(self._inputspectrum.GetXaxis().GetBinLowEdge(i))
      edges.append(self._inputspectrum.GetXaxis().GetBinUpEdge(i))
    return edges
  
  def _CreateMiniHist(self, inputbinning):
    """
    Create histogram in restricted range
    @param inputbinning: Binning of the new histogram in python list format
    @return: reduced histogram
    """
    minihist = TH1D("%srestricted" %(self._inputspectrum.GetName()), self._inputspectrum.GetTitle(), len(inputbinning)-1, array(inputbinning))
    for ibin in range(1, minihist.GetXaxis().GetNbins()+1):
      oldbin = self._inputspectrum.GetXaxis().FindBin(minihist.GetXaxis().GetBinCenter(ibin))
      minihist.SetBinContent(ibin, self._inputspectrum.GetBinContent(oldbin))
      minihist.SetBinError(ibin, self._inputspectrum.GetBinError(oldbin))
    return minihist

class HistogramRebinner(HistoHandler):
  """
  Histogram rebinning tool
  """
 
  def __init__(self, inputspectrum, xmin, xmax):
    """
    Constructor
    @param inputspectrum: input spectrum
    @param xmin: Min. limit of the range
    @param xmax: Max. limit of the range
    """
    HistoHandler.__init__(self, inputspectrum, xmin, xmax)
    
  def Rebin(self, rebinfactor):
    """
    Perform rebinning
    @param rebinfactor: Rebinning factor
    @return: reduced rebinned histogram
    """
    return self._CreateMiniHist(self._GuessBinning()).Rebin(rebinfactor)
  
class HistoRangeCutter(HistoHandler):
  """
  Helper class making histogram in a restricted range
  """
  
  def __init__(self, inputspectrum, xmin, xmax): 
    """
    Constructor
    @param inputspectrum: input spectrum
    @param xmin: Min. limit of the range
    @param xmax: Max. limit of the range
    """
    HistoHandler.__init__(self, inputspectrum, xmin, xmax)
    
  def GetReducedHistogram(self):
    """
    Create histogram in restricted area
    @return: reduced histogram
    """
    return self._CreateMiniHist(self._GuessBinning())
   
class HistoCombiner():
  """
  Tool combining histograms in adjacent ranges
  """
  
  def __init__(self):
    """
    Constructor
    """
    self.__components = []
    
  def AddComponent(self, component):
    """
    Add new partial histogram to the sources for the combined histogram
    @param component: Component 
    """
    self.__components.append(component)
    
  def __GetCombinedBinning(self):
    """
    Build combined binning from the different histograms
    using upper edge of a bin to identify bin edges
    @return: Combined binning
    """
    binning = []
    for comp in self.__components:
      for ib in range(1, comp.GetXaxis().GetNbins()+1):
    if not len(binning):
      binning.append(comp.GetXaxis().GetBinLowEdge(ib))
    binning.append(comp.GetXaxis().GetBinUpEdge(ib))
    return binning
  
  def MakeCombinedHistogram(self, name = "combinedHistogram", title = "combined histogram"):
    """
    Build combined histogram of several input histograms
    @param name: Name of the new histogram
    @param title: Title of the new histogram
    @return: The combined histogram
    """
    binning = self.__GetCombinedBinning()
    result = TH1D(name, title, len(binning)-1, array(binning))
    for comp in self.__components:
      for ib in range(1, comp.GetXaxis().GetNbins()+1):
    combinedbin = result.GetXaxis().FindBin(comp.GetXaxis().GetBinCenter(ib))
    result.SetBinContent(combinedbin, comp.GetBinContent(ib))
    result.SetBinError(combinedbin, comp.GetBinError(ib))
    return result

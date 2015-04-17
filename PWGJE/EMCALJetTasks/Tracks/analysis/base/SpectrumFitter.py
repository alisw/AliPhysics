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

from ROOT import TF1, TGraph, TH1F

class FitModel:
    
    def __init__(self):
        self._model = None
        
    def GetFunction(self):
        return self._model
    
    def SetFunctionName(self, name):
        self._model.SetName(name)
        
    def GetFunctionName(self):
        return self._model.GetName()
    
    def SetParLimits(self, parnum, parmin, parmax):
        self._model.SetParLimits(parnum, parmin, parmax)
        
    def SetParameter(self, parnum, parval):
        self._model.SetParameter(parnum, parval)
    
class PowerLawModel(FitModel):
    
    def __init__(self):
        self._model = TF1("fitfunctionPowerlaw", "[0] * TMath::Power(x,[1])", 0., 100.)
    
class ModifiedHagedornModel(FitModel):
    
    def __init__(self):
        self._model = TF1("fitfunctionModHagedorn", "[0]/TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])", 0., 100.)
        self._model.SetParName(0, "A")
        self._model.SetParName(1, "a")
        self._model.SetParName(2, "b")
        self._model.SetParName(3, "p0")
        self._model.SetParName(4, "n")
        # Force all parameters to be positive 
        self._model.SetParLimits(0, 1e-5, 10000)
        self._model.SetParLimits(1, 1e-3, 1)
        self._model.SetParLimits(2, 1e-3, 1)
        self._model.SetParLimits(3, 1e-3, 1)
        self._model.SetParLimits(4, 1e-5, 10)
        self._model.SetParameter(0, 100)
        self._model.SetParameter(1, 0.5)
        self._model.SetParameter(2, 0.5)
        self._model.SetParameter(3, 0.5)
        self._model.SetParameter(4, 5.)

class SpectrumFitter:
    '''
    Class fitting a raw spectrum
    '''
    
    class SpectrumFitterException(Exception):
        
        def __init__(self):
            pass
        
        def __str_(self):
            return "Fit of the spectrum not yet performed"

    def __init__(self, name, spectrum, fitmodel = None):
        '''
        Constructor
        '''
        self.__name = name
        self.__data = spectrum
        self.__model = fitmodel
        if not self.__model:
            self.__model = PowerLawModel()
        self.__fitDone = False
        
    def SetFitModel(self, fitmodel):
        self.__model = fitmodel
        
    def DoFit(self, rangemin, rangemax = 50):
        self.__data.Fit(self.__model.GetFunction(), "N", "", rangemin, rangemax)
        self.__fitDone = True
        
    def GetParameterisation(self):
        if not self.__fitDone:
            raise self.SpectrumFitterException()
        return self.__model.GetFunction()
    
    def GetParameterisedValueAt(self, x):
        if not self.__fitDone:
            raise self.SpectrumFitterException()
        return self.__model.GetFunction().Eval(x)
   
    def CalculateIntegralAbove(self, x):
        if not self.__fitDone:
            raise self.SpectrumFitterException()
        return self.__model.GetFunction().Integral(x, 10000)
   
    def CalculateNormalisedIntegralAbove(self, x):
        """
        Calculate per-event yield above a certain pt, done as sum in 1 GeV bins of the 
        binned content from a min. pt to a max. pt.
        """
        if not self.__fitDone:
            raise self.SpectrumFitterException()
        maxint = 1000.
        ptstart = x
        yieldval = 0
        while ptstart < maxint:
            yieldval += self.__model.GetFunction().Integral(ptstart, ptstart+10)/10.
            ptstart += 10
        return yieldval
            
    def MakeBinnedParameterisation(self, nbins, xmin, xmax, normBinWidth = False):
        """
        Make binned parameterisation. If normBinWith is set to True, the integral is 
        normalised by the bin width
        """
        result = TH1F("binned%s" %(self.__name), "", nbins, xmin, xmax)
        for mybin in range(2, result.GetXaxis().GetNbins()+1):
            value = 0
            if normBinWidth:
                value = self.CalculateBinMean(result.GetXaxis().GetBinLowEdge(mybin),result.GetXaxis().GetBinUpEdge(mybin))
            else:
                value = self.CalculateIntegral(result.GetXaxis().GetBinLowEdge(mybin),result.GetXaxis().GetBinUpEdge(mybin))
            #print "min %f, max %f, value %e" %(result.GetXaxis().GetBinLowEdge(mybin),result.GetXaxis().GetBinUpEdge(mybin), value)
            result.SetBinContent(mybin, value)
        return result
    
    def MakeBinnedParameterisationDefault(self, normBinWidth = False):
        """
        Make binned parameterisation. If normBinWith is set to True, the integral is 
        normalised by the bin width
        """
        result = TH1F("binned%s" %(self.__name), "", self.__data.GetXaxis().GetNbins(), self.__data.GetXaxis().GetXbins().GetArray())
        for mybin in range(2, result.GetXaxis().GetNbins()+1):
            value = 0
            if normBinWidth:
                value = self.CalculateBinMean(result.GetXaxis().GetBinLowEdge(mybin),result.GetXaxis().GetBinUpEdge(mybin))
            else:
                value = self.CalculateIntegral(result.GetXaxis().GetBinLowEdge(mybin),result.GetXaxis().GetBinUpEdge(mybin))
            #print "min %f, max %f, value %e" %(result.GetXaxis().GetBinLowEdge(mybin),result.GetXaxis().GetBinUpEdge(mybin), value)
            result.SetBinContent(mybin, value)
            result.SetBinError(mybin, 0)
        return result
            
    def CalculateIntegral(self, xmin, xmax):
        if not self.__fitDone:
            raise self.SpectrumFitterException()
        return self.__model.GetFunction().Integral(xmin, xmax)
    
    def CalculateBinMean(self, xmin, xmax):
        if not self.__fitDone:
            raise self.SpectrumFitterException()
        return self.__model.GetFunction().Integral(xmin, xmax)/(xmax - xmin)
                
    def GetFitFunction(self):
        return self.__model.GetFunction()
        
class MinBiasFitter(SpectrumFitter):
    
    def __init__(self, name, data, fitmin = 15., model = None):
        SpectrumFitter.__init__(self, name, data)
        self.DoFit(fitmin, 50.)
        
class TriggeredSpectrumFitter(SpectrumFitter):
    
    def __init__(self, name, data, model = None):
        SpectrumFitter.__init__(self, name, data)
        self.DoFit(50., 90.)

        
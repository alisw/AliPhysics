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
Toolset module for histogram operations

Original author (ROOT macro):
@author: Jacek Otwinowski
@organization: ALICE Collaboration
Translated into PYTHON by  
@author: Markus Fasel
@contact: markus.fasel@cern.ch
@organization: ALICE Collaboration
@organization: Lawrence Berkeley National Laboratory
@copyright: 1998-2014, ALICE Experiment at CERN, All rights reserved.
"""

from ROOT import TF1, TGraph, TGraphAsymmErrors, TGraphErrors, TMultiGraph
from PWGJE.EMCALJetTasks.Tracks.analysis.util.Interpolator import Interpolator
import math
from copy import deepcopy
from numpy import array
from scipy.optimize import fsolve 
import functools
import operator

        
class SpectrumTools(object):
    """
    Toolset for histogram operations 
    """


    def __init__(self):
        """
        Constructor
        """
        pass 
   
        
    def DivideGraphs(self, g1, g2, interpolationMethod = "lin"):
        """
        Divide 2 graphs g1 and g2, by each other applying an interpolation between the points for the denominator
        @param g1: graph 1
        @param g2: graph 2
        @param interpolationMethod: method used for interpolation between the two graphs 
        @return: The ratio graphs
        """
    
        #calculate x range
        xmin = max(g1.GetX()[0], g2.GetX()[0])
        xmax = min(g1.GetX()[g1.GetN()-1], g2.GetX()[g2.GetN()-1])

        
        # create new TGraphErors for result1
        result1 = self.MakeGraphDivision(g1, g2, xmin, xmax, interpolationMethod)
        result1.SetMarkerColor(1);
        result1.SetMarkerStyle(20);

        # create new TGraphErors for result2
        result2 = self.MakeGraphDivision(g2, g1, xmin, xmax, interpolationMethod)
        result2.SetMarkerColor(2)
        result2.SetMarkerStyle(20)

        result = TMultiGraph("result"," ")
        result.Add(result1)
        result.Add(result2)
        result.Draw("AP")
        result.GetXaxis().SetRangeUser(xmin,xmax)
        result.GetXaxis().SetTitle("p_{T} [Gev/c]")
        result.GetXaxis().SetTitleOffset(1.2)
        result.GetYaxis().SetTitleOffset(1.2)

        return result
    
    def MakeGraphDivision(self, numerator, denominator, xmin, xmax, interpolationMethod = "lin"):
        """
        Divide 2 graphs by each other.
        For the denominator graph use interpolation
        @param numerator: Numerator graph
        @param denominator: Denominator graph
        @param xmin: lower boundary for calculation
        @param xmax: upper boundary for calculation
        @param interpolationMethod: method used for the interpolation   
        @return: the ratio of the graphs
        """
        startBin = numerator.GetN()
        endBin = -1
        nBins = 0

        for i in range(0, numerator.GetN()+1):
            x = numerator.GetX()[i]
            if x >= xmin and x <= xmax: 
                nBins += 1 
                if i > endBin:
                    endBin = i
                if i < startBin:
                    startBin = i

        print "startBin %d" %(startBin)
        print "endBin %d" %(endBin) 
        print "nBins %d" %(nBins)
        print "xmin %f" %(xmin)
        print "xmax %f" %(xmax)
        
        # create new TGraphErors for result1
        result = TGraphErrors(nBins);

        interpolator = Interpolator()
        for i in range(startBin, endBin+1):
            x = numerator.GetX()[i]
            lower, upper = self.__FindNeighbors(x, denominator)
            # after break the two neighboring bins are found: k-1 and k
            y1  = denominator.GetY()[lower]
            y2  = denominator.GetY()[upper]
            ey1 = denominator.GetEY()[lower] 
            ey2 = denominator.GetEY()[upper] 
            if x == lower: 
                y  = y1
                ey = ey1
            elif x == upper: 
                y  = y2
                ey = ey2
            else:
                y  = interpolator.Interpolate(x,lower,y1,upper,y2, method = interpolationMethod)
                ey = interpolator.Interpolate(x,lower,ey1,upper,ey2, method=interpolationMethod)
            numy = numerator.GetY()[i]
            numey = numerator.GetEY()[i]
            yr = 0
            eyr = 0
            if y == 0:
                yr = 0
            else:
                yr = numy / y  
            if y != 0  and numy != 0:
                eyr = (yr * math.sqrt((ey/y)*(ey/y) + (numey/numy)*(numey/numy)))
            else:
                eyr = 0
            result.SetPoint(i-startBin,x,yr);
            result.SetPointError(i-startBin,numerator.GetEX()[i],eyr)
        return result

    
    def __FindNeighbors(self, x, graph):
        """
        Find the upper and lower neighbor of point x in the graph
        @param x: point to evaluate
        @param graph: graph to check  
        """
        lower = -1
        upper = -1
        for i in range(0, graph.GetN()-1):
            if x >= graph.GetX()[i] and x <= graph.GetX()[i+1]:
                lower = i 
                upper = i+1
                break
        return lower, upper
        
    def RebinPtSpectrum(self, h, nBins = 0, xBins = None):
        """
        Apply rebinning of the spectrum
        @param h: Input histogram
        @param nBins: Number of bins
        @param xbins: Binning of the new histogram
        @return: The rebinned histogram
        """
        if not h:
            return None
        if not nBins:
            return h
        if not xBins:
            return h
        h1 = deepcopy(h)
        h1.Sumw2()
        for i in range(1, h1.GetNbinsX()+1):
            value = h1.GetBinContent(i)
            width = h1.GetBinWidth(i)
            center = h1.GetBinCenter(i)
            error = h1.GetBinError(i)
            h1.SetBinContent(i,value*center*width)
            h1.SetBinError(i,error*center*width)

        h2 = h1.Rebin(nBins,"hnew",array(xBins))   
        for i in range(1, h2.GetNbinsX()+1):
            value = h2.GetBinContent(i)
            width = h2.GetBinWidth(i)
            center = h2.GetBinCenter(i)
            error = h2.GetBinError(i)
            h2.SetBinContent(i,value/(center*width))
            h2.SetBinError(i,error/(center*width))
        return h2

    def ApplyBinShiftCorrection(self, hist):
        """
        Apply bin-shift correction to the input spectrum using an iterative procedure
        @param hist: Input spectrum
        @return: Bin-shift corrected spectrum 
        """
    
        h = deepcopy(hist)
        h.SetName("htemp")    
    
        # Bin shift correction performed in model specturm * pt
        for i in range(1, h.GetNbinsX()+1):
            pt = h.GetBinCenter(i)
            h.SetBinContent(i, h.GetBinContent(i)*pt)
            h.SetBinError(i, h.GetBinError(i)*pt)
   
        result = TGraphErrors(h)
        for i in range(0, result.GetN()): 
            result.GetEX()[i] = 0.   
   
        fitfun = TF1("fitfun","([0]*(1.+x/[1])^(-[2])*x)-[3]",0.15,100.0)
        fitfun.SetParameter(0,1000)
        fitfun.SetParameter(1,1)
        fitfun.SetParameter(2,5)
        fitfun.FixParameter(3,0)
        h.Fit(fitfun,"") 
        self.__StableFit(h, fitfun, True)
   
        # Iterative approach:
        # - Use model to get the mean of the function inside the bin
        # - Get the X where the mean is found
        # - Use the new coordinate (x,y) for the next iteration of the fit
        # for now 10 iterations fixed
        for k in range(1, 11):
            for i in range(1, h.GetNbinsX()+1):
                y = fitfun.Integral(h.GetBinLowEdge(i), h.GetBinUpEdge(i)) / h.GetBinWidth(i)
                result.GetX()[i-1] = self.FindX(y, fitfun, h.GetBinLowEdge(i), h.GetBinUpEdge(i))
            self.__StableFit(result, fitfun, False)
   
        # Undo multiplication with pt
        for i in range(0, result.GetN()):
            pt = result.GetX()[i]
            result.GetY()[i]  /=  pt
            result.GetEY()[i] /=  pt
   
        #remove points that are 0
        while result.GetY()[0] < 1.e-99: 
            result.RemovePoint(0) 
   
        bval = 0
        for mybin in range(0, result.GetN()+1):
            if result.GetY()[bin] < 1.e-99: 
                bval = mybin
                break
            
        while result.RemovePoint(bval) > 0: 
            continue 
        return result
    
    def ApplyBinShiftCorrectionGeneral(self, hist, fit):
        """
        Alternative method for bin shift correction:
        - Apply user-default model for bin-shift correction
        - don't multiply by pt
        @param hist: Input spectrum for the bin shift correction
        @param fit: Model for the bin-shift correction
        @return: The bin-shift corrected spectrum as graph
        """
        h = deepcopy(hist)    
        hist.SetName("htemp")
   
        result = TGraphErrors(h);
        for i in range(0, result.GetN()):
            result.GetEX()[i] = 0.   
        y = 0
   
        #for now 10 iterations fixes
        for k in range(0, 10):
            for i in range(1, h.GetNbinsX()+1):
                y = fit.Integral(h.GetBinLowEdge(i),h.GetBinUpEdge(i)) / h.GetBinWidth(i)
                x = self.FindX(y, fit, h.GetBinLowEdge(i), h.GetBinUpEdge(i)) 
                result.GetX()[i-1] = x
  
        # remove points that are 0
        while result.GetY()[0] < 1.e-99:
            result.RemovePoint(0) 
   
        mybin = 0
        for biniter in range(0, result.GetN()):
            if result.GetY()[biniter] < 1.e-99:
                mybin = biniter
                break
        while result.RemovePoint(mybin) > 0:
            continue 
      
        return result

    def __StableFit(self, inputdata, model, doIntegral = False):
        """
        Perform stable fit: Fit until parameters don't change anymore with iteration
        @param inputdata: Input data to fit (TGraph or TH)
        @param model: Fit model
        @param doIntegral: if true we perform integration during fit
        """
        last = 0
        while True:
            inputdata.Fit(model,"IMR" if doIntegral else "MR")
            params = []
            for ipar in range(0, model.GetNumberOfParameters()):
                params.append(model.GetParameter(ipar))
            current = functools.reduce(operator.mul, params)    
            if current == last:
                break
            last= current
    
    def FindX(self, y, function, xmin, xmax):
        """
        ROOT-finding in PYTHON style
        @param y: y-value
        @param function: input of the equation
        @param xmin: min. x (for initial guess)
        @param xmax: max. x (for initial guess)
        @return: Solution of the equation f(x) = y    
        """
        return fsolve(lambda x : function.Eval(x)-y, (xmax - xmin)/2.)
   

    def PerformSpectrumNormalization(self, h, nevents, etarange):
        """
        Perform normalization
        1/(2 pi pt deleatEta Nevents)
        @param h: Spectrum histogram
        @param nevents: Number of events to scale
        @param etarange: Delta eta 
        @deprecated: Use Normalization class in the correction package
        """
        for i in range(1, h.GetNbinsX()+1):
            pt = h.GetBinCenter(i)
            width = h.GetBinWidth(i)
            val = h.GetBinContent(i)
            err = h.GetBinError(i) 
            cval = 0
            cerr = 0
            if  (etarange*nevents>0):
                cval = (val)/(width * 2.0 * math.pi * etarange * nevents * pt)
                cerr = (err)/(width * 2.0 * math.pi * etarange * nevents * pt)
            h.SetBinContent(i,cval);
            h.SetBinError(i,cerr);

    def MakeHisrogramFromGraph(self, g, prototype, options = "lin"):
        """
        Create a histogram from a TGraphErrors:
        options can be "lin" (default), "log", "exp", "pow" for the functional shape
        and "I" for intergal
        if g is a TGraphErrors errors will be calculated assuming uncorrelated errors
        if option "E" is given errors are calculated
        if option "EX" is given errors only from X errors
        "EY" is give only from Y
        "E0" supresses calculation of errors
        "EC" calculates correlated errors (for systematics)
        "I2" uses integral for value, but eval for errors
        @param g: inputgraph
        @param prototype: prototype histogram
        @param option: see above
        @return: Histogram created from the graph 
        """
        errx = False                # flag for x errors
        erry = False                # flag for y errors
        errc = False                # flag for syst errors  
        if isinstance(g, TGraphErrors):
            errx = True 
            erry = True
        
        if "E" in options:
            errx = True 
            erry = True
        if "EX" in options: 
            erry = False
        if "EY" in options:
            errx = False  
        if "EC" in options:
            errx = False 
            erry = True 
            errc = True
        if "E0" in options:
            errx = False
            erry = False
            errc = False
        if isinstance(g, TGraphErrors) and not g.GetEX(): 
            errx = False
        if isinstance(g, TGraphErrors) and not g.GetEY(): 
            erry = False
    
        xmin = g.GetX()[0]
        xmax = g.GetX()[g.GetN()-1]
        h = deepcopy(prototype)
        h.Reset()
  
        dx1 = 0
        dx2 = 0 
        ey = 0
        for i in range(1, h.GetNbinsX()+1):
            x = h.GetBinCenter(i)   
            # check the x range
            if x < xmin: 
                continue
            if x > xmax: 
                break
            # find point k in g closest in x
            lower, upper = self.__FindNeighbors(x, g)
            # now x1 and x2 are the points next to x
            x1 = g.GetX()[lower]
            x2 = g.GetX()[upper]
            y1 = g.GetY()[lower]
            y2 = g.GetY()[upper]
            y  = self.__GetInterpolatedValue(x,x1,y1,x2,y2,options,h.GetXaxis().GetBinLowEdge(i),h.GetXaxis().GetBinUpEdge(i))    
            if errx:
                exlow, exhigh = self.__GetXerrors(g, lower)
                dx1 = max(exlow, exhigh) if exlow and exhigh else 0.
                exlow, exhigh = self.__GetXerrors(g, upper)
                dx2 = max(exlow, exhigh) if exlow and exhigh else 0.
            if erry:
                eylow, eyhigh = self.__GetYerrors(g, lower)
                dy1 = max(eylow, eyhigh) if eylow and eyhigh else 0.
                eylow, eyhigh = self.__GetYerrors(g, upper)
                dy2 = max(eylow, eyhigh) if eylow and eyhigh else 0.
            if errx or erry:
                if errc:
                    ymax = self.__GetInterpolatedValue(x,x1,y1+dy1,x2,y2+dy2,options,h.GetXaxis().GetBinLowEdge(i),h.GetXaxis().GetBinUpEdge(i))
                    ymin = self.__GetInterpolatedValue(x,x1,y1-dy1,x2,y2-dy2,options,h.GetXaxis().GetBinLowEdge(i),h.GetXaxis().GetBinUpEdge(i))
                    ey = max(math.fabs(y-ymin),math.fabs(y-ymax))
                else:
                    ey = self.__GetInterpolatedUncertainty(x,x1,y1,x2,y2,dx1,dy1,dx2,dy2,options,h.GetXaxis().GetBinLowEdge(i),h.GetXaxis().GetBinUpEdge(i))
            h.SetBinContent(i,y)
            h.SetBinError(i,ey)
        h.SetName(g.GetName()) 
        return h
    
    def __GetXerrors(self, inputgraph, pointID):
        """
        Get the x-errors (low, up) in a transparent way for TGraph and TGraphErrors
        @param inputgraph: input for the graph
        @return: tuple of lower and upper x error
        """
        if isinstance(inputgraph, TGraphAsymmErrors):
            return inputgraph.GetEXlow()[pointID], inputgraph.GetEXhigh()[pointID]
        elif isinstance(inputgraph, TGraphErrors):
            return inputgraph.GetEX()[pointID], inputgraph.GetEX()[pointID]
        else:
            return None
        
    def __GetYerrors(self, inputgraph, pointID):
        """
        Get the y-errors (low, up) in a transparent way for TGraph and TGraphErrors
        @param inputgraph: input for the graph
        @return: tuple of lower and upper y error
        """
        if isinstance(inputgraph, TGraphAsymmErrors):
            return inputgraph.GetEYlow()[pointID], inputgraph.GetEYhigh()[pointID]
        elif isinstance(inputgraph, TGraphErrors):
            return inputgraph.GetEY()[pointID], inputgraph.GetEY()[pointID]
        else:
            return None

    def SetParameters(self, f, x1, y1, x2, y2):
        """
        Set Function parameters
        @param x1: lower x value
        @param y1: y-value at x1 
        @param x2: lower x value
        @param y2: y-value at x1 
        """
        if f.GetName() == "lin":
            f.SetParameter(1,(y1-y2)/(x1-x2))           # slope
            f.SetParameter(0,(x1*y2-x2*y1)/(x1-x2))     # const
        elif f.GetName() == "exp":
            f.SetParameter(1, math.log(y1/y2)/(x1-x2))        # exp
            f.SetParameter(0,y1*math.pow(y1/y2,x1/(x2-x1)))  # constfactor
        elif f.GetName() == "pow":        
            f.SetParameter(1,math.log(y1/y2)/math.log(x1/x2))     # power
            f.SetParameter(0,y1*math.pow(x1,-math.log(y1/y2)/math.log(x1/x2))); #constfactor
        else:
            gtemp = TGraph(2)
            gtemp.SetPoint(0,x1,y1)
            gtemp.SetPoint(1,x2,y2)
            gtemp.Fit(f,"")


    def MakeHistogramFromGraphSimple(self, g, prototype):
        """
        Create a histogram using TGraph's 
        interpolation between points
        @param g: input graph
        @param prototype: prototype for the 
        @return: the new histogram
        """
        h = deepcopy(prototype)
        h.SetName(g.GetName()) 
        h.Reset();
        for i in range(1, g.GetNbinsX()+1):
            h.SetBinContent(i, g.Eval(h.GetBinCenter(i)))
            h.SetBinError(i,0)
        return h

    def __GetInterpolatedValue(self, x, x1, y1, x2, y2, options = "lin", xmin = 0, xmax = 0):
        """
        Get value at x, interpolated using the points (x1,y1) and (x2,y2) as steps for the interpolation.
        Several models can be applied:
        - linear
        - exponential
        - power law
        @param x: Step at which to perform the interpolation
        @param x1: lower x coordinate
        @param y1: Function value at x1
        @param x2: upper x coordinate
        @param y2: Function value at x2
        @param option: Interpolatoin method
        @param xmin:
        @param xmax:
        @return: interpolation value    
        """
        integrate = True if "I" in options else False
        if "lin" in options:
            if integrate: 
                return (-2*x2*y1 + (xmax + xmin)*(y1 - y2) + 2*x1*y2)/(2.*(x1 - x2))
            else:
                return (x*y1 - x2*y1 - x*y2 + x1*y2)/(x1 - x2)
            # end of lin
        elif "exp" in options:
            if integrate: 
                return ((x1 - x2)*y1*(math.pow(y1/y2,xmax/(x1 - x2)) - math.pow(y1/y2,xmin/(x1 - x2)))* math.pow(y1/y2,x1/(-x1 + x2)))/((xmax - xmin)*math.log(y1/y2));
            else:
                return y1*math.pow(y1/y2,(x - x1)/(x1 - x2));
            #end of "exp"
        elif "pow" in options:
            c = math.pow(x1,-math.log(y1/y2)/(math.log(x1) - math.log(x2)))
            n = math.log(y1/y2)/(math.log(x1) - math.log(x2))
            if integrate: 
                if math.fabs(n+1.) < 1e-6:
                    return (c*(math.log(xmax) - math.log(xmin)))/(xmax - xmin)
                else: 
                    (c*(math.pow(xmax,1 + n) - math.pow(xmin,1 + n)))/((1 + n)*(xmax - xmin))
            else:
                return c*math.pow(x,n)
            # end of "pow"
        else:
            return 0

    def __GetInterpolatedUncertainty(self, x, x1, y1, x2, y2, dx1, dy1, dx2, dy2, options = "lin", xmin = 0, xmax = 0):
        """
        Get error at x, interpolated using the points (x1,y1) and (x2,y2) as steps for the interpolation.
        Several models can be applied:
        - linear
        - exponential
        - power law
        @param x: Step at which to perform the interpolation
        @param x1: lower x coordinate
        @param y1: Function value at x1
        @param x2: upper x coordinate
        @param y2: Function value at x2
        @param dx1: uncertainty in x at point 1 
        @param dy1: uncertainty in y at point 1 
        @param dx2: uncertainty in x at point 2 
        @param dy2: uncertainty in y at point 2 
        @param option: Interpolation method
        @param xmin:
        @param xmax:
        @return: error value
        """  
        integrate = True if "I" in options or "I2" in options else False 

        if "lin" in options:
            if integrate: 
                return math.sqrt((math.pow(dy2,2)*math.pow(x1 - x2,2)*math.pow(-2*x1 + xmax + xmin,2) +  
                                  math.pow(dy1,2)*math.pow(x1 - x2,2)*math.pow(-2*x2 + xmax + xmin,2) + 
                                  (math.pow(dx2,2)*math.pow(-2*x1 + xmax + xmin,2) + 
                                   math.pow(dx1,2)*math.pow(-2*x2 + xmax + xmin,2))*math.pow(y1 - y2,2))/math.pow(x1 - x2,4))/2.;
            else:
                return math.sqrt((math.pow(dy2,2)*math.pow(x - x1,2)*math.pow(x1 - x2,2) + \
                                  math.pow(dy1,2)*math.pow(x - x2,2)*math.pow(x1 - x2,2) + \
                                  (math.pow(dx2,2)*math.pow(x - x1,2) + math.pow(dx1,2)*math.pow(x - x2,2))*math.pow(y1 - y2,2))/math.pow(x1 - x2,4));
    
            # end of "lin"
        elif "exp" in options:
            if integrate: 
                return math.sqrt((math.pow(y1/y2,(2*x1)/(-x1 + x2))*(math.pow(dy1*math.pow(y1/y2,xmax/(x1 - x2))* 
                            (x1 - x2 + (x2 - xmax)*math.log(y1/y2)) - dy1*math.pow(y1/y2,xmin/(x1 - x2))* 
                            (x1 - x2 + (x2 - xmin)*math.log(y1/y2)),2) + (math.pow(dy2,2)*math.pow(y1,2)* 
                            math.pow(math.pow(y1/y2,xmax/(x1 - x2))*(x1 - x2 + (x1 - xmax)*math.log(y1/y2)) + 
                            math.pow(y1/y2,xmin/(x1 - x2))*(-x1 + x2 + (-x1 + xmin)*math.log(y1/y2)),2))/math.pow(y2,2) + 
                            (math.pow(dx2,2)*math.pow(y1,2)*math.pow(math.log(y1/y2),2)*math.pow(math.pow(y1/y2,xmax/(x1 - x2))* 
                            (x1 - x2 + (x1 - xmax)*math.log(y1/y2)) + math.pow(y1/y2,xmin/(x1 - x2))*(-x1 + x2 + (-x1 + xmin)* 
                            math.log(y1/y2)),2))/math.pow(x1 - x2,2) + (math.pow(dx1,2)*math.pow(y1,2)*math.pow(math.log(y1/y2),2)* 
                            math.pow(math.pow(y1/y2,xmax/(x1 - x2))*(x1 - x2 + (x2 - xmax)*math.log(y1/y2)) + math.pow(y1/y2,xmin/(x1 - x2))* 
                            (-x1 + x2 + (-x2 + xmin)*math.log(y1/y2)),2))/math.pow(x1 - x2,2)))/(math.pow(xmax - xmin,2)*math.pow(math.log(y1/y2),4)))
            else:
                return math.sqrt((math.pow(y1/y2,(2*(x - x1))/(x1 - x2))*((math.pow(x1 - x2,2)*(math.pow(dy2,2)*math.pow(x - x1,2)* 
                                math.pow(y1,2) + math.pow(dy1,2)*math.pow(x - x2,2)*math.pow(y2,2)))/math.pow(y2,2) + (math.pow(dx2,2)*math.pow(x - x1,2)  
                                + math.pow(dx1,2)*math.pow(x - x2,2))*math.pow(y1,2)*math.pow(math.log(y1/y2),2)))/math.pow(x1 - x2,4))
            # end of "exp"
        elif "pow" in options:
            c = math.pow(x1,-math.log(y1/y2)/(math.log(x1) - math.log(x2)))            
            n = math.log(y1/y2)/(math.log(x1) - math.log(x2))
            if integrate: 
                if math.fabs(n+1.) < 1e-6:
                    return math.sqrt((math.pow(c,2)*math.pow(math.log(xmax) - math.log(xmin),2)* 
                            ((math.pow(dx2,2)*math.pow(-2*math.log(x1) + math.log(xmax) + math.log(xmin),2))/ 
                            math.pow(x2,2) + (math.pow(dy2,2)*math.pow(-2*math.log(x1) + math.log(xmax) + math.log(xmin),2))/math.pow(y2,2) +  
                            (math.pow(dx1,2)*math.pow(-2*math.log(x2) + math.log(xmax) + math.log(xmin),2))/math.pow(x1,2) + (math.pow(dy1,2)* 
                            math.pow(-2*math.log(x2) + math.log(xmax) + math.log(xmin),2))/math.pow(y1,2)))/ 
                            (math.pow(xmax - xmin,2)*math.pow(math.log(x1) - math.log(x2),2)))/2.
                else: 
                    return math.sqrt(((math.pow(c,2)*math.pow(dx2,2)*math.pow(n,2)*math.pow(math.pow(xmax,1 + n)* 
                                (-1 - (1 + n)*math.log(x1) + (1 + n)*math.log(xmax)) + math.pow(xmin,1 + n)*(1 + (1 + n)*math.log(x1) - 
                                (1 + n)*math.log(xmin)),2))/math.pow(x2,2) + (math.pow(dy2,2)*math.pow(c*math.pow(xmax,1 + n)* 
                                (1 + (1 + n)*math.log(x1) - (1 + n)*math.log(xmax)) - c*math.pow(xmin,1 + n)*(1 + (1 + n)*math.log(x1) - 
                                (1 + n)*math.log(xmin)),2))/math.pow(y2,2) + (math.pow(dy1,2)*math.pow(c*math.pow(xmax,1 + n)* 
                                (1 + (1 + n)*math.log(x2) - (1 + n)*math.log(xmax)) - c*math.pow(xmin,1 + n)*(1 + (1 + n)*math.log(x2) - 
                                (1 + n)*math.log(xmin)),2))/math.pow(y1,2) + (math.pow(c,2)*math.pow(dx1,2)*math.pow(n,2)* 
                                math.pow(math.pow(xmax,1 + n)*(1 + (1 + n)*math.log(x2) - (1 + n)*math.log(xmax)) +  math.pow(xmin,1 + n)* 
                                (-1 - (1 + n)*math.log(x2) + (1 + n)*math.log(xmin)),2))/math.pow(x1,2))/ 
                                (math.pow(1 + n,4)*math.pow(xmax - xmin,2)*math.pow(math.log(x1) - math.log(x2),2)));
            else:        
                return math.sqrt((math.pow(c,2)*math.pow(x,2*n)*((math.pow(dx2,2)*math.pow(n,2)*math.pow(math.log(x) - 
                            math.log(x1),2))/math.pow(x2,2) + (math.pow(dy2,2)*math.pow(math.log(x) - math.log(x1),2))/math.pow(y2,2) +  
                            (math.pow(dx1,2)*math.pow(n,2)*math.pow(math.log(x) - math.log(x2),2))/math.pow(x1,2) + 
                            (math.pow(dy1,2)*math.pow(math.log(x) - math.log(x2),2))/math.pow(y1,2)))/math.pow(math.log(x1) - math.log(x2),2));
            # end of "pow"
        else: 
            return 0    
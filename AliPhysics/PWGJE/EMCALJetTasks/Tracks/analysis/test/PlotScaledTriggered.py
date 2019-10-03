'''
Created on 22.09.2014

@author: markusfasel
'''

from PWGJE.EMCALJetTasks.Tracks.analysis.base.Graphics import SinglePanelPlot, GraphicsObject, Style, Frame
from PWGJE.EMCALJetTasks.Tracks.analysis.correction.TriggeredSpectrumScaler import TriggeredSpectrumScaler
from PWGJE.EMCALJetTasks.Tracks.analysis.correction.SpectrumCombiner import SpectrumCombiner
from ROOT import kRed, kBlack, kBlue

class PlotScaledTriggeredToMinBias(SinglePanelPlot):
    '''
    classdocs
    '''


    def __init__(self, minbiasspectrum, triggeredSpectrum):
        '''
        Constructor
        '''
        SinglePanelPlot.__init__(self)
        self.__minbiasSpectrum = GraphicsObject(minbiasspectrum, Style(kRed,25))
        triggeredSpectrumMaker = TriggeredSpectrumScaler(minbiasspectrum, triggeredSpectrum)
        self.__triggeredSpectrum = GraphicsObject(triggeredSpectrumMaker.GetScaledTriggeredSpectrum(), Style(kBlue, 24))
        combinedSpectrumMaker = SpectrumCombiner(minbiasspectrum, self.__triggeredSpectrum.GetData())
        self.__combinedSpectrum = GraphicsObject(combinedSpectrumMaker.MakeCombinedSpectrum(50.), Style(kBlack, 20))
        self.__labeltext = None
        
    def SetLabel(self, label):
        self.__labeltext = label
        
    def Create(self):
        self._OpenCanvas("triggerSpectrumScalerPlot", "Compare scaled trigger to minbias")
        pad = self._GetFramedPad()
        #pad.GetPad().SetLogx()
        pad.GetPad().SetLogy()
        frame = Frame("framecomp", 0.1, 100, 1e-10, 2)
        frame.SetXtitle("p_{t} (GeV/c)")
        frame.SetYtitle("1/N_{ev} dN/dp_{t} ((GeV/c)^{-1})")
        pad.DrawFrame(frame)
        pad.DrawGraphicsObject(self.__combinedSpectrum, True, "Combined")
        pad.DrawGraphicsObject(self.__minbiasSpectrum, True, "MinBias")
        pad.DrawGraphicsObject(self.__triggeredSpectrum, True, "Triggered")
        pad.CreateLegend(0.55, 0.75, 0.89, 0.89)
        if self.__labeltext:
            pad.CreateLabel(0.15, 0.15, 0.45, 0.2, self.__labeltext)
        
#! /usr/bin/env python

from Helper import ReadHistList, MakeRatio, HistNotFoundException
from Graphics import Frame, Style
from SpectrumContainer import DataContainer
from ROOT import TCanvas, TLegend, TPaveText
from ROOT import kBlack, kRed, kBlue, kGreen, kOrange

gPlot = None

class Plot:
    """
    Class representation of the resulting plot
    The plot has n rows and 2 columns, one row per cut and within a column for the spectrum and a column for the ratio
    The first (upper) row will get a legend.
    """
    
    def __init__(self):
        self.__data = {}
        self.__canvas = None
        self.__legend = None
        self.__labels = []
        self.__axes = { "Spectrum" : Frame("tmplSpectrum", 0., 100., 1e-10, 100.), \
                        "Ratio" : Frame("tmplRatio", 0., 100., 0., 1000.)}
        self.__axes["Spectrum"].SetXtitle("p_{t} (GeV/c)")
        self.__axes["Ratio"].SetXtitle("p_{t} (GeV/c)")
        self.__axes["Spectrum"].SetYtitle("1/N_{event} 1/(#Delta p_{t}) dN/dp_{t} ((GeV/c)^{-2})")
        self.__axes["Ratio"].SetYtitle("Ratio to Min. Bias")
        self.__cutnames = {0 : "No cuts" , 1 : "Standard track cuts", 2 : "Hybrid track cuts"}
            
    def SetData(self, data, cutid):
        """
        Add data for a given cut id to the plot
        """
        self.__data[cutid] = data
        
    def GetCanvas(self):
        """ 
        Return resulting canvas
        """
        return self.__canvas
    
    def SaveAs(self, filenamebase):
        """
        Save plot as image file
        """
        types = ["eps", "pdf", "jpeg", "gif", "png"]
        for t in types:
            self.__canvas.SaveAs("%s.%s" %(filenamebase, t))
        
    def Create(self):
        """
        Create the final canvas
        """
        self.__canvas = TCanvas("plot", "Raw spectra comparison", 600, len(self.__data) * 300)
        self.__canvas.Divide(2, len(self.__data))
        row = 0
        for cut in self.__data.keys():
            self.__DrawCut(cut, row)
            row = row + 1
        return self.__canvas
    
    def __DrawCut(self, cutID, row):
        """
        Draw row with spectra comparison and ratios for a given cut combination
        """
        spectrumpad = self.__canvas.cd(row * 2 + 1)
        spectrumpad.SetGrid(False, False)
        spectrumpad.SetLogx()
        spectrumpad.SetLogy()
        self.__axes["Spectrum"].Draw()
        drawlegend = False
        if not self.__legend:
            self.__legend = TLegend(0.65, 0.55, 0.89, 0.89)
            self.__legend.SetBorderSize(0)
            self.__legend.SetFillStyle(0)
            self.__legend.SetTextFont(42)
            drawlegend = True
        for trg in self.__data[cutID]["Spectra"].keys():
            self.__data[cutID]["Spectra"][trg].Draw("epsame")
            if drawlegend:
                self.__legend.AddEntry(self.__data[cutID]["Spectra"][trg], trg, "lep")
        if drawlegend:
            self.__legend.Draw()
        
        cutlab = TPaveText(0.15, 0.15, 0.55, 0.22, "NDC")
        cutlab.SetBorderSize(0)
        cutlab.SetFillStyle(0)
        cutlab.SetTextFont(42)
        cutlab.AddText(self.__cutnames[cutID])
        cutlab.Draw()
        self.__labels.append(cutlab)
        
        ratiopad = self.__canvas.cd(row * 2 + 2)
        ratiopad.SetGrid(False, False)
        ratiopad.SetLogx()
        self.__axes["Ratio"].Draw()
        for trg in self.__data[cutID]["Ratios"].keys():
            self.__data[cutID]["Ratios"][trg].Draw("epsame")
        self.__canvas.cd()

def ReadSpectra(filename, triggers):
    """
    Read the spectra for different trigger classes from the root file
    Returns a dictionary of triggers - spectrum container
    """
    hlist = ReadHistList(filename, "PtEMCalTriggerTask")
    result = {}
    for trg in triggers:
        result[trg] = DataContainer(eventHist = hlist.FindObject("hEventHist%s" %(trg)), trackHist = hlist.FindObject("hTrackHist%s" %(trg)))
    return result

def MakeSpectraCut(inputdata, cutid):
    """
    Create for all trigger classes rawspectra for a given cut id and for events within +-10 cm in z-Vertex with 
    pileup rejection on and the ratios to min. bias events
    """
    styles = {"MinBias" : Style(kBlack, 20), "EMCJHigh" : Style(kRed, 24), "EMCJLow" : Style(kOrange, 26), "EMCGHigh" : Style(kBlue, 25), "EMCGLow" : Style(kGreen, 27)}
    rawspectra = {} 
    for trg in inputdata.keys():
        inputdata[trg].SetVertexRange(-10., 10.)
        inputdata[trg].SetPileupRejection(True)
        inputdata[trg].SelectTrackCuts(cutid)
        rawspectra[trg] = inputdata[trg].MakeProjection(0, "ptSpectrum%s", "p_{t} (GeV/c)", "1/N_{event} 1/(#Delta p_{t}) dN/dp_{t} ((GeV/c)^{-2})")
        rawspectra[trg].SetMarkerColor(styles[trg].GetColor())
        rawspectra[trg].SetLineColor(styles[trg].GetColor())
        rawspectra[trg].SetMarkerStyle(styles[trg].GetMarker())
        inputdata[trg].Reset()
    ratios = {}
    for trg in rawspectra.keys():
        if trg == "MinBias":
            continue
        ratios[trg] = MakeRatio(rawspectra[trg], rawspectra["MinBias"])
        ratios[trg].SetMarkerColor(styles[trg].GetColor())
        ratios[trg].SetLineColor(styles[trg].GetColor())
        ratios[trg].SetMarkerStyle(styles[trg].GetMarker())
    result = {"Spectra" : rawspectra, "Ratios" : ratios}
    return result


def MakeRawSpectraComparisonPlot(filename, doSave = False):
    """
    Create the final comparison plot
    """
    triggers = ["MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow"]
    data = ReadSpectra(filename, triggers)
    plot = Plot()
    for cut in range(0, 3):
        plot.SetData(MakeSpectraCut(data, cut), cut)
    plot.Create()
    gPlot = plot
    if doSave:
        plot.SaveAs("TriggercomparisonCuts")
    return plot
        
def main():
    """
    Main function: Delegate drawing 
    """
    inputfile = sys.argv[1]
    MakeSpectraComparisonPlot(inputfile, True)
    
if __name__ == "__main__":
    main()
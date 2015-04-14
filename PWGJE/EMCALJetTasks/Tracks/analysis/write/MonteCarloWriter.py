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
import os
from ROOT import TFile,TH1D,TList,TObject

from base.MonteCarloFileHandler import MonteCarloFileHandler

class BinContent(object):
    
    def __init__(self):
        self.__MCtruth = None
        self.__triggers = {}
        
    def SetMCtruth(self, mctruth):
        self.__MCtruth = mctruth
        
    def AddTrigger(self, triggername, content):
        self.__triggers[triggername] = content
        
    def MakeROOTPrimitive(self, name):
        result = TList()
        result.SetName(name)
        if self.__MCtruth:
            result.Add(self.__MCtruth)
        for hist in self.__triggers.itervalues():
            result.Add(hist)
        return result
    
class MonteCarloWriter(object):
    
    def __init__(self, isNew):
        '''
        Constructor
        '''
        self._weights = TH1D("weights", "Pythia weights", 11, -0.5, 10.5)
        self._nevents = {}
        self.__CreateEventContainers()
        self._isNew = isNew 
        self._pthardbins = {}
        self._listofbins = []
        self._inputcol = self.ReadData()
        
    def __CreateEventContainers(self):
        for trigger in ["MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow"]:
            self._nevents[trigger] = TH1D("nevents%s" %(trigger), "nevent for %s events" %(trigger), 11, -0.5, 10.5);
        
    def ReadData(self):
        reader = MonteCarloFileHandler(True)
        entries = os.listdir(os.getcwd())
        for mypthardbin in entries:
            if not str(mypthardbin).isdigit():
                continue
            print "Processing %s" %(mypthardbin)
            pthardbin = int(mypthardbin)
            reader.AddFile("%02d/AnalysisResults.root" %(pthardbin), pthardbin, isNew = self._isNew)
            self._listofbins.append(pthardbin)
        return reader.GetCollection()
    
    def SetNumberOfEvents(self, trigger, pthardbin, nevents):
        if not trigger in self._nevents.keys():
            return
        self._nevents[trigger].SetBinContent(self._nevents[trigger].GetXaxis().FindBin(pthardbin), nevents)

    def Convert(self):
        for mybin in self._listofbins:
            self._weights.SetBinContent(mybin+1, self._inputcol.GetWeigthHandler().GetWeight(mybin))
            self._pthardbins[mybin] = self.ProcessBin(int(mybin))
            
    def WriteResults(self):
        outputfile = TFile(self.CreateOutputFilename(), "RECREATE")
        outputfile.cd()
        self._weights.Write(self._weights.GetName(), TObject.kSingleKey)
        for trigger in self._nevents.values():
            trigger.Write(trigger.GetName(), TObject.kSingleKey)
        for mybin in self._pthardbins:
            bindata = self._pthardbins[mybin].MakeROOTPrimitive("bin%d" %(mybin))
            bindata.Write(bindata.GetName(), TObject.kSingleKey)
        outputfile.Close()
        
    # pure virtual methods:
    def CreateOutputFilename(self):
        pass
    
    def ProcessBin(self, mybin):
        pass

class TrackWriter(MonteCarloWriter):
    '''
    Class Writing projected raw spectrum and MC truth 
    '''

    def __init__(self, isNew):
        '''
        Constructor
        '''
        MonteCarloWriter.__init__(self, isNew)
        self.__inAcceptance = False
        self.__MCKine = False
        self.__etacut = None
        self.__phicut = None
        
    def SetInAcceptance(self):
        self.__inAcceptance = True
        
    def SetAll(self):
        self.__inAcceptance = False
        
    def SetMCKine(self):
        self.__MCKine = True
        
    def SetRecKine(self):
        self.__MCKine = False
        
    def SetEtaCut(self, etaMin, etaMax, tag):
        self.__etacut = {"etaMin":etaMin, "etaMax":etaMax, "tag":tag}
        
    def SetPhiCut(self, phiMin, phiMax, tag):
        self.__phicut = {"phiMin":phiMin, "phiMax":phiMax, "tag":tag}
        
    def ProcessBin(self, mybin):
        results = BinContent()
        bindata = self._inputcol.GetData(mybin)
        results.SetMCtruth(self.ProjectMCtruth(bindata.GetMCTruth(), "MCTruthbin%d" %(mybin)))
        kinestring = "tracksMCKine" if self.__MCKine else "tracks"
        acceptancestring = "WithClusters" if self.__inAcceptance else "All"
        for trigger in ["MinBias","EMCJHigh","EMCJLow","EMCGHigh","EMCGLow"]:
            histname = "%s%s" %(kinestring, acceptancestring)
            print "histname: %s" %(histname)
            tc = bindata.GetData(trigger).FindTrackContainer(histname)
            if self.__etacut:
                print "Using eta range %f %f" %(self.__etacut["etaMin"], self.__etacut["etaMax"])
                tc.SetEtaRange(self.__etacut["etaMin"], self.__etacut["etaMax"])
            if self.__phicut:
                tc.SetPhiRange(self.__phicut["phiMin"], self.__phicut["phiMax"])
            sn = "%sbin%d" %(trigger, mybin)
            spectrum = self.Project(tc, sn)
            results.AddTrigger(trigger, spectrum)
            self.SetNumberOfEvents(trigger, mybin, tc.GetEventCount())
        return results
    
    def CreateOutputFilename(self):
        kinestring="etaall"
        if self.__etacut:
            kinestring = "eta%s" %(self.__etacut["tag"])
        if self.__phicut:
            kinestring +=  self.__phicut["tag"]
        return "MonteCarloProjected%s%sKine%s.root" %("Acc" if self.__inAcceptance else "All", "MC" if self.__MCKine else "Rec", kinestring)
    
    def ProjectMCtruth(self, inputcontainer, outputname):
        print "Projecting MC-truth"
        inputcontainer.ApplyCut(inputcontainer.GetAxisDefinition().GetAxisName(3),-10., 10.)
        #inputcontainer.ApplyCut(inputcontainer.GetAxisDefinition().GetAxisName(4), 1, 1)
        print "finish applying cut"
        if self.__etacut:
            print "apply eta cut"
            inputcontainer.ApplyCut("eta", self.__etacut["etaMin"], self.__etacut["etaMax"])
        return inputcontainer.Projection1D(outputname, inputcontainer.GetAxisDefinition().GetAxisName(0))

    def Project(self, inputcontainer, outputname):
        inputcontainer.SetVertexRange(-10., 10.)
        inputcontainer.SetPileupRejection(True)
        inputcontainer.SelectTrackCuts(1)
        return inputcontainer.MakeProjection(0, outputname, doNorm = False)
        
        
class ClusterWriter(MonteCarloWriter):
    
    def __init__(self, isNew):
        MonteCarloWriter.__init__(self, isNew)
        self.__calibrated = True
        
    def SetCalibrated(self):
        self.__calibrated = True
    
    def SetUncalibrated(self):
        self.__calibrated = False
        
    
    def ProcessBin(self, pthatbin):
        """
        Make projections of the different cluster hists
        """
        results = BinContent()
        bindata = self._inputcol.GetData(pthatbin)
        for trigger in ["MinBias","EMCJHigh","EMCJLow","EMCGHigh","EMCGLow"]:
            clustercont = bindata.GetData(trigger).FindClusterContainer("Calib" if self.__calibrated else "Uncalib")
            spectrum = self.ProjectContainer(clustercont, "%sbin%d" %(trigger, pthatbin))
            results.AddTrigger(trigger, spectrum)
            self.SetNumberOfEvents(trigger, pthatbin, clustercont.GetEventCount())
        return results

    def ProjectContainer(self, inputcontainer, outputname):
        inputcontainer.SetVertexRange(-10, 10)
        inputcontainer.SetPileupRejection(True)
        return inputcontainer.MakeProjection(0, outputname, doNorm = False)
    
    def CreateOutputFilename(self):
        return "MonteCarloClusterProjection%s.root" %("Calib" if self.__calibrated else "Uncalib")
    
class JetData(object):
    """
    Container object for jet based histogram with a given minimum jet pt
    """
    
    def __init__(self, jetpt, trigger):
        """
        Constructor
        """
        self.__jetpt = jetpt
        self.__trigger = trigger
        self.__spectra = []
        
    def ROOTify(self):
        """
        Create simple primitive ROOT object structure
        """
        outputlist = TList()
        outputlist.SetName("JetSpectraPt%03d" %(int(self.__jetpt)))
        for spec in self.__spectra:
            outputlist.Add(spec)
        return outputlist
        
    def AddSpectrum(self, spectrum, isMCkine):
        """
        Add new object to data set
        """
        histname = "spectrumTrackJetData%s%s" %("MCKine" if isMCkine else "RecKine", self.__trigger)
        spectrum.SetName(histname)
        self.__spectra.append(spectrum)
        
        
class JetWriter(MonteCarloWriter):
    
    def __init__(self, isNew):
        MonteCarloWriter.__init__(self, isNew)
        
    def ProcessBin(self, mybin):
        results = BinContent()
        bindata = self._inputcol.GetData(mybin)
        for trigger in ["MinBias","EMCJHigh","EMCJLow","EMCGHigh","EMCGLow"]:
            print "Doing trigger %s" %(trigger)
            jetcont =  bindata.GetData(trigger).GetJetContainer()
            outputcont = TList()
            outputcont.SetName(trigger)
            for jetpt in jetcont.GetListOfJetPts():
                print "Inspecting jet pt %f" %(jetpt)
                jetdat = JetData(jetpt,trigger)
                projectedRec = jetcont.MakeProjectionRecKine(jetpt, 0, "projectedPtRec")
                projectedMC = jetcont.MakeProjectionMCKine(jetpt, 0, "projectedPtMC")
                jetdat.AddSpectrum(projectedRec, False)
                jetdat.AddSpectrum(projectedMC, True)
                outputcont.Add(jetdat.ROOTify())
            results.AddTrigger(trigger, outputcont)
        return results
    
    def CreateOutputFilename(self):
        return "MCTracksInJets.root"
    
def RunTrackProjection(doAcc = False, doMCKine = False, etaSel = "all", phiSel = False, isNew = False):
    writer = TrackWriter(isNew)
    if doAcc:
        writer.SetInAcceptance()
    else:
        writer.SetAll()
    if doMCKine:
        writer.SetMCKine()
    else:
        writer.SetRecKine()
    if etaSel == "centcms":
        writer.SetEtaCut(-0.79999, -0.3001, "centcms")
    elif etaSel == "emc":
        writer.SetEtaCut(-0.5, 0.5, "emc")
    if phiSel:
        writer.SetPhiCut(1.5, 3.1, "emcphi") 
    writer.Convert()
    writer.WriteResults()
    
def RunClusterProjection(doCalib = True, isNew = False):
    writer = ClusterWriter(isNew)
    if doCalib:
        writer.SetCalibrated()
    else:
        writer.SetUncalibrated()
    writer.Convert()
    writer.WriteResults()
    
def RunJetProjection(isNew = False):
    writer = JetWriter(isNew)
    writer.Convert()
    writer.WriteResults()
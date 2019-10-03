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
from PWGJE.EMCALJetTasks.Tracks.analysis.base.FileHandler import LegoTrainFileReader
from ROOT import TFile, TList, TObject, TH1F

class DataSpectra(object):
    
    class TriggerData(object):

        def __init__(self, triggername, events, spec):
            self.__triggername = triggername
            self.__eventcount = events
            self.__projectedSpectrum = spec

        def get_triggername(self):
            return self.__triggername

        def get_eventcount(self):
            return self.__eventcount

        def get_projected_spectrum(self):
            return self.__projectedSpectrum
            
        def set_triggername(self, value):
            self.__triggername = value

        def set_eventcount(self, value):
            self.__eventcount = value

        def set_projected_spectrum(self, value):
            self.__projectedSpectrum = value
    
        def del_eventcount(self):
            del self.__eventcount

        def del_projected_spectrum(self):
            del self.__projectedSpectrum
        
        def del_triggername(self):
            del self.__triggername
    
        def MakeROOTPrimitive(self):
            result = TList()
            result.SetName(self.triggername)
            self.__eventcount.SetName("events")
            self.__projectedSpectrum.SetName("spectrum")
            result.Add(self.__eventcount)
            result.Add(self.__projectedSpectrum)
            return result
    
        eventcount = property(get_eventcount, set_eventcount, del_eventcount, "Histogram with the event selection")
        projectedSpectrum = property(get_projected_spectrum, set_projected_spectrum, del_projected_spectrum, "Histogram with the not-normalised projected spectrum")
        triggername = property(get_triggername, set_triggername, del_triggername, "Name of the trigger class")
    
    def __init__(self):
        self.__triggers = {}

    def get_triggers(self):
        return self.__triggers

    def set_triggers(self, value):
        self.__triggers = value

    def del_triggers(self):
        del self.__triggers

    triggers = property(get_triggers, set_triggers, del_triggers, "Trigger data list")
        
    def AddTrigger(self, triggerName, spectrum, events):
        self.__triggers[triggerName] = self.TriggerData(triggerName, events, spectrum)
        
    def GetListOfROOTPrimitives(self):
        rootprim = []
        for trigger in self.__triggers.itervalues():
            rootprim.append(trigger.MakeROOTPrimitive())
        return rootprim

class DataWriter(object):
    
    def __init__(self, filename, isNew):
        self._inputdata = self.__ReadFile(filename, isNew)
        self._outputdata = DataSpectra()
        
    def __ReadFile(self, filename, isNewStruct):
        reader = LegoTrainFileReader(filename, isNew = isNewStruct, trackCuts="")
        return reader.ReadFile()
    
    def Convert(self):
        for trigger in ["MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh", "EMCGLow"]:
            trdata = self._ProcessTrigger(trigger, self._inputdata.GetData(trigger))
            self._outputdata.AddTrigger(trigger, trdata["spectrum"], trdata["events"])
        
    def WriteOutput(self):
        myoutput = TFile(self._GetOutputFile(), "RECREATE")
        myoutput.cd()
        for obj in self._outputdata.GetListOfROOTPrimitives():
            obj.Write(obj.GetName(), TObject.kSingleKey)
        myoutput.Close()
        
    def _MakeProjection(self, trackContainer):
        return trackContainer.MakeProjection(0, "spectrum", doNorm = False)
    
    def _GetNumberOfEvents(self, trackContainer):
        eventHist = TH1F("events", "Events", 1, 0.5, 1.5)
        eventHist.SetBinContent(1, trackContainer.GetEventCount())
        return eventHist
         
    # Pure virtual functions
    def _ProcessTrigger(self, trigger, dset):
        print "Method needs to be implemented by inheriting classes"
        
    def _GetOutputFile(self):
        print "Method needs to be implemented by inheriting classes"
        
class DataTrackWriter(DataWriter):
    
    def __init__(self, filename, isNew):
        DataWriter.__init__(self, filename, isNew)
        self.__etacut = None
        self.__phicut = None
        self.__inAcceptance = False
        
    def SetEtaCut(self, tag, emin, emax):
        self.__etacut = {"etamin":emin, "etamax":emax, "tag":tag}
        
    def SetPhiCut(self, tag, phimin, phimax):
        self.__phicut = {"phimin":phimin, "phimax":phimax, "tag": tag}
        
    def SetInAcceptance(self, inAcceptance = True):
        self.__inAcceptance = inAcceptance
        
    def _GetOutputFile(self):
        kinestring="etaall"
        if self.__etacut:
            kinestring = "eta%s" %(self.__etacut["tag"])
        if self.__phicut:
            kinestring +=  self.__phicut["tag"]
        accString = "All"
        if self.__inAcceptance:
            accString = "InAcceptance"
        return "DataTracks%s%s.root" %(accString, kinestring)
    
    def __DefineTracks(self, tc):
        tc.SetVertexRange(-10, 10)
        tc.SetPileupRejection(True)
        tc.SelectTrackCuts(1)
        if self.__etacut:
            tc.SetEtaRange(self.__etacut["etamin"], self.__etacut["etamax"])
        if self.__phicut:
            tc.SetPhiRange(self.__phicut["phimin"], self.__phicut["phimax"])
    
    def _ProcessTrigger(self, trigger, dset):
        containerName = "tracksAll"
        if self.__inAcceptance:
            containerName = "tracksWithClusters" 
        tc = dset.FindTrackContainer(containerName)
        self.__DefineTracks(tc)
        projected = self._MakeProjection(tc)
        nevents = self._GetNumberOfEvents(tc)
        return {"spectrum":projected, "events":nevents}

        
class DataClusterWriter(DataWriter): 
    
    def __init__(self, filename, isNew):
        DataWriter.__init__(self, filename, isNew)
        self.__useCalibrated = True
        
    def SetUseCalibratedClusters(self, doUse = True):
        self.__useCalibrated = doUse
        
    def _GetOutputFile(self):
        calibstring = "Uncalib"
        if self.__useCalibrated:
            calibstring = "Calib"
        return "DataCluster%s.root" %(calibstring)
    
    def _ProcessTrigger(self, trigger, dset):
        dset.Print()
        containerName = "Uncalib"
        if self.__useCalibrated:
            containerName = "Calib"
        cc = dset.FindClusterContainer(containerName)
        self.__DefineClusters(cc)
        projected = self._MakeProjection(cc)
        nevents = self._GetNumberOfEvents(cc)
        return {"spectrum":projected, "events":nevents}
    
    def __DefineClusters(self, cc):
        cc.SetVertexRange(-10, 10)
        cc.SetPileupRejection(True)
    
def WriteTracks(filename, inAcceptance = False, etaSel = "all", isNew = True):
    writer = DataTrackWriter(filename, isNew)
    writer.SetInAcceptance(inAcceptance)
    if etaSel == "centcms":
        writer.SetEtaCut("centcms", -0.7999, -0.3001)
    writer.Convert()
    writer.WriteOutput()

def WriteClusters(filename, calib = True, isNew = True):
    writer = DataClusterWriter(filename, isNew)
    writer.SetUseCalibratedClusters(calib)
    writer.Convert()
    writer.WriteOutput()
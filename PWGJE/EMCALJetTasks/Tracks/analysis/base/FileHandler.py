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
from ROOT import TFile,TIter,TObject,gDirectory,gROOT
from PWGJE.EMCALJetTasks.Tracks.analysis.base.DataSet import DataSet
from PWGJE.EMCALJetTasks.Tracks.analysis.base.FileResults import ResultData
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.ParticleTHnSparse import ParticleTHnSparse
from PWGJE.EMCALJetTasks.Tracks.analysis.base.struct.DataContainerFactory import DataContainerFactory
    
class FileReader(object):
    
    class FileReaderException(Exception):
        """
        Exception class handling root files which are
        either not found or not readable.
        """

        def __init__(self, filename):
            """
            Constructor, assigning the name of the file which 
            failed to be read.
            """
            self.filename = filename

        def __str__(self):
            """
            Create string representation of the error message
            """
            return "Could not open file %s" %(self.filename)

    
    def __init__(self, filename, isMC = False):
        """
        Initialise file reader with name of the file to read
        """
        self.__filename = filename
        self.__directory = None
        self.__isMC = isMC
        self.__isReadWeights = False
        self.__weightlist = None
        self.__datafactory = DataContainerFactory("new")
        self._histlist = "histosPtEMCalTriggerHistograms"
        self._trackCutsTag = ""
        
    def SetHistList(self, histlist):
        self._histlist = histlist
        
    def SetTrackCuts(self, tc):
        self._trackCutsTag = tc
        
    def SetReadWeights(self):
        """
        Read also histograms for the weight calculation
        """
        self.__isReadWeights = True
        
    def GetDataFormat(self):
        return "old" if self._histlist == "histosPtEMCalTriggerHistograms" else "new"
        
    def GetWeightHistograms(self):
        """
        Access to weight histograms
        """
        return self.__weightlist
        
    def SetDirectory(self, dirname):
        """
        Set directory inside the rootfile
        """
        self.__directory = dirname
    
    def ReadFile(self):
        """
        Read rootfile and create a ResultData structure with all the histograms sorted according
        to trigger classes. Raising FileReaderExceptions if the file can't be opened, doesn't contain
        the directory or list, or has an empty histogram list
        """
        filecontent = self.__ReadHistList()
        if self.__isReadWeights:
            self.__weightlist = filecontent["weights"]
        hlist = filecontent["spectra"]
        if not hlist.GetEntries():
            raise self.FileReaderException("Empty list of histograms in file %s" %(self.__filename))
        result = ResultData("result")
        
        # build list of histograms and extract trigger names
        histnames = []
        for oID in range(0, hlist.GetEntries()):
            histnames.append(hlist.At(oID).GetName())    
        triggers = []
        for hname in histnames:
            if not "hEventHist" in hname:
                continue
            triggers.append(hname.replace("hEventHist",""))
        print "Found the following triggers:"
        print "================================="
        isFirst = True
        trgstring = ""
        for trg in triggers:
            if isFirst:
                trgstring += trg 
                isFirst = False
            else:
                trgstring += ", %s" %(trg)
        print trgstring
        
        self.__datafactory.SetDataFormat(self.GetDataFormat())
        
        # Handle MC-truth data
        if self.__isMC:
            result.SetMCTruth(self.__datafactory.CreateParticleContainer(hlist.FindObject("hMCtrueParticles")))
        
        # Add the result hists to the result container
        for trigger in triggers:
            eventhist = hlist.FindObject("hEventHist%s" %(trigger))
            trackhist = hlist.FindObject("hTrackHist%s" %(trigger))
            #eventhist.Sumw2()
            #trackhist.Sumw2()
            triggerdata = DataSet()
            triggerdata.AddEventHistForJets(eventhist)
            triggerdata.AddTrackContainer("tracksAll", self.__datafactory.CreateTrackContainer(eventhist, trackhist))
            tracksWithClusters = hlist.FindObject("hTrackInAcceptanceHist%s" %(trigger)) 
            if tracksWithClusters:
                triggerdata.AddTrackContainer("tracksWithClusters", self.__datafactory.CreateTrackContainer(eventhist, tracksWithClusters))
            tracksMCKine = hlist.FindObject("hMCTrackHist%s" %(trigger))
            if tracksMCKine:
                triggerdata.AddTrackContainer("tracksMCKineAll", self.__datafactory.CreateTrackContainer(eventhist, tracksMCKine))
            tracksMCKineWithClusters = hlist.FindObject("hMCTrackInAcceptanceHist%s" %(trigger)) 
            if tracksMCKineWithClusters:
                triggerdata.AddTrackContainer("tracksMCKineWithClusters", self.__datafactory.CreateTrackContainer(eventhist, tracksMCKineWithClusters))
            clusterhists = ["hClusterCalibHist","hClusterUncalibHist"]
            for clust in clusterhists:
                clhist = hlist.FindObject("%s%s" %(clust, trigger))
                if clhist:
                    tag = clust.replace("hCluster","").replace("Hist","")
                    #clhist.Sumw2()
                    triggerdata.AddClusterContainer(tag, self.__datafactory.CreateClusterContainer(eventhist, clhist))
            self.ProcessJets(trigger, triggerdata, hlist)
            result.SetData(trigger, triggerdata)
        return result
    
    def ProcessJets(self, triggerclass, dataset, histlist):
        """
        Fill jet hists to the histogram container
        
        1. find all histograms for the given trigger class that contain the trigger class name and Jet
        2. Group them according to jet pt and histogram type
        """
        histiter = TIter(histlist)
        histfound = histiter.Next()
        histlist = []
        while histfound:
            histname = str(histfound.GetName())
            if triggerclass in histname and "TrackJetHist" in histname:
                histlist.append(histfound)
            histfound = histiter.Next()
            
        for jethist in histlist:
            histname = str(jethist.GetName())
            jetpt = self.__GetJetPt(histname)
            dataset.AddJetSpectrum(jethist,jetpt, True if "hMC" in histname else False)
            
    def __GetJetPt(self, histname):
        start = histname.index("jetPt") + 5
        ptstring = histname[start:start + 3]
        return int(ptstring)
        
    def __ReadHistList(self):
        """
        Read the list of histograms from a given rootfile
        optionally the list can be wihtin a directory (i.e. when analysing lego train output)
        """
        result = {"spectra":None, "weights":None}
        inputfile = TFile.Open(self.__filename)
        if not inputfile or inputfile.IsZombie():
            raise self.FileReaderException(self.__filename)
        mydirectory = None
        path = self.__filename
        if self.__directory:
            path += "#%s" %(self.__directory)
            if not inputfile.cd(self.__directory):
                inputfile.Close()
                raise self.FileReaderException(path)
            else:
                mydirectory = gDirectory
        else:
            mydirectory = inputfile
            path += "#"
            path += "#"
        rlist = mydirectory.Get("results")  # old structure
        if not rlist:
            rlist = mydirectory.Get("TriggerTracksResults%s" %(self._trackCutsTag))
        if self.__isReadWeights:
            result["weights"] = {"crosssection":rlist.FindObject("fHistXsection"), "trials":rlist.FindObject("fHistTrials")}
        hlist = rlist.FindObject(self._histlist)
        inputfile.Close()
        if not hlist:
            raise self.FileReaderException("%s/%s" %(path, self._histlist))
        result["spectra"] = hlist
        return result
    
class LegoTrainFileReader(FileReader):
    """
    File reader adapted to the file format in the lego train
    """
    
    def __init__(self, filename, trackCuts = "standard", isMC = False, isNew = True):
        """
        Initialise file reader with filename and set the directory according to the definition in
        the lego train
        """
        FileReader.__init__(self, filename, isMC)
        self.SetDirectory("PtEMCalTriggerTask%s" %(trackCuts))
        self.SetTrackCuts(trackCuts)
        if isNew:
            self.SetHistList("histosptemcaltriggertask%s" %(trackCuts))
        
class ResultStructureReader(object):
    
    def __init__(self, filename):
        self.__filename = filename
        
    def ReadFile(self):
        return ResultData.BuildFromRootFile(self.__filename, "read")

def TestFileReader(filename):
    """
    Test procedure for the lego train file reader
    """
    testreader = LegoTrainFileReader(filename)
    return testreader.ReadFile()

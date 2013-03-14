/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$

///
/// \ingroup macros
///
/// \file CreateWeightedRejectList.C
///
/// \brief Macro test for the creation of a RejectList object (for simulations) in the OCDB, for the MUON Tracke 
///
/// Usage:
///
/// root[0] .L CreateWeightedRejectList.C+
/// root[1] CreateWeightedRejectList("runlist.txt","local://$HOME/myLocalOCDB");
///
/// where runlist.txt has 2 integers per line = "run nevents"
/// where nevents is the number of events where there's (for instance) a CMUS1B trigger
///
///
/// Assuming the file coming from the Export feature of the logbook is runlist.logbook.txt =>
///
/// 115521;PHYSICS_1;1357;0.9
/// 115516;PHYSICS_1;944;0.9
///
/// The awk command below will output the needed format for this macro :
///
/// awk '{split ($0,a,";"); print a[1] " " a[3];}' runlist.logbook.txt =>
///
/// 115521 1357
/// 115516 944
///
/// \author Matthieu Lenhardt, Subatech
///

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliDAQ.h"
#include "AliDCSValue.h"
#include "AliMUONCDB.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONPadStatusMapMaker.h"
#include "AliMUONRecoParam.h"
#include "AliMUONRejectList.h"
#include "AliMUONRejectList.h"
#include "AliMUONTrackerData.h"
#include "AliMUONVCalibParam.h"
#include "AliMpBusPatch.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDCSNamer.h"
#include "AliMpDDLStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliSysInfo.h"
#include "TEnv.h"
#include <Riostream.h>
#include <TFile.h>
#include <TGrid.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TTree.h>
#include <vector>
#include "AliCounterCollection.h"

#include "AliMUONTrackerDataWrapper.h"
#include "AliMUONPainterDataRegistry.h"

class RunInfo
{
public:
  RunInfo(Int_t run, Long64_t nevents) : fNumber(run), fNevents(nevents) { }

  Int_t Number() const { return fNumber; }

  Long64_t Nevents() const { return fNevents; }

  void SetNevents(Long64_t n) { fNevents = n; }

private:
  Int_t fNumber;
  Long64_t fNevents;
};

std::ostream& operator<<(std::ostream& out, const RunInfo& run)
{
  out << Form("RUN %09d NEVENTS %lld",run.Number(),run.Nevents());
  return out;
}

AliMUONRejectList* CheckDE_BP_ManuPedestals(Int_t rejectMask, AliMUONPadStatusMaker& status);

Int_t AddEventsSingleRun(Int_t run_number, Long64_t nEvents,AliMUONRejectList& rejectedEvents);

std::vector<RunInfo> runs;

//______________________________________________________________________________
AliMUONRejectList* CreateWeightedRejectList(AliCounterCollection& cc,
                                            const char* ocdbPath="local:///Users/laurent/Alice/OCDBcopy2013",
                                            const char* weigthWithTrigger="CMSL7-B-NOPF-MUON")
{
  // create the reject using the runs in the collection
  //
  // the counter collection is expected to have "trigger","event","run" rubrics
  // and event:PSALL should be a valid one.
  //
  
  if (!gGrid)
  {
    TGrid::Connect("alien://");
  }
  
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  
  AliMUONRejectList* weightedRejectList = new AliMUONRejectList;
  
  AliMUONRejectList rejectedEvents;

  TString sruns = cc.GetKeyWords("run");
  TObjArray* truns = sruns.Tokenize(",");
  
  TIter next(truns);
  TObjString* s;
  
  Double_t nEventsTotal(0.0);
  
  while ( ( s = static_cast<TObjString*>(next())) )
  {
    Int_t runNumber = s->String().Atoi();
    AliSysInfo::Instance()->AddStamp(Form("RUN%d",runNumber));
    Long64_t nEvents = cc.GetSum(Form("trigger:%s/event:PSALL/run:%d",weigthWithTrigger,runNumber));
    if (!nEvents)
    {
      std::cout << "No events for run " << runNumber << " ??" << std::endl;
      continue;
    }
    nEventsTotal += AddEventsSingleRun(runNumber,nEvents,rejectedEvents);
  }
  
  delete truns;
  
  AliMpDEIterator DEiter;
  
  for (DEiter.First(); !DEiter.IsDone(); DEiter.Next())
    if (DEiter.CurrentDEId() < 1100)
    {
      Int_t DEid = DEiter.CurrentDEId();
      
      Double_t nBadEventsDE = rejectedEvents.DetectionElementProbability(DEid);
      Double_t nEventsTotalDE = static_cast<Double_t>(nEventsTotal);
      Float_t probaDE = 0.0;
      if (nEventsTotalDE != 0.0)
        probaDE = nBadEventsDE / nEventsTotalDE;
      weightedRejectList->SetDetectionElementProbability(DEid, probaDE);
      
      Int_t nBusPatches = DEiter.CurrentDE()->GetNofBusPatches();
      for (Int_t ii = 0; ii < nBusPatches; ii++)
      {
        Int_t BPid = DEiter.CurrentDE()->GetBusPatchId(ii);
        
        Double_t nBadEventsBP = rejectedEvents.BusPatchProbability(BPid);
        Double_t nEventsTotalBP = nEventsTotalDE - nBadEventsDE;
        Float_t probaBP = 0.0;
        if (nEventsTotalBP != 0.0)
          probaBP = nBadEventsBP / nEventsTotalBP;
        weightedRejectList->SetBusPatchProbability(BPid, probaBP);
        
        Int_t nManus = AliMpDDLStore::Instance(kFALSE)->GetBusPatch(BPid, kFALSE)->GetNofManus();
        for (Int_t jj = 0; jj < nManus; jj++)
	      {
          Int_t manuId = AliMpDDLStore::Instance(kFALSE)->GetBusPatch(BPid, kFALSE)->GetManuId(jj);
          
          Double_t nBadEventsManu = rejectedEvents.ManuProbability(DEid, manuId);
          Double_t nEventsTotalManu = nEventsTotalBP - nBadEventsBP;
          Float_t probaManu = 0.0;
          if (nEventsTotalManu != 0.0)
            probaManu = nBadEventsManu / nEventsTotalManu;
          weightedRejectList->SetManuProbability(DEid, manuId, probaManu);
          
          for (Int_t channel = 0; channel < AliMpConstants::ManuNofChannels(); channel++)
          {
            Double_t nBadEventsChannel = rejectedEvents.ChannelProbability(DEid, manuId, channel);
            Double_t nEventsTotalChannel = nEventsTotalManu - nBadEventsManu;
            Float_t probaChannel = 0.0;
            if (nEventsTotalChannel != 0.0)
            {
              probaChannel = nBadEventsChannel / nEventsTotalChannel;
            }
            
            weightedRejectList->SetChannelProbability(DEid, manuId, channel, probaChannel);
          }
	      }
      }
    }
  
  
  AliMUONTrackerData* td = new AliMUONTrackerData("WeightedRejectList","RejectList",*weightedRejectList);
  
  AliMUONVTrackerDataMaker* dw = new AliMUONTrackerDataWrapper(td);
  
  AliMUONPainterDataRegistry::Instance()->Register(dw);

  return weightedRejectList;
}

//______________________________________________________________________________
bool CreateWeightedRejectList(const char* runlistfile="runlist.txt", 
                              const char* rejectListPath="local://$HOME/OCDB",
                              const char* author="Matthieu Lenhardt")
{
  /// Create a weighted RejectList for the runs included in the run list
  /// The cuts are the same that the one applied in the RecoParam used to create the ESDs.
  ///
  /// \param runlistfile : an ASCII file where each line is a pair "run nevents"
  /// \param rejectListPath : set this to the folder where the RejectList will be created (must have right access to it, of course)
  
  TGrid::Connect("alien://");
  
  ifstream in(runlistfile);
  int run, nevents;
  
  while ( in >> run >> nevents ) 
  {
    runs.push_back(RunInfo(run,nevents));
  };
  
  for ( std::vector<RunInfo>::size_type i = 0; i < runs.size(); ++i )
  {
    cout << runs[i] << endl;
  }
  
  if (runs.empty()) return false;
  
  Int_t firstRun = runs[0].Number();
  Int_t lastRun = runs[runs.size()-1].Number();
  
  const char* rawOCDB = "raw://";
  
  AliCDBManager* manager = AliCDBManager::Instance();

  manager->SetCacheFlag(kFALSE);
  
  manager->SetDefaultStorage(rawOCDB);

  manager->SetSpecificStorage("MUON/Calib/RejectList", rejectListPath);
  
  AliMUONRejectList weightedRejectList;
  AliMUONRejectList rejectedEvents;
  
  Long64_t nEventsTotal = 0;
  
  for (std::vector<RunInfo>::size_type ii = 0; ii < runs.size(); ++ii)
  {
    Int_t runNumber = runs[ii].Number();
    AliSysInfo::Instance()->AddStamp(Form("RUN%d",runNumber));
    Long64_t nEvents = runs[ii].Nevents();
    nEventsTotal += AddEventsSingleRun(runNumber, nEvents, rejectedEvents);
  }  
  
  AliMpDEIterator DEiter;
  
  for (DEiter.First(); !DEiter.IsDone(); DEiter.Next())
    if (DEiter.CurrentDEId() < 1100)
    {
      Int_t DEid = DEiter.CurrentDEId();
      
      Double_t nBadEventsDE = rejectedEvents.DetectionElementProbability(DEid);
      Double_t nEventsTotalDE = static_cast<Double_t>(nEventsTotal);
      Float_t probaDE = 0.0;
      if (nEventsTotalDE != 0.0)
        probaDE = nBadEventsDE / nEventsTotalDE;
      weightedRejectList.SetDetectionElementProbability(DEid, probaDE);
      
      Int_t nBusPatches = DEiter.CurrentDE()->GetNofBusPatches();
      for (Int_t ii = 0; ii < nBusPatches; ii++)
      {
        Int_t BPid = DEiter.CurrentDE()->GetBusPatchId(ii);
        
        Double_t nBadEventsBP = rejectedEvents.BusPatchProbability(BPid);
        Double_t nEventsTotalBP = nEventsTotalDE - nBadEventsDE;
        Float_t probaBP = 0.0;
        if (nEventsTotalBP != 0.0)
          probaBP = nBadEventsBP / nEventsTotalBP;
        weightedRejectList.SetBusPatchProbability(BPid, probaBP);
        
        Int_t nManus = AliMpDDLStore::Instance(kFALSE)->GetBusPatch(BPid, kFALSE)->GetNofManus();
        for (Int_t jj = 0; jj < nManus; jj++)
	      {
          Int_t manuId = AliMpDDLStore::Instance(kFALSE)->GetBusPatch(BPid, kFALSE)->GetManuId(jj);
          
          Double_t nBadEventsManu = rejectedEvents.ManuProbability(DEid, manuId);
          Double_t nEventsTotalManu = nEventsTotalBP - nBadEventsBP;
          Float_t probaManu = 0.0;
          if (nEventsTotalManu != 0.0)
            probaManu = nBadEventsManu / nEventsTotalManu;
          weightedRejectList.SetManuProbability(DEid, manuId, probaManu);
          
          for (Int_t channel = 0; channel < AliMpConstants::ManuNofChannels(); channel++)
          {
            Double_t nBadEventsChannel = rejectedEvents.ChannelProbability(DEid, manuId, channel);
            Double_t nEventsTotalChannel = nEventsTotalManu - nBadEventsManu;
            Float_t probaChannel = 0.0;
            if (nEventsTotalChannel != 0.0)
            {
              probaChannel = nBadEventsChannel / nEventsTotalChannel;
            }
            
            weightedRejectList.SetChannelProbability(DEid, manuId, channel, probaChannel);
          }
	      }
      }
    }
  
  
  AliMUONCDB::WriteToCDB(&weightedRejectList, "MUON/Calib/RejectList",
                         firstRun , lastRun, 
                         "Weighted reject List for MCH, for simulations only", 
                         author);
  
  return true;
}

//____________________________________________________________________________________________
AliMUONRejectList* CheckDE_BP_ManuPedestals(Int_t rejectMask, AliMUONPadStatusMaker& status)
{
  AliMUONRejectList* nBadChannels = new AliMUONRejectList;
  
  AliMpDEIterator DEiter;
  for (DEiter.First(); !DEiter.IsDone(); DEiter.Next())
    if (DEiter.CurrentDEId() < 1100)
    {
      Int_t DEid = DEiter.CurrentDEId();
      Int_t nBusPatches = DEiter.CurrentDE()->GetNofBusPatches();
      for (Int_t ii = 0; ii < nBusPatches; ii++)
      {
        Int_t BPid = DEiter.CurrentDE()->GetBusPatchId(ii);
        
        Int_t nManus = AliMpDDLStore::Instance(kFALSE)->GetBusPatch(BPid, kFALSE)->GetNofManus();
        for (Int_t jj = 0; jj < nManus; jj++)
	      {
          Int_t manuId = AliMpDDLStore::Instance(kFALSE)->GetBusPatch(BPid, kFALSE)->GetManuId(jj);
          
          Int_t nChannels = DEiter.CurrentDE()->NofChannelsInManu(manuId);
          for (Int_t channel = 0; channel < nChannels; channel++)
          {
            Int_t padStatus = status.PadStatus(DEid, manuId, channel);
            if ((rejectMask == 0 && padStatus != 0) || (padStatus & rejectMask) != 0)
            {
              nBadChannels->SetDetectionElementProbability(DEid, nBadChannels->DetectionElementProbability(DEid) + 1.0);
              nBadChannels->SetBusPatchProbability(BPid, nBadChannels->BusPatchProbability(BPid) + 1.0);
              nBadChannels->SetManuProbability(DEid, manuId, nBadChannels->ManuProbability(DEid, manuId) + 1.0);
              nBadChannels->SetChannelProbability(DEid, manuId, channel, 1.0);
            }
          }
	      }
      }
    }
  
  return nBadChannels;
}


//____________________________________________________________________________________________
Int_t AddEventsSingleRun(Int_t runNumber, Long64_t nEvents, AliMUONRejectList& rejectedEvents)
{
  AliCDBManager::Instance()->SetRun(runNumber);
  
  AliMpCDB::LoadAll();

  Double_t maxBadChannels = 80.0;
  
  AliMUONRecoParam *recoParam = AliMUONCDB::LoadRecoParam();
  AliMUONCalibrationData calibrationData(runNumber);
  AliMUONPadStatusMaker status(calibrationData);
  
  status.SetLimits(*recoParam);
  
  Int_t rejectMask = recoParam->PadGoodnessMask();
  
  // Functions to compute the number of bad channels at the DE, BP and Manu level
  AliMUONRejectList* nBadChannels = CheckDE_BP_ManuPedestals(rejectMask, status);
  
  AliMpDEIterator DEiter;
  for (DEiter.First(); !DEiter.IsDone(); DEiter.Next())
    if (DEiter.CurrentDEId() < 1100)
    {
      Int_t DEid = DEiter.CurrentDEId();
      Int_t nBusPatches = DEiter.CurrentDE()->GetNofBusPatches();
      
      Double_t nChannelsDE = static_cast<Double_t>(AliMpDDLStore::Instance(kFALSE)->GetDetElement(DEid)->NofChannels());
      Double_t nBadChannelsDE = nBadChannels->DetectionElementProbability(DEid);
      Double_t ratioBadChannelsDE = 0.0;
      if (nChannelsDE != 0.0)
        ratioBadChannelsDE = 100.0 * nBadChannelsDE / nChannelsDE;
      
      
      Bool_t goodDE = kTRUE;
      if (ratioBadChannelsDE >= maxBadChannels)
      {
        Double_t newRejectedEventsDE = rejectedEvents.DetectionElementProbability(DEid) + static_cast<Double_t>(nEvents);
        rejectedEvents.SetDetectionElementProbability(DEid, newRejectedEventsDE);
        goodDE = kFALSE;
      }
	  	
      if (goodDE)      
        for (Int_t ii = 0; ii < nBusPatches; ii++)
        {
          Int_t BPid = DEiter.CurrentDE()->GetBusPatchId(ii);
          Int_t nManus = AliMpDDLStore::Instance(kFALSE)->GetBusPatch(BPid, kFALSE)->GetNofManus();
          
          Double_t nChannelsBP = static_cast<Double_t>(nManus * AliMpConstants::ManuNofChannels());
          Double_t nBadChannelsBP = nBadChannels->BusPatchProbability(BPid);
          Double_t ratioBadChannelsBP = 100.0 * nBadChannelsBP / nChannelsBP;
          
          Bool_t goodBP = kTRUE;
          
          if (goodBP)
            if (ratioBadChannelsBP >= maxBadChannels)
            {
              Double_t newRejectedEventsBP = rejectedEvents.BusPatchProbability(BPid) + static_cast<Double_t>(nEvents);
              rejectedEvents.SetBusPatchProbability(BPid, newRejectedEventsBP);
              goodBP = kFALSE;
            }
          
          if (goodBP)
            for (Int_t jj = 0; jj < nManus; jj++)
            {
              Int_t manuId = AliMpDDLStore::Instance(kFALSE)->GetBusPatch(BPid, kFALSE)->GetManuId(jj);
              
              Double_t nChannelsManu = DEiter.CurrentDE()->NofChannelsInManu(manuId);;
              Double_t nBadChannelsManu = nBadChannels->ManuProbability(DEid, manuId);
              Double_t ratioBadChannelsManu = 100.0 * nBadChannelsManu / nChannelsManu;
              
              Bool_t goodManu = kTRUE;
              
              if (ratioBadChannelsManu >= maxBadChannels)
              {
                Double_t newRejectedEventsManu = rejectedEvents.ManuProbability(DEid, manuId) + static_cast<Double_t>(nEvents);
                rejectedEvents.SetManuProbability(DEid, manuId, newRejectedEventsManu);
                goodManu = kFALSE;
              }
              
              
              Int_t nChannels = DEiter.CurrentDE()->NofChannelsInManu(manuId);
              
              if (goodManu)
                for (Int_t channel = 0; channel < nChannels; channel++)
                  if (nBadChannels->ChannelProbability(DEid, manuId, channel) > 0.0)
                  {
                    Double_t newRejectedEventsChannel = rejectedEvents.ChannelProbability(DEid, manuId, channel) + static_cast<Double_t>(nEvents);
                    rejectedEvents.SetChannelProbability(DEid, manuId, channel, newRejectedEventsChannel);
                  } // end of Channels loop
            } // end of MANUs loop	    
        } // end of BPs loop
    } // end of DEs loop
  
  delete nBadChannels;
  
  return nEvents;
}

//______________________________________________________________________________
//
//
// The code below might be used someday to retrieve the number of events
// with at least one muon directly from a runlist, instead of requiring
// as input both the run numbers and the number of events.
//
// This approach though is not without its problems :
//
// - even if we are using tags, as the tags are chunk-based (for the moment at
//   least), it's pretty slow to loop on all the chunks (it is becoming usual
//   to get runs with some 600 chunks or so...)
// - if part of a run is not reconstructed (for whatever reason) with the pass
//   we use to compute the number of events with at least one muon, we might
//   underestimate its weight in the final rejectlist for the next passes (
//   if for the next passes the full run is correctly reconstructed for instance)
// - a muon in the tag = either a trigger track or a tracker track, i.e. we don't
//   know for sure it's a good track...
//
// Given all those reasons, we for the moment use a dumb approach : 
//
// - go to alice-logbook.cern.ch . Select the runs you want and the trigger class
//   you want (most probably CMUS1B), and export the list as (run,nevents)
// - run this macro with this (run,nevents) list
//
//______________________________________________________________________________

#if 0

#include "TChain.h"
#include "AliDetectorTagCuts.h"
#include "AliEventTagCuts.h"
#include "AliLHCTagCuts.h"
#include "AliRunTagCuts.h"
#include "AliTagAnalysis.h"
#include "TGridResult.h"

//______________________________________________________________________________
bool CreateXMLCollectionFromRunList(const char* collectionName, 
                                    const char* runlist,
                                    const char* type,
                                    int passNumber)
{
  /// From a list of run, create a collection with only events containing at least one muon
  /// \param collectionName : output name of the collection (without extension, which will be .xml)
  /// \param runlist : text file containing one integer per line = one run number per line
  /// \param type : files to consider, either ESD or AD
  /// \param passNumber : 1 or 2 most probably (to distinguish which reco pass is used)
  
  if (!gGrid) TGrid::Connect("alien://");
  if (!gGrid) 
  {
    return 0x0;
  }
  
  TString stype(type);
  stype.ToUpper();
  
  if ( stype != "ESD" && stype != "AOD" )
  {
    cout << "Only ESD or AOD type supported" << endl;
    return false;
  }
  
  ifstream in(gSystem->ExpandPathName(runlist));
  Int_t runNumber;
  Int_t ntagfiles(0);
  
  AliTagAnalysis tagAnalysis("ESD");
  
  while ( in >> runNumber ) 
  {
    TGridResult *res = gGrid->Query("/alice/data",
                                    Form("%09d/%ss/pass%d/*%d*/Run%d.Event*.ESD.tag.root",
                                         runNumber,stype.Data(),passNumber,runNumber,runNumber));
    Int_t nFiles = res->GetEntries();
    if (!nFiles) 
    {
      continue;
    }
    
    for (Int_t i = 0; i < nFiles; ++i) 
    {
      TString filename = res->GetKey(i, "turl");
      if(filename == "") continue;
      tagAnalysis.AddTagsFile(filename.Data(),kFALSE);
      ++ntagfiles;
    }
    delete res;
  }
  
  cout << ntagfiles << " tag files added" << endl;
  
  AliRunTagCuts runCuts;
  AliEventTagCuts eventCuts;
  AliLHCTagCuts lhcCuts;
  AliDetectorTagCuts detCuts;
  
  eventCuts.SetNMuonRange(1,99999);
  
  return tagAnalysis.CreateXMLCollection(collectionName,&runCuts,&lhcCuts,&detCuts,&eventCuts);  
  
}

//______________________________________________________________________________
TChain* CreateChainFromXMLCollection(const char* collectionName, const char* type)
{
  /// Create a chain from an XML collection file.
  
  if (!gGrid) TGrid::Connect("alien://");
  if (!gGrid) 
  {
    return 0x0;
  }
  
  TString stype(type);
  stype.ToUpper();
  
  if ( stype != "ESD" && stype != "AOD" )
  {
    cout << "Only ESD or AOD type supported" << endl;
    return 0x0;
  }
  
  AliTagAnalysis tagAnalysis(stype.Data());
  
  return tagAnalysis.CreateChainFromCollection(collectionName,stype=="ESD" ? "esdTree":"aodTree");
}

#endif

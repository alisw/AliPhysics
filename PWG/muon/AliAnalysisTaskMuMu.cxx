#include "AliAnalysisTaskMuMu.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliAnalysisMuonUtility.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODTZERO.h"
#include "AliCentrality.h"
#include "AliCodeTimer.h"
#include "AliCounterCollection.h"
#include "AliESDEvent.h"
#include "AliESDTZERO.h"
#include "AliInputEventHandler.h"
#include "AliLog.h" 
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMergeableCollection.h"
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TH1.h"
#include "TH2.h"
#include "THashList.h"
#include "TList.h"
#include "TMath.h"
#include "TObjString.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TRegexp.h"
#include "TROOT.h"
#include <algorithm>
#include <cassert>
#include "AliAnalysisMuMuBase.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisMuMuCutRegistry.h"
#include "AliAnalysisMuMuCutElement.h"
#include "AliAnalysisMuMuCutCombination.h"
#include <set>

/**
 * \class AliAnalysisTaskMuMu
 * 
 * This class steers the work of one or more sub-analysis deriving from AliAnalysisMuMuBase
 * The output contains an AliHistogramCollection, an AliCounterCollection 
 * and an AliAnalysisMuMuBinning
 * This task must be configured a bit before being used. For instance
 * you can select various event cuts, single muon track cuts and
 *  muon pairs cut, as well as defining various bins (for minv and mean pt
 * histograms) in pt,y,phi etc...
 *
 * Note that it's also possible to disable some (or all) histograms
 * (to save speed/memory), using DisableHistograms() method.
 *
 * For an example of such configuration, \see AddTaskMuMu.C
 */

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMuMu)

//_____________________________________________________________________________
AliAnalysisTaskMuMu::AliAnalysisTaskMuMu()
: AliAnalysisTaskSE("AliAnalysisTaskMuMu"),
fHistogramCollection(0),
fEventCounters(0),
fBinning(0x0),
fCutRegistry(0x0),
fBeamYear(""),
fHistogramToDisable(0x0),
fSubAnalysisVector(0x0)
{
  /// Constructor with a predefined list of triggers to consider
  /// Note that we take ownership of cutRegister
  ///
  
//  fBranchNames = "AOD:header,tracks,vertices,tracklets,AliAODTZERO,AliAODVZERO";

  DefineOutput(1,AliMergeableCollection::Class());
  DefineOutput(2,AliCounterCollection::Class());
  DefineOutput(3,AliAnalysisMuMuBinning::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskMuMu::~AliAnalysisTaskMuMu()
{
  /// dtor

  if (fHistogramCollection && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode())
  {
    delete fHistogramCollection;
  }

  if (fEventCounters && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode())
  {
    delete fEventCounters;
  }

  if (fBinning && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode())
  {
    delete fBinning;
  }

  delete fHistogramToDisable;
  
  delete fCutRegistry;
  
  delete fSubAnalysisVector;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::AdoptSubAnalysis(AliAnalysisMuMuBase* analysis)
{
  if (!fSubAnalysisVector)
  {
    fSubAnalysisVector = new TObjArray;
    fSubAnalysisVector->SetOwner(kTRUE);
  }
  if ( !fSubAnalysisVector->FindObject(analysis) )
  {
    fSubAnalysisVector->Add(analysis);
  }
}

//_____________________________________________________________________________
AliAnalysisMuMuCutRegistry* AliAnalysisTaskMuMu::CutRegistry() const
{
    /// Return (and create if not yet there) our cut registry
  if (!fCutRegistry)
  {
    fCutRegistry = new AliAnalysisMuMuCutRegistry;
  }
  return fCutRegistry;
}

//_____________________________________________________________________________
const char* 
AliAnalysisTaskMuMu::DefaultCentralityName() const
{
  /// Get default centrality name
  if ( !fBeamYear.Contains("pp") ) return "CENTX";
  else return "PP";
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::DisableHistograms(const char* pattern)
{
  /// Disable the histogramming of all the histograms matching the pattern
  
  TIter next(fSubAnalysisVector);
  AliAnalysisMuMuBase* a;
  
  while ( ( a = static_cast<AliAnalysisMuMuBase*>(next()) ) )
  {
    a->DisableHistograms(pattern);
  }
}

//_____________________________________________________________________________
AliVEvent*
AliAnalysisTaskMuMu::Event() const
{
  // some const-dirty-dancing
  return const_cast<AliAnalysisTaskMuMu*>(this)->InputEvent();
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::Fill(const char* eventSelection, const char* triggerClassName)
{
  // Fill one set of histograms (only called for events which pass the eventSelection cut)
  
  TString seventSelection(eventSelection);
  seventSelection.ToLower();
  
  fEventCounters->Count(Form("event:%s/trigger:%s/centrality:%s/run:%d", seventSelection.Data(), triggerClassName, "ALL", fCurrentRunNumber));

  if ( !IsHistogrammingDisabled() )
  {
    TObjArray* centralities = fBinning->CreateBinObjArray("centrality");
    
    TIter next(centralities);
    AliAnalysisMuMuBinning::Range* r;
    
    while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
    {
      TString estimator = r->Quantity();
      
      Float_t fcent = Event()->GetCentrality()->GetCentralityPercentile(estimator.Data());
      if ( fcent < 0.) FillHistos(eventSelection,triggerClassName,"MV0");
      if ( fcent == 0.) FillHistos(eventSelection,triggerClassName,"0V0");
      if ( r->IsInRange(fcent) )
      {
        FillHistos(eventSelection,triggerClassName,r->AsString());
      }
    }
    delete centralities;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillHistos(const char* eventSelection,
                                     const char* triggerClassName,
                                     const char* centrality)
{
  /// Fill histograms for /physics/triggerClassName/centrality
  
  AliCodeTimerAuto("",0);
  
  TIter nextAnalysis(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;
  
  Int_t nTracks = AliAnalysisMuonUtility::GetNTracks(Event());
  
  fEventCounters->Count(Form("event:%s/trigger:%s/centrality:%s/run:%d", eventSelection, triggerClassName, centrality, fCurrentRunNumber));
  
  TIter nextTrackCut(fCutRegistry->GetCutCombinations(AliAnalysisMuMuCutElement::kTrack));
  TIter nextPairCut(fCutRegistry->GetCutCombinations(AliAnalysisMuMuCutElement::kTrackPair));
  
  // loop on single tracks (whatever the type of tracks
  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(nextAnalysis()) ) )
  {
    analysis->DefineHistogramCollection(eventSelection,triggerClassName,centrality);
    
    AliCodeTimerAuto(Form("%s (FillHistosForEvent)",analysis->ClassName()),1);

    if ( MCEvent() != 0x0 )
    {
      analysis->FillHistosForMCEvent(eventSelection,triggerClassName,centrality);
    }

    analysis->FillHistosForEvent(eventSelection,triggerClassName,centrality);
    
    for (Int_t i = 0; i < nTracks; ++i)
    {
      AliVParticle* tracki = AliAnalysisMuonUtility::GetTrack(i,Event());
      
      nextTrackCut.Reset();
      AliAnalysisMuMuCutCombination* trackCut;
      
      while ( ( trackCut = static_cast<AliAnalysisMuMuCutCombination*>(nextTrackCut()) ) )
      {
        if ( trackCut->Pass(*tracki) )
        {
          analysis->FillHistosForTrack(eventSelection,triggerClassName,centrality,trackCut->GetName(),*tracki);
        }
      }
      
      if (!AliAnalysisMuonUtility::IsMuonTrack(tracki) ) continue;
      
      // loop on track pairs (here we only consider muon pairs)
      
      for (Int_t j = i+1; j < nTracks; ++j)
      {
        AliVParticle* trackj = AliAnalysisMuonUtility::GetTrack(j,Event());
        
        if (!AliAnalysisMuonUtility::IsMuonTrack(trackj) ) continue;
        
        nextPairCut.Reset();
        AliAnalysisMuMuCutCombination* pairCut;
        
        while ( ( pairCut = static_cast<AliAnalysisMuMuCutCombination*>(nextPairCut()) ) )
        {
          Bool_t testi = (pairCut->IsTrackCutter()) ? pairCut->Pass(*tracki) : kTRUE;
          Bool_t testj = (pairCut->IsTrackCutter()) ? pairCut->Pass(*trackj) : kTRUE;
          Bool_t testij = pairCut->Pass(*tracki,*trackj);
          
          if ( ( testi || testj ) && testij )
          {
            analysis->FillHistosForPair(eventSelection,triggerClassName,centrality,pairCut->GetName(),*tracki,*trackj);
          }
        }
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FinishTaskOutput()
{
  /// prune empty histograms BEFORE mergin, in order to save some bytes...
  
  if ( fHistogramCollection )
  {
    fHistogramCollection->PruneEmptyObjects();
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::GetSelectedTrigClassesInEvent(const AliVEvent* event, TObjArray& array)
{
  /// Fills the array with a list of TObjString of the trigger classes that the various
  /// cuts accept for this event
  
  array.Clear();
  
  if (!event)
  {
    AliError("Will get a hard time selecting trigger classes with an empty event...");
    return;
  }
  
  TString firedTriggerClasses = event->GetFiredTriggerClasses();
  UInt_t l0 = AliAnalysisMuonUtility::GetL0TriggerInputs(event);
  UInt_t l1 = AliAnalysisMuonUtility::GetL1TriggerInputs(event);
  UInt_t l2 = AliAnalysisMuonUtility::GetL2TriggerInputs(event);

  std::set<std::string> tmpArray;
  
  TIter nextCutCombination(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kTriggerClass));
  AliAnalysisMuMuCutCombination* cutCombination;
  
  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextCutCombination()) ) )
  {
    TString acceptedTriggerClasses;
    
    if ( cutCombination->Pass(firedTriggerClasses,acceptedTriggerClasses,l0,l1,l2) )
    {
      TObjArray* split = acceptedTriggerClasses.Tokenize(" ");
      TIter next(split);
      TObjString* str;
      while ( ( str = static_cast<TObjString*>(next()) ) )
      {
        tmpArray.insert(str->String().Data());
      }
      delete split;
    }
  }
  
  std::set<std::string>::const_iterator it;
  
  for ( it = tmpArray.begin(); it != tmpArray.end(); ++it )
  {
    array.Add(new TObjString(it->c_str()));
  }
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::IsHistogramDisabled(const char* hname) const
{
  /// Whether or not a given histogram (identified by its name)
  /// is disabled or not
  
  TIter next(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;
  
  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(next()) ) )
  {
    if ( analysis->IsHistogramDisabled(hname) )
    {
      return kTRUE;
    }
  }
  
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::IsHistogrammingDisabled() const
{
  /// Whether or not *all* histograms are disabled
  
  Bool_t disabled(kTRUE);
  
  TIter next(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;

  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(next()) ) )
  {
    disabled = disabled && analysis->IsHistogrammingDisabled();
  }

  return disabled;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::IsPP() const
{
  // whether we're dealing with proton proton collisions
  return fBeamYear.Contains("pp");
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::NotifyRun()
{
  /// Called at each change of run 
  
  AliDebug(1,Form("Run %09d File %s",fCurrentRunNumber,CurrentFileName()));
 
  TIter next(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;
  
  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(next()) ) )
  {
    analysis->SetRun(fInputHandler);
  }
}

//_____________________________________________________________________________
void 
AliAnalysisTaskMuMu::Print(Option_t* opt) const
{
  /// Print the definition of this analysis
  
  cout << ClassName() << " - " << GetName() << " - " << fBeamYear.Data() << endl;

  TIter next(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;
  
  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(next()) ) )
  {
    analysis->Print(opt);
  }

  fCutRegistry->Print("ALL");
  
  if ( fBinning )
  {
    cout << "Binning" << endl;
    fBinning->Print();
  }
}

//_____________________________________________________________________________
void
AliAnalysisTaskMuMu::Terminate(Option_t* opt)
{
  /// Called once at the end of the query
  /// Just a simple printout of the stat we analyse and how many histograms
  /// we got
  
  TIter next(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;
  
  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(next()) ) )
  {
    analysis->Terminate(opt);
  }

  fHistogramCollection = dynamic_cast<AliMergeableCollection*>(GetOutputData(1));

  if (!fHistogramCollection)
  {
    AliError("Could not find back histogram collection in output...");
  }
  else
  {
    // Removes empty objects and also the event histos of the Nch task
    fHistogramCollection->PruneEmptyObjects();
    
    UInt_t size2 = fHistogramCollection->EstimateSize();

    TIter nextHistogram(fHistogramCollection->CreateIterator());
    TObject* object;
    
    while ( ( object = nextHistogram() ) )
    {
      if ( object->IsA()->InheritsFrom(TH1::Class()) )
      {
        TH1* h = static_cast<TH1*>(object);
        if ( h->GetXaxis()->GetLabels() )
        {
          h->LabelsDeflate("X");
        }
      }
    }
    
    AliInfo(Form("size after prune histograms = %5.1f MB",size2/1024.0/1024.0));
  
    fHistogramCollection->Print("-");
  }
  
  fEventCounters = dynamic_cast<AliCounterCollection*>(GetOutputData(2));
  
  if (!fEventCounters)
  {
    AliError("Could not find back counters in output...");
  }
  else
  {
    fEventCounters->Print("trigger/event");
  }
  
  // post param container(s)
  PostData(3,fBinning);
}

//_____________________________________________________________________________
AliAnalysisMuMuBinning* AliAnalysisTaskMuMu::Binning() const
{
  // Return our binning (making a default one if not already created
  if ( fBinning ) return fBinning;
  
  fBinning = new AliAnalysisMuMuBinning("BIN");
  
  return fBinning;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::UserExec(Option_t* /*opt*/)
{
  /// Executed at each event
  
//  static Int_t n(0);
//  AliInfo(Form("EVENT %10d Event()=%p MCEvent()=%p",n,Event(),MCEvent()));
//  ++n;
//  
  AliCodeTimerAuto("",0);
  
  Binning(); // insure we have a binning...
  
  //  if ( MCEvent() )
  //  {
  TIter nextAnalysis(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;  
  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(nextAnalysis()) ) )
  {
    if ( MCEvent() ) // Set the MC flag for all analysis (prior to call anything from them
      // (e.g. any trigger class selection that might behave differently for
      // MC and real trigger classes)
    {
      analysis->SetMC();
    }
    analysis->SetEvent(Event(),MCEvent()); // Set the new event properties derived in the analysis
  }
  //  }
  
  
  TString firedTriggerClasses(AliAnalysisMuonUtility::GetFiredTriggerClasses(Event()));
  
  // first loop to count things not associated to a specific trigger
  TIter nextEventCutCombination(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kEvent));
  AliAnalysisMuMuCutCombination* cutCombination;
  
  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextEventCutCombination())))
  {
    if ( cutCombination->Pass(*fInputHandler) )
    {
      fEventCounters->Count(Form("event:%s/trigger:%s/centrality:%s/run:%d", cutCombination->GetName(), "EVERYTHING",  "ALL", fCurrentRunNumber));

      if ( firedTriggerClasses == "" )
      {
        fEventCounters->Count(Form("event:%s/trigger:%s/centrality:%s/run:%d", cutCombination->GetName(), "EMPTY", "ALL", fCurrentRunNumber));
      }
    }
  }

  // second loop to count only the triggers we're interested in
  TObjArray selectedTriggerClasses;

  GetSelectedTrigClassesInEvent(Event(),selectedTriggerClasses);
  
  TIter next(&selectedTriggerClasses);
  TObjString* tname;
//  Bool_t hasSetEventBeenCalled(kFALSE);

  while ( ( tname = static_cast<TObjString*>(next()) ) )
  {
    nextEventCutCombination.Reset();

    while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextEventCutCombination())) )
    {
      if ( cutCombination->Pass(*fInputHandler) )
      {
//        if (!hasSetEventBeenCalled)
//        {
//          TIter nextAnalysis(fSubAnalysisVector);
//          AliAnalysisMuMuBase* analysis;
//          
//          while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(nextAnalysis()) ) )
//          {
//            analysis->SetEvent(Event(),MCEvent());
//          }
//          hasSetEventBeenCalled = kTRUE;
//        }
        Fill(cutCombination->GetName(),tname->String().Data());
      }
    }
  }
  
  // Post output data.
  PostData(1, fHistogramCollection);
  PostData(2, fEventCounters);
  PostData(3, fBinning);
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::UserCreateOutputObjects()
{
  /// Create histograms
  /// Called once
  
  OpenFile(1);
  
  AliInfo(Form("fCutRegistry=%p",fCutRegistry));
  
  if ( fCutRegistry )
  {
    fCutRegistry->Print();
  }
  
  fHistogramCollection = new AliMergeableCollection("OC");

  fEventCounters = new AliCounterCollection("CC");

  // initialize event counters

  TString eventRubric;
  TIter next(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kEvent));
  AliAnalysisMuMuCutCombination* cutCombination;
  
  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(next())) )
  {
    TString cutName = cutCombination->GetName();
    if ( eventRubric.Length() > 0 ) eventRubric += "/";
    eventRubric += cutName;
  }
  
  fEventCounters->AddRubric("event", eventRubric.Data());
  
  fEventCounters->AddRubric("trigger", 100);
  
  fEventCounters->AddRubric("centrality", 100);
  
  fEventCounters->AddRubric("run", 1000000);
  
  // Initialize our subtasks, if any...
  
  TIter nextAnalysis(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;
  
  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(nextAnalysis()) ) )
  {
    analysis->Init(*fEventCounters,*fHistogramCollection,*fBinning,*fCutRegistry);
  }

  // finally end the counters initialization
  fEventCounters->Init();
  
  // Post output data.
  PostData(1,fHistogramCollection);
  PostData(2,fEventCounters);
  PostData(3,fBinning);
}

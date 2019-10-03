#include "AliAnalysisTaskMuMu.h"


#include "AliAnalysisManager.h"
#include "AliMultEstimator.h"
#include "AliAnalysisMuMuBase.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliAnalysisMuMuCutCombination.h"
#include "AliAnalysisMuMuCutElement.h"
#include "AliAnalysisMuMuCutRegistry.h"
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
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMergeableCollection.h"
#include "AliMultSelection.h"
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
#include "TParameter.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TRegexp.h"
#include "TROOT.h"
#include <algorithm>
#include <cassert>
#include <set>
///
/// \ class AliAnalysisTaskMuMu
///
/// This class stores several analysis subclasses based on AliAnalysisMuMuBase.
/// The output contains an AliHistogramCollection and
/// an AliCounterCollection filled with tthe different subclasses.
///
/// \author: L. Aphecetche (Subatech)
///
/// This task must be configured a bit before being used. For instance
/// you can select various event cuts, single muon track cuts and
/// muon pairs cut, as well as defining various bins (for minv and mean pt
/// histograms) in pt,y,phi etc...
///
/// Note that it's also possible to disable some (or all) histograms
/// (to save speed/memory), using DisableHistograms() method. Also, this task can run both on single and mixed events.
/// In the case of mixed event, both the UserExec() and UserExecMix() are executed, but in the later case the AliCounterCollection is not filled.
///
/// For an example of such configuration, \see AddTaskMuMu.C
///

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
fCutRegistryMix(0x0),
fBeamYear(""),
fHistogramToDisable(0x0),
fSubAnalysisVector(0x0),
fCountInBins(kFALSE),
fbinWhat(""),
fbinQuantity(""),
fbinFlavor(""),
fDisableHistoLoop(kFALSE),
fLegacyCentrality(kFALSE),
fPool(0x0),
fMaxPoolSize(0),
fMix(kFALSE)
{
  /// Constructor with a predefined list of triggers to consider
  /// Note that we take ownership of cutRegister
  ///

//  fBranchNames = "AOD:header,tracks,vertices,tracklets,AliAODTZERO,AliAODVZERO";

  // Create the pool
  fPool = new TObjArray(2);
  fPool->SetOwner(kTRUE);

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

  if (fPool) delete fPool;

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
AliAnalysisMuMuCutRegistry* AliAnalysisTaskMuMu::CutRegistryMix() const
{
  /// Return (and create if not yet there) our cut registry
  if (!fCutRegistryMix)
  {
    fCutRegistryMix = new AliAnalysisMuMuCutRegistry;
  }
  return fCutRegistryMix;
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
void AliAnalysisTaskMuMu::SetCountInBins( const char* binWhat, const char* binQuantity, const char* binFlavor, Bool_t disableHistoLoop )
{
  /// fCountInBins serve to add a rubric for bins in the Event counter collection
  /// Bin to count, can be set like in AliAnalysisMuMuBinning class, and has to be the same as one of the binnings we give to the task through this class
  /// Only one kind of binning can be used in the counters, since otherwise the bin integrated counts will not be correct (events counted several times)
  /// ONLY FOR EVENT PROPERTIES !
  ///
  ///  FIXME: make a new protection

  if ( !fCountInBins) {
    fCountInBins      = kTRUE;
    fbinWhat          = binWhat;
    fbinQuantity      = binQuantity;
    fbinFlavor        = binFlavor;
    fDisableHistoLoop = disableHistoLoop;
  }
  else AliFatal("Can't be called twice");
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::SetPoolSize (Int_t size)
{
  /// Create pool and set pool size.
  fMix = kTRUE;
  fMaxPoolSize = size;
}


//_____________________________________________________________________________
float AliAnalysisTaskMuMu::CentralityFromCentrality(const char* estimator) const
{
  /// Estimate Centrality from old centrality framework

  AliCentrality* centrality = Event()->GetCentrality();
  if ( centrality ) return centrality->GetCentralityPercentile(estimator);
  else {
    AliWarning("Did not find Centrality !");
    return -9999.0;
  }
}

//_____________________________________________________________________________
float AliAnalysisTaskMuMu::CentralityFromMultSelection(const char* estimator) const
{
  /// Estimate Centrality from new centrality framework

  AliMultSelection* multSelection = static_cast<AliMultSelection*>(Event()->FindListObject("MultSelection"));
  if ( multSelection ) return multSelection->GetMultiplicityPercentile(estimator);
  else {
    AliWarning("Did not find MultSelection !");
    return -9999.0;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::CreateCentralityPools( const char* poolName ) const
{
  /// Create pool according to binnging

  AliInfo( "Creating pools" );
  TObjArray* centralities = fBinning->CreateBinObjArray("centrality");
  if( !centralities )  return;
  Int_t PoolSize = centralities->GetEntries();

  TObjArray* list = new TObjArray(PoolSize);
  list->SetOwner(kTRUE);
  list->SetName(poolName);
  fPool->Add( list );

  for( Int_t iPool = 0; iPool < PoolSize; ++iPool ){
    TList* listbis = new TList();
    listbis->SetOwner(kTRUE);
    static_cast<TObjArray*>(fPool->FindObject(poolName))->Add( listbis );
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::Fill(const char* eventSelection, const char* triggerClassName)
{
  /// Fill one set of histograms (only called for events which pass the eventSelection cut) for a given trigger/event .

  TString seventSelection(eventSelection);
  seventSelection.ToLower();

  // Fill counter collections (only for UserExec() )
  FillCounters(seventSelection.Data(), triggerClassName, "ALL", fCurrentRunNumber);

  TObjArray* centralities = fBinning->CreateBinObjArray("centrality");

  TIter next(centralities);
  AliAnalysisMuMuBinning::Range* r;

  next.Reset();
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) ){

    Float_t fcent     = -42.0;
    TString estimator = r->Quantity();
    if(estimator.Contains("V0MPLUS05")) estimator ="V0Mplus05";
    if(estimator.Contains("V0MMINUS05")) estimator ="V0Mminus05";
    Bool_t isPP(kFALSE);

    // select centrality
    if ( estimator.CompareTo("pp",TString::kIgnoreCase) == 0 ) isPP = kTRUE;
    else {
      if  (fLegacyCentrality)fcent = CentralityFromCentrality(estimator.Data());
      else fcent                   = CentralityFromMultSelection(estimator.Data());
    }

    // Fill histo
    if ( isPP || r->IsInRange(fcent) ){
      if ( !isPP  && !r->IsInRange(fcent) ) continue;

      FillHistos(eventSelection,triggerClassName,r->AsString(),fcent);

      // FIXME: this filling of global centrality histo is misplaced somehow...
      TH1* hcent = fHistogramCollection->Histo(Form("/%s/%s/V0M/Centrality",eventSelection,triggerClassName));
      if (hcent) hcent->Fill(fcent);
    }
  }
  delete centralities;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillPools(const char* eventSelection, const char* triggerClassName)
{
  /// Fill the different centrality with event pools for each combination of eventSelection/triggerClassName

  TString seventSelection(eventSelection);
  seventSelection.ToLower();

  TObjArray* centralities = fBinning->CreateBinObjArray("centrality");

  TIter next(centralities);
  AliAnalysisMuMuBinning::Range* r;

  next.Reset();
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) ){

    Float_t fcent     = -42.0;
    TString estimator = r->Quantity();
    Bool_t isPP(kFALSE);

    // select centrality
    if ( estimator.CompareTo("pp",TString::kIgnoreCase) == 0 ) isPP = kTRUE;
    else {
      if  (fLegacyCentrality)fcent = CentralityFromCentrality(estimator.Data());
      else fcent                   = CentralityFromMultSelection(estimator.Data());
    }

    // Fill histo
    if ( isPP || r->IsInRange(fcent) ){
      if ( !isPP  && !r->IsInRange(fcent) ) continue;

      FillPoolsWithTracks(eventSelection,triggerClassName,fcent);
    }
  }
  delete centralities;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillHistos(const char* eventSelection,
                                     const char* triggerClassName,
                                     const char* centrality,
                                     Float_t cent)
{
  /// Fill histograms

  // Fill counter collections (only for UserExec() )
  FillCounters( eventSelection, triggerClassName, centrality, fCurrentRunNumber);

  // timer
  AliCodeTimerAuto(Form("/%s/%s/%s",eventSelection,triggerClassName,centrality),0);

  // prepare iterators
  TIter nextAnalysis(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;
  TIter nextTrackCut(fCutRegistry->GetCutCombinations(AliAnalysisMuMuCutElement::kTrack));
  TIter nextPairCut(fCutRegistry->GetCutCombinations(AliAnalysisMuMuCutElement::kTrackPair));

  // Get number of tracks
  Int_t nTracks   = AliAnalysisMuonUtility::GetNTracks(Event());

  // The main part, loop over subanalysis and fill histo
  if ( !IsHistogrammingDisabled() && !fDisableHistoLoop ){

    while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(nextAnalysis()) ) )
    {

      // Create proxy for the Histogram collections
      analysis->DefineHistogramCollection(eventSelection,triggerClassName,centrality,fMix);

      if ( MCEvent() != 0x0 )
      {
        AliCodeTimerAuto(Form("%s (FillHistosForMCEvent)",analysis->ClassName()),1);
        if(!fMix) analysis->FillHistosForMCEvent(eventSelection,triggerClassName,centrality);// Implemented in AliAnalysisMuMuNch and AliAnalysisMuMuMinv at the moment
      }

      AliCodeTimerAuto(Form("%s (FillHistosForEvent)",analysis->ClassName()),1);
      analysis->FillHistosForEvent(eventSelection,triggerClassName,centrality); // Implemented in AliAnalysisMuMuNch at the moment

      // --- Loop on all event tracks ---
      for (Int_t i = 0; i < nTracks; ++i){

        // Get track
        AliVParticle* tracki = AliAnalysisMuonUtility::GetTrack(i,Event());
       if (!AliAnalysisMuonUtility::IsMuonTrack(tracki) ) continue;

        nextTrackCut.Reset();
        AliAnalysisMuMuCutCombination* trackCut;

        // Loop on all track selections and fill histos for track that pass it
        while ( ( trackCut = static_cast<AliAnalysisMuMuCutCombination*>(nextTrackCut()) ) )
        {
          if ( trackCut->Pass(*tracki) )
          {
            AliCodeTimerAuto(Form("%s (FillHistosForTrack)",analysis->ClassName()),2);
            analysis->FillHistosForTrack(eventSelection,triggerClassName,centrality,trackCut->GetName(),*tracki);
          }
        }

        // --- loop on muon track pairs (no mix) ---

        for (Int_t j = i+1; j < nTracks; ++j){
          // Get track
          AliVParticle    * trackj = 0x0;
          trackj = AliAnalysisMuonUtility::GetTrack(j,Event());

          if (!AliAnalysisMuonUtility::IsMuonTrack(trackj) ) continue;

          nextPairCut.Reset();
          AliAnalysisMuMuCutCombination* pairCut;

          // Fill pair histo
          while ( ( pairCut = static_cast<AliAnalysisMuMuCutCombination*>(nextPairCut()) ) )
          {
            // Weither or not the pairs pass the tests
            Bool_t testi  = (pairCut->IsTrackCutter()) ? pairCut->Pass(*tracki) : kTRUE;
            Bool_t testj  = (pairCut->IsTrackCutter()) ? pairCut->Pass(*trackj) : kTRUE;
            Bool_t testij = pairCut->Pass(*tracki,*trackj);

            if ( ( testi && testj ) && testij )
            {
              AliCodeTimerAuto(Form("%s (FillHistosForPair)",analysis->ClassName()),3);
              analysis->FillHistosForPair(eventSelection,triggerClassName,centrality,pairCut->GetName(),*tracki,*trackj,kFALSE);
            }
          }
        }

        // --- mix part ---

        if(!fMix) continue;

        TList* currentPool  =0x0;
        nextPairCut.Reset();
        nextTrackCut.Reset();

        AliAnalysisMuMuCutCombination* pairCut;

        // Loop over pair cut
        while ( ( pairCut = static_cast<AliAnalysisMuMuCutCombination*>(nextPairCut()) ) )
        {
          // Loop over single track cut from mixing configuration
          while ( ( trackCut = static_cast<AliAnalysisMuMuCutCombination*>(nextTrackCut()) ) )
          {
            currentPool = FindPool(cent,Form("%s/%s/%s",eventSelection,triggerClassName,trackCut->GetName()));
            if(!currentPool) continue;

            for (Int_t iTrack2 = 0; iTrack2 < currentPool->GetSize(); ++iTrack2)
            {
              // Get track
              AliVParticle    * trackj = 0x0;
              trackj = static_cast<AliVParticle*>(currentPool->At(iTrack2));

              // Weither or not the pairs pass the tests
              Bool_t testi  = trackCut->Pass(*tracki);
              Bool_t testj  = trackCut->Pass(*trackj);
              Bool_t testij = pairCut->Pass(*tracki,*trackj);

              if ( testij && testi && testj ) analysis->FillHistosForPair(eventSelection,triggerClassName,centrality,pairCut->GetName(),*tracki,*trackj,fMix);
            }
          }
        }
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillPoolsWithTracks(const char* eventSelection,
                                             const char* triggerClassName,
                                             Float_t cent)
{
  /// Fill Pools with event track
  Int_t nTrackRemoved =0;
  TIter nextTrackCut(fCutRegistryMix->GetCutCombinations(AliAnalysisMuMuCutElement::kTrack));
  AliAnalysisMuMuCutCombination* trackCut;

  // Get number of tracks
  Int_t nTracks   = AliAnalysisMuonUtility::GetNTracks(Event());
  TList* currentPool(0x0);

  for (Int_t j = 0; j < nTracks; ++j){

    // Get track
    AliVParticle    * trackj = 0x0;
    trackj = AliAnalysisMuonUtility::GetTrack(j,Event());
    if( !AliAnalysisMuonUtility::IsMuonTrack(trackj) ) continue;

    // Fill pools
    nextTrackCut.Reset();
    while ( ( trackCut = static_cast<AliAnalysisMuMuCutCombination*>(nextTrackCut()) ) ){
      if(!trackCut->Pass(*trackj)) continue;

      TString poolName = Form("%s/%s/%s",eventSelection,triggerClassName,trackCut->GetName());
      if( !FindPool( cent,poolName.Data() ) ) CreateCentralityPools(poolName.Data());
      currentPool = FindPool(cent,poolName.Data());
      currentPool->AddFirst( trackj->Clone() );

      // truncate pool if necessary
      while( currentPool->GetSize() > fMaxPoolSize )
      {
        delete currentPool->Last();
        currentPool->RemoveLast();
        ++nTrackRemoved;
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillCounters(const char* eventSelection, const char* triggerClassName, const char* centrality, Int_t currentRun)
{
  /// Fill the AliCounterCollection for a given event/trigger/centrality/run combination
  /// The binning has to be an already existing event property or one (like i.e. <dNch/dEta>) which we can compute in the SetEvent() method and attach it to the event list
  /// We can generalize this method (if needed), now it is only valid for multiplicity

  AliCodeTimerAuto("",0);

  if( fCountInBins ){
    TParameter<Double_t>* p(0x0);
    TObjArray* bin = fBinning->CreateBinObjArray(fbinWhat.Data(),fbinQuantity.Data(),fbinFlavor.Data());

    TString sfbinQuantity(fbinQuantity);
    TString parToFind("");
    if ( !sfbinQuantity.CompareTo("ntrcorr") )      parToFind = "NtrCorr";
    else if ( !sfbinQuantity.CompareTo("ntr") )     parToFind = "Ntr";
    else if ( !sfbinQuantity.CompareTo("nch") )     parToFind = "Nch";
    else if ( !sfbinQuantity.CompareTo("v0a") )     parToFind = "V0ARaw";
    else if ( !sfbinQuantity.CompareTo("v0acorr") ) parToFind = "V0ACorr";
    else if ( !sfbinQuantity.CompareTo("v0ccorr") ) parToFind = "V0CCorr";
    else if ( !sfbinQuantity.CompareTo("v0mcorr") ) parToFind = "V0MCorr";
    else AliError(Form("%s bin quantity not found in event",sfbinQuantity.Data())); //FIXME: Not all the possible binnings are implemented here

    if ( !bin ) AliError(Form("%s,%s,%s binning does not exist",fbinWhat.Data(),fbinQuantity.Data(),fbinFlavor.Data()));
    else{

      TList* list = static_cast<TList*>(Event()->FindListObject("NCH"));
      if (list){

        Int_t i(-1);
        Bool_t parFound(kFALSE);
        while ( i < list->GetEntries() - 1 && !parFound ){
          i++;
          while ( list->At(i)->IsA() != TParameter<Double_t>::Class()  && i < list->GetEntries() - 1 ) i++;// In case there is a diferent object, just to skip it
          p = static_cast<TParameter<Double_t>*>(list->At(i));
          if ( TString(p->GetName()).Contains(parToFind.Data()) )
            parFound = kTRUE;
        }
      } else AliFatal("No multiplicity info on Event");

      TIter next(bin);
      AliAnalysisMuMuBinning::Range* r;
      while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) ){

        if ( r->IsInRange(p->GetVal()) )
          fEventCounters->Count(Form("event:%s/trigger:%s/centrality:%s/run:%d/bin:%s",eventSelection, triggerClassName,centrality, currentRun,r->AsString().Data()));
      }
      delete bin;
    }
  } else fEventCounters->Count(Form("event:%s/trigger:%s/centrality:%s/run:%d", eventSelection, triggerClassName,  centrality, currentRun));
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FinishTaskOutput()
{
  /// prune empty histograms BEFORE mergin, in order to save some bytes...
  if ( fHistogramCollection ) fHistogramCollection->PruneEmptyObjects();
}

//________________________________________________________________________
TList*  AliAnalysisTaskMuMu::FindPool( Float_t cent, const char* poolName ) const
{

  // define number of pools and boundary
  // in principle one could also use vertex range

  if(!fPool->FindObject(poolName)) return 0x0;

  TObjArray* centralities = fBinning->CreateBinObjArray("centrality");
  TIter next(centralities);
  AliAnalysisMuMuBinning::Range* r;

  next.Reset();
  TList* o;
  Int_t iPool =0;
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    if( r->IsInRange(cent) && fPool->FindObject( poolName ) ){
      o = static_cast<TList*>( fPool->FindObject( poolName ) );
      return dynamic_cast<TList*>( o->At( iPool ) );
    }
    else ++iPool;
  }
  return 0x0;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::GetSelectedTrigClassesInEvent(const AliVEvent* event, TObjArray& array)
{
  /// Fills the array with a list of TObjString of the trigger classes that the various
  /// cuts accept for this event

  array.Clear();

  if (!event){
    AliError("Will get a hard time selecting trigger classes with an empty event...");
    return;
  }

  TString firedTriggerClasses = event->GetFiredTriggerClasses();
  UInt_t l0                   = event->GetHeader()->GetL0TriggerInputs();
  UInt_t l1                   = event->GetHeader()->GetL1TriggerInputs();
  UInt_t l2                   = event->GetHeader()->GetL2TriggerInputs();

  std::set<std::string> tmpArray;

  TIter nextCutCombination(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kTriggerClass));
  AliAnalysisMuMuCutCombination* cutCombination;

  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextCutCombination()) ) ){

    TString acceptedTriggerClasses;

    if ( cutCombination->Pass(firedTriggerClasses,acceptedTriggerClasses,l0,l1,l2) ){

      TObjArray* split = acceptedTriggerClasses.Tokenize(" ");
      TIter next(split);
      TObjString* str;
      while ( ( str = static_cast<TObjString*>(next()) ) ) tmpArray.insert(str->String().Data());

      delete split;
    }
  }

  std::set<std::string>::const_iterator it;

  for ( it = tmpArray.begin(); it != tmpArray.end(); ++it ) array.Add(new TObjString(it->c_str()));
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::GetSelectedTrigClassesInEventMix(const AliVEvent* event, TObjArray& array)
{
  /// Fills the array with a list of TObjString of the trigger classes that the various
  /// cuts accept for this event

  array.Clear();

  if (!event){
    AliError("Will get a hard time selecting trigger classes with an empty event...");
    return;
  }

  TString firedTriggerClasses = event->GetFiredTriggerClasses();
  UInt_t l0                   = event->GetHeader()->GetL0TriggerInputs();
  UInt_t l1                   = event->GetHeader()->GetL1TriggerInputs();
  UInt_t l2                   = event->GetHeader()->GetL2TriggerInputs();

  std::set<std::string> tmpArray;

  TIter nextCutCombination(CutRegistryMix()->GetCutCombinations(AliAnalysisMuMuCutElement::kTriggerClass));
  AliAnalysisMuMuCutCombination* cutCombination;

  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextCutCombination()) ) ){

    TString acceptedTriggerClasses;

    if ( cutCombination->Pass(firedTriggerClasses,acceptedTriggerClasses,l0,l1,l2) ){

      TObjArray* split = acceptedTriggerClasses.Tokenize(" ");
      TIter next(split);
      TObjString* str;
      while ( ( str = static_cast<TObjString*>(next()) ) ) tmpArray.insert(str->String().Data());

      delete split;
    }
  }

  std::set<std::string>::const_iterator it;

  for ( it = tmpArray.begin(); it != tmpArray.end(); ++it ) array.Add(new TObjString(it->c_str()));
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::IsHistogramDisabled(const char* hname) const
{
  /// Whether or not a given histogram (identified by its name)
  /// is disabled or not

  TIter next(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;

  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(next()) ) ){
    if ( analysis->IsHistogramDisabled(hname) ) return kTRUE;
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

  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(next()) ) ) disabled = disabled && analysis->IsHistogrammingDisabled();

  return disabled;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::IsPP() const
{
  /// whether we're dealing with proton proton collisions. Affect the centrality selection.
  return fBeamYear.Contains("pp");
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::NotifyRun()
{
  /// Called at each change of run

  AliDebug(1,Form("Run %09d File %s",fCurrentRunNumber,CurrentFileName()));

  TIter next(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;

  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(next()) ) ) analysis->SetRun(fInputHandler);
}

//_____________________________________________________________________________
void
AliAnalysisTaskMuMu::Print(Option_t* opt) const
{
  /// Print the definition of this analysis

  cout << ClassName() << " - " << GetName() << " - " << fBeamYear.Data() << endl;

  TIter next(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;

  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(next()) ) )analysis->Print(opt);

  fCutRegistry->Print("ALL");

  if ( fBinning ){
    cout << "Binning" << endl;
    fBinning->Print();
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::PrintPools() const
{
  /// Print the definition of this analysis

   if(!fPool || fPool->IsEmpty())return;

  printf("\n --- Centrality pools --- \n\n");
  printf(" -> Number of pools : %d \n",fPool->GetEntries());
  printf(" -------------------------- \n");
  for (int i = 0; i < fPool->GetEntries(); ++i){

    TObjArray* pool = static_cast<TObjArray*>(fPool->At(i));

    printf(" ---> pool n°%i \n\n",i);
    printf(" - name : %s\n",fPool->At(i)->GetName() );
    for (int j = 0; j < pool->GetEntries(); ++j){

      TList* list = static_cast<TList*>(pool->At(j));

      printf(" - number of muons stored in centrality bins n°%d : %d\n",j, list->GetEntries() );
    }
    printf("\n");
  }
}

//_____________________________________________________________________________
void
AliAnalysisTaskMuMu::Terminate(Option_t *)
{
  /// Called once at the end of the query
  /// Just a simple printout of the stat we analyse and how many histograms
  /// we got


  fHistogramCollection = dynamic_cast<AliMergeableCollection*>(GetOutputData(1));

  TIter nextAnalysis(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;

  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(nextAnalysis()) ) ){
    analysis->SetHistogramCollection(fHistogramCollection);
    analysis->Terminate();
  }

  if (!fHistogramCollection) AliError("Could not find back histogram collection in output...");
  else{
    // Removes empty objects
    fHistogramCollection->PruneEmptyObjects();

    UInt_t size2 = fHistogramCollection->EstimateSize();

    TIter nextHistogram(fHistogramCollection->CreateIterator());
    TObject* object;

    while ( ( object = nextHistogram() ) ){
      if ( object->IsA()->InheritsFrom(TH1::Class()) ){
        TH1* h = static_cast<TH1*>(object);
        if ( h->GetXaxis()->GetLabels() ) h->LabelsDeflate("X");
      }
    }

    AliInfo(Form("size after prune histograms = %5.1f MB",size2/1024.0/1024.0));

    fHistogramCollection->Print("-");
  }

  fEventCounters = dynamic_cast<AliCounterCollection*>(GetOutputData(2));

  if (!fEventCounters) AliError("Could not find back counters in output...");
  else fEventCounters->Print("trigger/event","centrality:all");

  // post param container(s)
  PostData(3,fBinning);
}

//_____________________________________________________________________________
AliAnalysisMuMuBinning* AliAnalysisTaskMuMu::Binning() const
{
  /// Return our binning (making a default one if not already created
  if ( fBinning ) return fBinning;

  fBinning = new AliAnalysisMuMuBinning("BIN");

  return fBinning;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::UserExec(Option_t* /*opt*/)
{
  /// Executed at each event

 // static Int_t n(0);
 // AliInfo(Form("EVENT %10d Event()=%p MCEvent()=%p",n,Event(),MCEvent()));
 // ++n;

  AliCodeTimerAuto("",0);

  Binning(); // insure we have a binning...

  TIter nextAnalysis(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;

  // Loop over each subanalysis
  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(nextAnalysis()) ) ) {
    // Set the MC flag for all analysis (prior to call anything from them
    // (e.g. any trigger class selection that might behave differently for
    // MC and real trigger classes)
    if ( MCEvent() ) analysis->SetMC();
    analysis->SetEvent(Event(),MCEvent()); // Set the new event properties derived in the analysis
  }

  TString firedTriggerClasses(Event()->GetFiredTriggerClasses());

  TIter nextEventCutCombination(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kEvent));
  AliAnalysisMuMuCutCombination* cutCombination;

  TIter nextEventCutCombinationMix(CutRegistryMix()->GetCutCombinations(AliAnalysisMuMuCutElement::kEvent));
  AliAnalysisMuMuCutCombination* cutCombinationMix;

  // loop over cut combination on event level. Fill counters
  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextEventCutCombination()))){
    if ( cutCombination->Pass(*fInputHandler) )
    {
      // Fill counters
      FillCounters(cutCombination->GetName(), "EVERYTHING",  "ALL", fCurrentRunNumber);
      // Default counter
      if ( firedTriggerClasses == "" ) FillCounters(cutCombination->GetName(),"EMPTY","ALL",fCurrentRunNumber);
    }
  }

  // loop over trigger selected list and cut combination on event level. Fill histos
  TObjArray selectedTriggerClasses;
  selectedTriggerClasses.SetOwner(kTRUE);

  GetSelectedTrigClassesInEvent(Event(),selectedTriggerClasses);

  TIter next(&selectedTriggerClasses);
  TObjString* tname;

  while ( ( tname = static_cast<TObjString*>(next()) ) ){
    nextEventCutCombination.Reset();

    while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextEventCutCombination())) ){
      if ( cutCombination->Pass(*fInputHandler) ) Fill(cutCombination->GetName(),tname->String().Data());
    }
  }

  if(fMix){

    GetSelectedTrigClassesInEventMix(Event(),selectedTriggerClasses);
    TIter nextmix(&selectedTriggerClasses);
    nextmix.Reset();

    while ( ( tname = static_cast<TObjString*>(nextmix()) ) ){
      nextEventCutCombinationMix.Reset();

      while ( ( cutCombinationMix = static_cast<AliAnalysisMuMuCutCombination*>(nextEventCutCombinationMix())) ){
        if ( cutCombinationMix->Pass(*fInputHandler) ) FillPools(cutCombinationMix->GetName(),tname->String().Data());
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

  AliDebug(1,Form("fCutRegistry=%p",fCutRegistry));

  if ( fCutRegistry ) StdoutToAliDebug(1,fCutRegistry->Print());

  fHistogramCollection = new AliMergeableCollection("OC");
  fEventCounters       = new AliCounterCollection("CC");

  // initialize event counters
  TString eventRubric;
  TIter next(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kEvent));
  AliAnalysisMuMuCutCombination* cutCombination;

  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(next())) ){
    TString cutName = cutCombination->GetName();
    if ( eventRubric.Length() > 0 ) eventRubric += "/";
    eventRubric += cutName;
  }

  fEventCounters->AddRubric("event", eventRubric.Data());
  fEventCounters->AddRubric("trigger", 100);
  fEventCounters->AddRubric("centrality", 100);
  fEventCounters->AddRubric("run", 1000000);

  if ( fCountInBins ) fEventCounters->AddRubric("bin", 1000000);

  // Initialize our subtasks, if any...
  TIter nextAnalysis(fSubAnalysisVector);
  AliAnalysisMuMuBase* analysis;

  while ( ( analysis = static_cast<AliAnalysisMuMuBase*>(nextAnalysis()) ) ) analysis->Init(*fEventCounters,*fHistogramCollection,*fBinning,*fCutRegistry);

  // finally end the counters initialization
  fEventCounters->Init();

  // Post output data.
  PostData(1,fHistogramCollection);
  PostData(2,fEventCounters);
  PostData(3,fBinning);
}


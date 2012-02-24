/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliVAnalysisMuon.cxx 47782 2011-02-24 18:37:31Z martinez $ */

//-----------------------------------------------------------------------------
/// \class AliVAnalysisMuon
/// Base class with utilities for muon analysis
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

#include "AliVAnalysisMuon.h"

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "THashList.h"
#include "TStyle.h"
//#include "TMCProcess.h"
#include "TLorentzVector.h"

// STEER includes
#include "AliInputEventHandler.h"
#include "AliCentrality.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
//#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCounterCollection.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

// CORRFW includes
#include "AliCFGridSparse.h"

// PWG3 includes
#include "AliMergeableCollection.h"
#include "AliMuonTrackCuts.h"
#include "AliMuonPairCuts.h"

/// \cond CLASSIMP
ClassImp(AliVAnalysisMuon) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliVAnalysisMuon::AliVAnalysisMuon() :
  AliAnalysisTaskSE(),
  fMuonTrackCuts(0x0),
  fMuonPairCuts(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fTerminateOptions(0x0),
  fChargeKeys(0x0),
  fSrcKeys(0x0),
  fPhysSelKeys(0x0),
  fTriggerClasses(0x0),
  fCentralityClasses(0x0),
  fEventCounters(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0),
  fSelectedTrigPattern(0x0),
  fRejectedTrigPattern(0x0),
  fSelectedTrigLevel(0x0),
  fOutputPrototypeList(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliVAnalysisMuon::AliVAnalysisMuon(const char *name, const AliMuonTrackCuts& trackCuts, const AliMuonPairCuts& pairCuts) :
  AliAnalysisTaskSE(name),
  fMuonTrackCuts(new AliMuonTrackCuts(trackCuts)),
  fMuonPairCuts(new AliMuonPairCuts(pairCuts)),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fTerminateOptions(0x0),
  fChargeKeys(0x0),
  fSrcKeys(0x0),
  fPhysSelKeys(0x0),
  fTriggerClasses(new THashList()),
  fCentralityClasses(0x0),
  fEventCounters(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0),
  fSelectedTrigPattern(new TObjArray()),
  fRejectedTrigPattern(new TObjArray()),
  fSelectedTrigLevel(new TObjArray()),
  fOutputPrototypeList(0x0)
{
  //
  /// Constructor.
  //
  
  fTriggerClasses->SetOwner();
  InitKeys();
  SetTrigClassPatterns();
  SetCentralityClasses();

  DefineOutput(1, TObjArray::Class());
}

//________________________________________________________________________
AliVAnalysisMuon::AliVAnalysisMuon(const char *name, const AliMuonTrackCuts& trackCuts) :
  AliAnalysisTaskSE(name),
  fMuonTrackCuts(new AliMuonTrackCuts(trackCuts)),
  fMuonPairCuts(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fTerminateOptions(0x0),
  fChargeKeys(0x0),
  fSrcKeys(0x0),
  fPhysSelKeys(0x0),
  fTriggerClasses(new THashList()),
  fCentralityClasses(0x0),
  fEventCounters(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0),
  fSelectedTrigPattern(new TObjArray()),
  fRejectedTrigPattern(new TObjArray()),
  fSelectedTrigLevel(new TObjArray()),
  fOutputPrototypeList(0x0)
{
  //
  /// Constructor.
  //
  
  fTriggerClasses->SetOwner();
  InitKeys();
  SetTrigClassPatterns();
  SetCentralityClasses();
  
  DefineOutput(1, TObjArray::Class());
}


//________________________________________________________________________
AliVAnalysisMuon::AliVAnalysisMuon(const char *name, const AliMuonPairCuts& pairCuts) :
  AliAnalysisTaskSE(name),
  fMuonTrackCuts(0x0),
  fMuonPairCuts(new AliMuonPairCuts(pairCuts)),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fTerminateOptions(0x0),
  fChargeKeys(0x0),
  fSrcKeys(0x0),
  fPhysSelKeys(0x0),
  fTriggerClasses(new THashList()),
  fCentralityClasses(0x0),
  fEventCounters(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0),
  fSelectedTrigPattern(new TObjArray()),
  fRejectedTrigPattern(new TObjArray()),
  fSelectedTrigLevel(new TObjArray()),
  fOutputPrototypeList(0x0)
{
  //
  /// Constructor.
  //
  fTriggerClasses->SetOwner();
  InitKeys();
  SetTrigClassPatterns();
  SetCentralityClasses();
    
  DefineOutput(1, TObjArray::Class());
}


//________________________________________________________________________
AliVAnalysisMuon::~AliVAnalysisMuon()
{
  //
  /// Destructor
  //

  delete fMuonTrackCuts;
  delete fMuonPairCuts;
  delete fTerminateOptions;
  delete fChargeKeys;
  delete fSrcKeys;
  delete fPhysSelKeys;
  delete fTriggerClasses;
  delete fCentralityClasses;
  delete fSelectedTrigPattern;
  delete fRejectedTrigPattern;
  delete fSelectedTrigLevel;
  delete fOutputPrototypeList;


  // For proof: do not delete output containers
  if ( ! AliAnalysisManager::GetAnalysisManager() || ! AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) {
    delete fOutputList;
  }
}

//___________________________________________________________________________
void AliVAnalysisMuon::FinishTaskOutput()
{
  //
  /// Remove empty histograms to reduce the number of histos to be merged
  //


  fMergeableCollection->PruneEmptyObjects();

  TString objectName = "";
  
  // Add stat. info from physics selection
  // (usefull when running on AODs)
  if ( fInputHandler ) {
    for ( Int_t istat=0; istat<2; istat++ ) {
      TString statType = ( istat == 0 ) ? "ALL" : "BIN0";
      TH2* hStat = dynamic_cast<TH2*>(fInputHandler->GetStatistics(statType.Data()));
      if ( hStat ) {
        objectName = Form("%s_%s", hStat->GetName(), GetName());
        TH2* cloneStat = static_cast<TH2*>(hStat->Clone(objectName.Data()));
        cloneStat->SetDirectory(0);
        fOutputList->Add(cloneStat);
      }
      else {
        AliWarning("Stat histogram not available");
        break;
      }
    } // loop on stat type
  }
}


//___________________________________________________________________________
void AliVAnalysisMuon::NotifyRun()
{
  /// Set run number for cuts
  if ( fMuonTrackCuts ) fMuonTrackCuts->SetRun(fCurrentRunNumber);
  if ( fMuonPairCuts ) fMuonPairCuts->SetRun(fCurrentRunNumber);
}

//___________________________________________________________________________
void AliVAnalysisMuon::UserCreateOutputObjects() 
{
  //
  /// Create output objects
  //
  AliInfo(Form("   CreateOutputObjects of task %s\n", GetName()));
  
  fOutputList = new TObjArray();
  fOutputList->SetOwner();

  fEventCounters = new AliCounterCollection("eventCounters");

  if ( ! fCentralityClasses ) SetCentralityClasses();
  TString centralityClasses = "";
  for ( Int_t icent=1; icent<=fCentralityClasses->GetNbins(); ++icent ) {
    if ( ! centralityClasses.IsNull() ) centralityClasses += "/";
    centralityClasses += fCentralityClasses->GetBinLabel(icent);
  }
  fEventCounters->AddRubric("selected", "yes/no");
  fEventCounters->AddRubric("trigger", 100);
  fEventCounters->AddRubric("centrality", centralityClasses);
  fEventCounters->Init();
  fOutputList->Add(fEventCounters);
 
  fMergeableCollection = new AliMergeableCollection("outputObjects");
  fOutputList->Add(fMergeableCollection);

  PostData(1, fOutputList);
  
  MyUserCreateOutputObjects();
}


//________________________________________________________________________
void AliVAnalysisMuon::UserExec(Option_t * /*option*/) 
{
  //
  /// Main loop
  /// Called for each event
  //

  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) 
    fESDEvent = dynamic_cast<AliESDEvent*> (InputEvent());

  if ( ! fAODEvent && ! fESDEvent ) {
    AliError ("AOD or ESD event not found. Nothing done!");
    return;
  }

  Int_t physSel = ( fInputHandler->IsEventSelected() & AliVEvent::kAny ) ? kPhysSelPass : kPhysSelReject;

  //
  // Global event info
  //

  TString firedTrigClasses = ( fAODEvent ) ? fAODEvent->GetFiredTriggerClasses() : fESDEvent->GetFiredTriggerClasses();
  firedTrigClasses.Prepend("ANY ");
  AliDebug(2, Form("Fired classes %s", firedTrigClasses.Data()));
  TObjArray* selectTrigClasses = BuildTriggerClasses(firedTrigClasses);
  if ( selectTrigClasses->GetEntries() == 0 ) {
    delete selectTrigClasses;
    return;
  }

  Double_t centrality = InputEvent()->GetCentrality()->GetCentralityPercentile("V0M");
  Int_t centralityBin = fCentralityClasses->FindBin(centrality);
  TString centralityBinLabel = fCentralityClasses->GetBinLabel(centralityBin);

  TString selKey = ( physSel == kPhysSelPass ) ? "yes" : "no";
  for ( Int_t itrig=0; itrig<selectTrigClasses->GetEntries(); ++itrig ) {
    TString trigName = selectTrigClasses->At(itrig)->GetName();
    fEventCounters->Count(Form("trigger:%s/selected:%s/centrality:%s", trigName.Data(), selKey.Data(), centralityBinLabel.Data()));
  }

  ProcessEvent(fPhysSelKeys->At(physSel)->GetName(), *selectTrigClasses, fCentralityClasses->GetBinLabel(centralityBin));

  delete selectTrigClasses;

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliVAnalysisMuon::Terminate(Option_t *)
{
  //
  /// Draw some histogram at the end.
  //
  
  if ( gROOT->IsBatch() ) return;
    
  fOutputList = dynamic_cast<TObjArray*>(GetOutputData(1));
  if ( ! fOutputList ) return;
  fEventCounters = static_cast<AliCounterCollection*>(fOutputList->FindObject("eventCounters"));
  fMergeableCollection = static_cast<AliMergeableCollection*>(fOutputList->FindObject("outputObjects"));
  
  if ( ! fTerminateOptions ) SetTerminateOptions();
  if ( ! fMergeableCollection ) return;
  AliInfo(Form("Histogram collection size %g MB", fMergeableCollection->EstimateSize()/1024.0/1024.0));
  if ( fTerminateOptions->At(3) ) {
    TString sopt = fTerminateOptions->At(3)->GetName();
    if ( sopt.Contains("verbose") ) fMergeableCollection->Print("*"); 
  }
}



//________________________________________________________________________
Int_t AliVAnalysisMuon::GetNTracks()
{
  //
  /// Return the number of tracks in event
  //
  return ( fAODEvent ) ? fAODEvent->GetNTracks() : fESDEvent->GetNumberOfMuonTracks();
}


//________________________________________________________________________
AliVParticle* AliVAnalysisMuon::GetTrack(Int_t itrack)
{
  //
  /// Get the current track
  //
  AliVParticle* track = 0x0;
  if ( fAODEvent ) track = fAODEvent->GetTrack(itrack);
  else track = fESDEvent->GetMuonTrack(itrack);
  return track;
}

//________________________________________________________________________
Double_t AliVAnalysisMuon::MuonMass2() const
{
  /// A usefull constant
  static Double_t m2 = 1.11636129640000012e-02; // using a constant here as the line below is a problem for CINT...
  return m2;
}

//________________________________________________________________________
TLorentzVector AliVAnalysisMuon::GetTrackPair(AliVParticle* track1, AliVParticle* track2) const
{
  //
  /// Get track pair
  //
  
  AliVParticle* tracks[2] = {track1, track2};
  
  TLorentzVector vec[2];
  for ( Int_t itrack=0; itrack<2; ++itrack ) {
    Double_t trackP = tracks[itrack]->P();
    Double_t energy = TMath::Sqrt(trackP*trackP + MuonMass2());
    vec[itrack].SetPxPyPzE(tracks[itrack]->Px(), tracks[itrack]->Py(), tracks[itrack]->Pz(), energy);
  }
  
  TLorentzVector vecPair = vec[0] + vec[1];
  return vecPair;
}


//________________________________________________________________________
Int_t AliVAnalysisMuon::GetNMCTracks()
{
  //
  /// Return the number of MC tracks in event
  //
  Int_t nMCtracks = 0;
  if ( fMCEvent ) nMCtracks = fMCEvent->GetNumberOfTracks();
  else if ( fAODEvent ) {
    TClonesArray* mcArray = (TClonesArray*)fAODEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if ( mcArray ) nMCtracks = mcArray->GetEntries();
  }
  return nMCtracks;
}

//________________________________________________________________________
AliVParticle* AliVAnalysisMuon::GetMCTrack(Int_t trackLabel)
{
  //
  /// MC information can be provided by the MC input handler
  /// (mostly when analyising ESDs) or can be found inside AODs
  /// This method returns the correct one
  //
  AliVParticle* mcTrack = 0x0;
  if ( fMCEvent ) mcTrack = fMCEvent->GetTrack(trackLabel);
  else if ( fAODEvent ) {
    TClonesArray* mcArray = (TClonesArray*)fAODEvent->FindListObject(AliAODMCParticle::StdBranchName());
    if ( mcArray ) mcTrack =  (AliVParticle*)mcArray->At(trackLabel);
  }
  if ( ! mcTrack ) AliWarning(Form("No track with label %i!", trackLabel));
  return mcTrack;
}

//________________________________________________________________________
Int_t AliVAnalysisMuon::GetMotherIndex(AliVParticle* mcParticle)
{
  //
  /// Return the mother index
  //
  Int_t imother = ( fMCEvent ) ? ((AliMCParticle*)mcParticle)->GetMother() : ((AliAODMCParticle*)mcParticle)->GetMother();
  return imother;
}

//________________________________________________________________________
Int_t AliVAnalysisMuon::GetDaughterIndex(AliVParticle* mcParticle, Int_t idaughter)
{
  //
  /// Return the daughter index
  /// idaughter can be:
  /// 0 -> first daughter
  /// 1 -> last daughter
  //
  if ( idaughter < 0 || idaughter > 1 ) {
    AliError(Form("Requested daughter %i Daughter index can be either 0 (first) or 1 (last)", idaughter));
    return -1;
  }
  
  if ( fMCEvent ) {
    if ( idaughter == 0 ) return ((AliMCParticle*)mcParticle)->GetFirstDaughter();
    else return ((AliMCParticle*)mcParticle)->GetLastDaughter();
  }
  
  return ((AliAODMCParticle*)mcParticle)->GetDaughter(idaughter);
}



//________________________________________________________________________
Bool_t AliVAnalysisMuon::IsMC()
{
  //
  /// Contains MC info
  //
  return ( fMCEvent || ( fAODEvent && fAODEvent->FindListObject(AliAODMCParticle::StdBranchName()) ) );
}


//________________________________________________________________________
Int_t AliVAnalysisMuon::GetParticleType(AliVParticle* track)
{
  //
  /// Get particle type from mathced MC track
  //
  
  Int_t trackSrc = kUnidentified;
  Int_t trackLabel = track->GetLabel();
  if ( trackLabel >= 0 ) {
    AliVParticle* matchedMCTrack = GetMCTrack(trackLabel);
    if ( matchedMCTrack ) trackSrc = RecoTrackMother(matchedMCTrack);
  } // track has MC label
  return trackSrc;
}


//________________________________________________________________________
Int_t AliVAnalysisMuon::RecoTrackMother(AliVParticle* mcParticle)
{
  //
  /// Find track mother from kinematics
  //
  
  Int_t recoPdg = mcParticle->PdgCode();
  
  // Track is not a muon
  if ( TMath::Abs(recoPdg) != 13 ) return kRecoHadron;
  
  Int_t imother = GetMotherIndex(mcParticle);
  
  Int_t den[3] = {100, 1000, 1};
  
  Int_t motherType = kDecayMu;
  while ( imother >= 0 ) {
    AliVParticle* part = GetMCTrack(imother);
    //if ( ! part ) return motherType;
    
    Int_t absPdg = TMath::Abs(part->PdgCode());
    
    Bool_t isPrimary = ( fMCEvent ) ? ( imother < fMCEvent->GetNumberOfPrimaries() ) : ((AliAODMCParticle*)part)->IsPrimary();
    
    if ( isPrimary ) {
      if ( absPdg == 24 ) return kWbosonMu;
      
      for ( Int_t idec=0; idec<3; idec++ ) {
        Int_t flv = (absPdg%100000)/den[idec];
        if ( flv > 0 && flv < 4 ) return kDecayMu;
        else if ( flv == 0 || flv > 5 ) continue;
        else {
          if ( den[idec] == 100 ) motherType = kQuarkoniumMu;
          else if ( flv == 4 ) motherType = kCharmMu;
          else motherType = kBeautyMu;
          break; // break loop on pdg code
          // but continue the global loop to find higher mass HF
        }
      } // loop on pdg code
      if ( absPdg < 10 ) break; // particle loop
    } // is primary
    else {
      if ( part->Zv() < -90. ) {
        // If hadronic process => secondary
        //if ( part->GetUniqueID() == kPHadronic ) {
        return kSecondaryMu;
        //}
      }
    } // is secondary
    
    imother = GetMotherIndex(part);
    
  } // loop on mothers
  
  return motherType;
}

//________________________________________________________________________
Bool_t AliVAnalysisMuon::AddObjectToCollection(TObject* object, Int_t index)
{
  //
  /// Add object to collection
  //
  
  if ( ! fOutputPrototypeList ) {
    fOutputPrototypeList = new TObjArray();
    fOutputPrototypeList->SetOwner();
  }
  if ( fOutputPrototypeList->FindObject(object->GetName() ) ) {
    AliWarning(Form("Object with name %s already in the list", object->GetName()));
    return kFALSE;
  }
  if ( index < 0 ) fOutputPrototypeList->Add(object);
  else fOutputPrototypeList->AddAtAndExpand(object, index);
  
  return kTRUE;
}

//________________________________________________________________________
TObject* AliVAnalysisMuon::GetMergeableObject(TString physSel, TString trigClassName, TString centrality, TString objectName)
{
  //
  /// Get mergeable object
  /// (create collection if necessary)
  //
  
  TString identifier = Form("/%s/%s/%s/", physSel.Data(), trigClassName.Data(), centrality.Data());
  
  TObject* obj = fMergeableCollection->GetObject(identifier.Data(), objectName.Data());
  if ( ! obj ) {
    CreateMergeableObjects(physSel, trigClassName, centrality);
    obj = fMergeableCollection->GetObject(identifier.Data(), objectName.Data());
    AliInfo(Form("Mergeable object collection size %g MB", fMergeableCollection->EstimateSize()/1024.0/1024.0));
  }
  return obj;
}

//________________________________________________________________________
TObject* AliVAnalysisMuon::GetSum(TString physSel, TString trigClassNames, TString centrality, TString objectPattern)
{
  //
  /// Sum objects
  /// - physSel, trigClassNames must be in the form: key1,key2
  /// - centrality must be in the form minValue_maxValue
  /// - objectPattern must be in the form match1&match2&match3,match4
  ///   meaning that the object name must contain match1 and match2 and either one of match3 and match4
  
  if ( ! fMergeableCollection ) return 0x0;
  
  // Get centrality range
  Int_t firstCentrality = 1;
  Int_t lastCentrality = fCentralityClasses->GetNbins();
  
  TObjArray* centralityRange = centrality.Tokenize("_");
  Float_t range[2] = {0., 100.};
  if ( centralityRange->GetEntries() >= 2 ) {
    for ( Int_t irange=0; irange<2; ++irange ) {
      range[irange] = ((TObjString*)centralityRange->At(irange))->GetString().Atof();
    }
    firstCentrality = fCentralityClasses->FindBin(range[0]+0.0001);
    lastCentrality = fCentralityClasses->FindBin(range[1]-0.0001);
  }
  delete centralityRange;
  
  TString sumCentralityString = "";
  for ( Int_t icent=firstCentrality; icent<=lastCentrality; ++icent ) {
    if ( ! sumCentralityString.IsNull() ) sumCentralityString += ",";
    sumCentralityString += fCentralityClasses->GetBinLabel(icent);
  }
  
  objectPattern.ReplaceAll(" ","");
  TObjArray* objPatternList = objectPattern.Tokenize("&");
  TObjArray objPatternMatrix(objPatternList->GetEntries());
  objPatternMatrix.SetOwner();
  for ( Int_t ikey=0; ikey<objPatternList->GetEntries(); ikey++ ) {
    TObjArray* subKeyList = ((TObjString*)objPatternList->At(ikey))->GetString().Tokenize(",");
    objPatternMatrix.AddAt(subKeyList, ikey);
  }
  delete objPatternList;
  

  TObjArray objectNameInCollection;
  objectNameInCollection.SetOwner();
  TObjArray* physSelList = physSel.Tokenize(",");
  TObjArray* trigClassList = trigClassNames.Tokenize(",");
  TObjArray* centralityList = sumCentralityString.Tokenize(",");
  for ( Int_t isel=0; isel<physSelList->GetEntries(); isel++ ) {
    for ( Int_t itrig = 0; itrig<trigClassList->GetEntries(); itrig++ ) {
      for ( Int_t icent=0; icent<centralityList->GetEntries(); icent++ ) {
        TString currId = Form("/%s/%s/%s/", physSelList->At(isel)->GetName(), trigClassList->At(itrig)->GetName(),centralityList->At(icent)->GetName());
        TList* objNameList = fMergeableCollection->CreateListOfObjectNames(currId.Data());
        for ( Int_t iobj=0; iobj<objNameList->GetEntries(); iobj++ ) {
          TString objName = objNameList->At(iobj)->GetName();
          if ( ! objectNameInCollection.FindObject(objName.Data()) ) objectNameInCollection.Add(new TObjString(objName.Data()));
        }
        delete objNameList;
      }
    }
  }
  delete physSelList;
  delete trigClassList;
  delete centralityList;

  TString matchingObjectNames = "";
  for ( Int_t iobj=0; iobj<objectNameInCollection.GetEntries(); iobj++ ) {
    TString objName = objectNameInCollection.At(iobj)->GetName();
    Bool_t matchAnd = kTRUE;
    for ( Int_t ikey=0; ikey<objPatternMatrix.GetEntries(); ikey++ ) {
      TObjArray*  subKeyList = (TObjArray*)objPatternMatrix.At(ikey);
      Bool_t matchOr = kFALSE;
      for ( Int_t isub=0; isub<subKeyList->GetEntries(); isub++ ) {
        TString subKeyString = ((TObjString*)subKeyList->At(isub))->GetString();
        if ( objName.Contains(subKeyString.Data()) ) {
          matchOr = kTRUE;
          break;
        }
      }
      if ( ! matchOr ) {
        matchAnd = kFALSE;
        break;
      }
    }
    if ( ! matchAnd ) continue;
    if ( ! matchingObjectNames.IsNull() ) matchingObjectNames.Append(",");
    matchingObjectNames += objName;
  }

  TString idPattern = Form("/%s/%s/%s/%s", physSel.Data(), trigClassNames.Data(), sumCentralityString.Data(), matchingObjectNames.Data());
  idPattern.ReplaceAll(" ","");
  
  AliDebug(1,Form("Sum pattern %s", idPattern.Data()));
  
  return fMergeableCollection->GetSum(idPattern.Data());
}

//___________________________________________________________________________
void AliVAnalysisMuon::CreateMergeableObjects(TString physSel, TString trigClassName, TString centrality)
{
  TObject* obj = 0x0;
  TString objectName = "";
  TString identifier = Form("/%s/%s/%s/", physSel.Data(), trigClassName.Data(), centrality.Data());
  for ( Int_t iobj=0; iobj<fOutputPrototypeList->GetEntries(); ++iobj ) {
    objectName = fOutputPrototypeList->At(iobj)->GetName();
    obj = fOutputPrototypeList->At(iobj)->Clone(objectName.Data());
    fMergeableCollection->Adopt(identifier, obj);
  } // loop on histos
}


//_______________________________________________________________________
Bool_t AliVAnalysisMuon::SetSparseRange(AliCFGridSparse* gridSparse,
                                        Int_t ivar, TString labelName,
                                        Double_t varMin, Double_t varMax,
                                        TString option)
{
  //
  /// Set range in a smart way.
  /// Allows to select a bin from the label.
  /// Check the bin limits.
  //
  
  option.ToUpper();
  Int_t minVarBin = -1, maxVarBin = -1;
  TAxis* axis = gridSparse->GetAxis(ivar);
  
  if ( ! axis ) {
    printf("Warning: Axis %i not found in %s", ivar, gridSparse->GetName());
    return kFALSE;
  }
  
  if ( ! labelName.IsNull() ) {
    minVarBin = axis->FindBin(labelName.Data());
    maxVarBin = minVarBin;
    if ( minVarBin < 1 ) {
      printf("Warning: %s: label %s not found. Nothing done", gridSparse->GetName(), labelName.Data());
      return kFALSE;
    }
  }
  else if ( option.Contains( "USEBIN" ) ) {
    minVarBin = (Int_t)varMin;
    maxVarBin = (Int_t)varMax;
  }
  else {
    minVarBin = axis->FindBin(varMin);
    maxVarBin = axis->FindBin(varMax);
  }
  
  if ( axis->GetFirst() == minVarBin && axis->GetLast() == maxVarBin ) return kFALSE;
  
  gridSparse->SetRangeUser(ivar, axis->GetBinCenter(minVarBin), axis->GetBinCenter(maxVarBin));
  
  return kTRUE;
}

//________________________________________________________________________
void AliVAnalysisMuon::SetTrigClassPatterns(TString pattern)
{
  /// Set trigger classes
  ///
  /// Classes are filled dynamically according to the pattern
  /// - if name contains ! (without spaces): reject it
  /// - in the matching pattern it is also possible to specify the
  ///   pt cut level associated to the trigger
  /// example:
  /// SetTrigClassPatterns("CMBAC CPBI1MSL:Lpt CPBI1MSH:Hpt !ALLNOTRD")
  /// keeps classes containing CMBAC, CPBI1MSL and CPBI1MSH and not containing ALLNOTRD.
  /// In addition, it knows that the class matching CPBI1MSL requires an Lpt trigger
  /// and the one with CPBI1MSH requires a Hpt trigger.
  /// Hence, in the analysis, the function
  /// TrackPtCutMatchTrigClass(track, "CPBIMSL") returns true if track match Lpt
  /// TrackPtCutMatchTrigClass(track, "CPBIMSL") returns true if track match Hpt
  /// TrackPtCutMatchTrigClass(track, "CMBAC") always returns true
  ///
  /// CAVEAT: if you use an fCFContainer and you want an axis to contain the trigger classes,
  ///         please be sure that each pattern matches only 1 trigger class, or triggers will be mixed up
  ///         when merging different chuncks.

  fSelectedTrigPattern->SetOwner();
  if ( fSelectedTrigPattern->GetEntries() > 0 ) fSelectedTrigPattern->Delete();
  fRejectedTrigPattern->SetOwner();
  if ( fRejectedTrigPattern->GetEntries() > 0 ) fRejectedTrigPattern->Delete();
  fSelectedTrigLevel->SetOwner();
  if ( fSelectedTrigLevel->GetEntries() > 0 ) fSelectedTrigLevel->Delete();

  pattern.ReplaceAll("  "," ");
  pattern.ReplaceAll("! ","!");
  pattern.ReplaceAll(" :",":");

  TObjArray* fullList = pattern.Tokenize(" ");

  for ( Int_t ipat=0; ipat<fullList->GetEntries(); ++ipat ) {
    TString currPattern = fullList->At(ipat)->GetName();
    if ( currPattern.Contains("!") ) {
      currPattern.ReplaceAll("!","");
      fRejectedTrigPattern->AddLast(new TObjString(currPattern));
    }
    else {
      TObjArray* arr = currPattern.Tokenize(":");
      fSelectedTrigPattern->AddLast(new TObjString(arr->At(0)->GetName()));
      TString selTrigLevel = ( arr->At(1) ) ? arr->At(1)->GetName() : "none";
      selTrigLevel.ToUpper();
      fSelectedTrigLevel->AddLast(new TObjString(selTrigLevel));
      delete arr;
    }
  }
  
  delete fullList;
}

//________________________________________________________________________
void AliVAnalysisMuon::SetCentralityClasses(Int_t nCentralityBins, Double_t* centralityBins)
{
  //
  /// Set centrality classes
  //
  Double_t* bins = centralityBins;
  Int_t nbins = nCentralityBins;
  
  Double_t defaultCentralityBins[] = {-5., 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 100., 105.};
  if ( ! centralityBins ) {
    bins = defaultCentralityBins;
    nbins = sizeof(defaultCentralityBins)/sizeof(defaultCentralityBins[0])-1;
  }

  if ( fCentralityClasses ) delete fCentralityClasses;
  fCentralityClasses = new TAxis(nbins, bins);
  TString currClass = "";
  for ( Int_t ibin=1; ibin<=fCentralityClasses->GetNbins(); ++ibin ){
    currClass = Form("%.0f_%.0f",fCentralityClasses->GetBinLowEdge(ibin),fCentralityClasses->GetBinUpEdge(ibin));
    fCentralityClasses->SetBinLabel(ibin, currClass.Data());
  }
}

//________________________________________________________________________
void AliVAnalysisMuon::SetTerminateOptions(TString physSel, TString trigClass, TString centralityRange, TString furtherOpts)
{
  //
  /// Set terminate options
  //
  if ( ! fTerminateOptions ) {
    fTerminateOptions = new TObjArray(4);
    fTerminateOptions->SetOwner();
  }
  fTerminateOptions->AddAt(new TObjString(physSel), 0);
  fTerminateOptions->AddAt(new TObjString(trigClass), 1);
  fTerminateOptions->AddAt(new TObjString(centralityRange),2);
  fTerminateOptions->AddLast(new TObjString(furtherOpts));
}

//________________________________________________________________________
void AliVAnalysisMuon::InitKeys()
{
  TString chargeKeys = "MuMinus MuPlus";
  fChargeKeys = chargeKeys.Tokenize(" ");
  
  TString srcKeys = "CharmMu BeautyMu QuarkoniumMu WbosonMu DecayMu SecondaryMu Hadron Unidentified";
  fSrcKeys = srcKeys.Tokenize(" ");
  
  TString physSelKeys = "PhysSelPass PhysSelReject";
  fPhysSelKeys = physSelKeys.Tokenize(" ");
}

//________________________________________________________________________
TObjArray* AliVAnalysisMuon::BuildTriggerClasses(TString firedTrigClasses)
{
  //
  /// Return the list of trigger classes to be considered
  /// for current event. Update the global list if necessary
  //

  TObjArray* selectedTrigClasses = new TObjArray(0);
  selectedTrigClasses->SetOwner();
  
  TObjArray* firedTrigClassesList = firedTrigClasses.Tokenize(" ");

  for ( Int_t itrig=0; itrig<firedTrigClassesList->GetEntries(); ++itrig ) {
    TString trigName = ((TObjString*)firedTrigClassesList->At(itrig))->GetString();
    Bool_t rejectTrig = kFALSE;
    for ( Int_t ipat=0; ipat<fRejectedTrigPattern->GetEntries(); ++ipat ) {
      if ( trigName.Contains(fRejectedTrigPattern->At(ipat)->GetName() ) ) {
        rejectTrig = kTRUE;
        break;
      }
    } // loop on reject pattern
    if ( rejectTrig ) continue;

    Int_t matchPatternIndex = -1;
    for ( Int_t ipat=0; ipat<fSelectedTrigPattern->GetEntries(); ++ipat ) {
      if ( trigName.Contains(fSelectedTrigPattern->At(ipat)->GetName() ) ) {
        matchPatternIndex = ipat;
        break;
      }
    } // loop on keep pattern
    if ( matchPatternIndex < 0 ) continue;

    selectedTrigClasses->AddLast(new TObjString(trigName));
    if ( fTriggerClasses->FindObject(trigName.Data()) ) continue;
    Int_t trigLevel = 0;
    TString trigLevelString = fSelectedTrigLevel->At(matchPatternIndex)->GetName();
    if ( trigLevelString.Contains("APT") ) trigLevel = 1;
    else if ( trigLevelString.Contains("LPT") ) trigLevel = 2;
    else if ( trigLevelString.Contains("HPT") ) trigLevel = 3;
    AliInfo(Form("Adding %s to considered trigger classes",trigName.Data()));
    TObjString* addTrig = new TObjString(trigName);
    UInt_t uniqueId = trigLevel;
    addTrig->SetUniqueID(uniqueId);
    fTriggerClasses->Add(addTrig);
  } // loop on trigger classes

  delete firedTrigClassesList;

  return selectedTrigClasses;
}


//________________________________________________________________________
Bool_t AliVAnalysisMuon::TrackPtCutMatchTrigClass(AliVParticle* track, TString trigClassName)
{
  /// Check if track passes the trigger pt cut level used in the trigger class
  Int_t matchTrig = ( fAODEvent ) ? ((AliAODTrack*)track)->GetMatchTrigger() : ((AliESDMuonTrack*)track)->GetMatchTrigger();
  Int_t classMatchLevel = GetTrigClassPtCutLevel(trigClassName);
  return matchTrig >= classMatchLevel;
}


//________________________________________________________________________
Int_t AliVAnalysisMuon::GetTrigClassPtCutLevel(TString trigClassName)
{
  /// Get trigger class pt cut level for tracking/trigger matching
  TObject* obj = fTriggerClasses->FindObject(trigClassName.Data());
  if ( ! obj ) {
    AliWarning(Form("Class %s not in the list!", trigClassName.Data()));
    return -1;
  }
  
  return obj->GetUniqueID();
}

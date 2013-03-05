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
#include "TRegexp.h"

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
#include "AliVVertex.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

// CORRFW includes
#include "AliCFGridSparse.h"

// PWG3 includes
#include "AliMergeableCollection.h"
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"
#include "AliMuonPairCuts.h"
#include "AliAnalysisMuonUtility.h"

/// \cond CLASSIMP
ClassImp(AliVAnalysisMuon) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliVAnalysisMuon::AliVAnalysisMuon() :
  AliAnalysisTaskSE(),
  fMuonEventCuts(0x0),
  fMuonTrackCuts(0x0),
  fMuonPairCuts(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fTerminateOptions(0x0),
  fChargeKeys(0x0),
  fSrcKeys(0x0),
  fPhysSelKeys(0x0),
  fWeights(0x0),
  fEventCounters(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0),
  fOutputPrototypeList(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliVAnalysisMuon::AliVAnalysisMuon(const char *name, const AliMuonTrackCuts& trackCuts, const AliMuonPairCuts& pairCuts) :
  AliAnalysisTaskSE(name),
  fMuonEventCuts(new AliMuonEventCuts("stdEventCuts","stdEventCuts")),
  fMuonTrackCuts(new AliMuonTrackCuts(trackCuts)),
  fMuonPairCuts(new AliMuonPairCuts(pairCuts)),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fTerminateOptions(0x0),
  fChargeKeys(0x0),
  fSrcKeys(0x0),
  fPhysSelKeys(0x0),
  fWeights(new THashList()),
  fEventCounters(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0),
  fOutputPrototypeList(0x0)
{
  //
  /// Constructor.
  //
  
  InitKeys();
  SetTrigClassPatterns("");
  SetCentralityClasses();
  fWeights->SetOwner();

  DefineOutput(1, TObjArray::Class());
}

//________________________________________________________________________
AliVAnalysisMuon::AliVAnalysisMuon(const char *name, const AliMuonTrackCuts& trackCuts) :
  AliAnalysisTaskSE(name),
  fMuonEventCuts(new AliMuonEventCuts("stdEventCuts","stdEventCuts")),
  fMuonTrackCuts(new AliMuonTrackCuts(trackCuts)),
  fMuonPairCuts(0x0),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fTerminateOptions(0x0),
  fChargeKeys(0x0),
  fSrcKeys(0x0),
  fPhysSelKeys(0x0),
  fWeights(new THashList()),
  fEventCounters(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0),
  fOutputPrototypeList(0x0)
{
  //
  /// Constructor.
  //
  
  InitKeys();
  SetTrigClassPatterns("");
  SetCentralityClasses();
  fWeights->SetOwner();
  
  DefineOutput(1, TObjArray::Class());
}


//________________________________________________________________________
AliVAnalysisMuon::AliVAnalysisMuon(const char *name, const AliMuonPairCuts& pairCuts) :
  AliAnalysisTaskSE(name),
  fMuonEventCuts(new AliMuonEventCuts("stdEventCuts","stdEventCuts")),
  fMuonTrackCuts(0x0),
  fMuonPairCuts(new AliMuonPairCuts(pairCuts)),
  fESDEvent(0x0),
  fAODEvent(0x0),
  fTerminateOptions(0x0),
  fChargeKeys(0x0),
  fSrcKeys(0x0),
  fPhysSelKeys(0x0),
  fWeights(new THashList()),
  fEventCounters(0x0),
  fMergeableCollection(0x0),
  fOutputList(0x0),
  fOutputPrototypeList(0x0)
{
  //
  /// Constructor.
  //
  InitKeys();
  SetTrigClassPatterns("");
  SetCentralityClasses();
  fWeights->SetOwner();
    
  DefineOutput(1, TObjArray::Class());
}


//________________________________________________________________________
AliVAnalysisMuon::~AliVAnalysisMuon()
{
  //
  /// Destructor
  //

  delete fMuonEventCuts;
  delete fMuonTrackCuts;
  delete fMuonPairCuts;
  delete fTerminateOptions;
  delete fChargeKeys;
  delete fSrcKeys;
  delete fPhysSelKeys;
  delete fWeights;
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
  
  // Add stat. info from physics selection
  // (usefull when running on AODs)
  if ( fInputHandler ) {
    for ( Int_t istat=0; istat<2; istat++ ) {
      TString statType = ( istat == 0 ) ? "ALL" : "BIN0";
      TH2* hStat = dynamic_cast<TH2*>(fInputHandler->GetStatistics(statType.Data()));
      if ( hStat ) {
        TString objectName = Form("%s_%s", hStat->GetName(), GetName());
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
  if ( fMuonTrackCuts ) fMuonTrackCuts->SetRun(fInputHandler);
  if ( fMuonPairCuts ) fMuonPairCuts->SetRun(fInputHandler);
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

  if ( ! GetCentralityClasses() ) SetCentralityClasses();
  TString centralityClasses = "";
  for ( Int_t icent=1; icent<=GetCentralityClasses()->GetNbins(); ++icent ) {
    if ( ! centralityClasses.IsNull() ) centralityClasses += "/";
    centralityClasses += GetCentralityClasses()->GetBinLabel(icent);
  }
  fEventCounters->AddRubric("selected", "yes/no");
  fEventCounters->AddRubric("trigger", 100);
  fEventCounters->AddRubric("centrality", centralityClasses);
  fEventCounters->AddRubric("run", 10000);
  fEventCounters->Init();
  fOutputList->Add(fEventCounters);
 
  fMergeableCollection = new AliMergeableCollection("outputObjects");
  fOutputList->Add(fMergeableCollection);

  PostData(1, fOutputList);
  
  fMuonEventCuts->Print();
  
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
  
  if ( ! fMuonEventCuts->IsSelected(fInputHandler) ) return;

  Int_t physSel = ( fInputHandler->IsEventSelected() & AliVEvent::kAny ) ? kPhysSelPass : kPhysSelReject;

  //
  // Global event info
  //
  const TObjArray* selectTrigClasses = fMuonEventCuts->GetSelectedTrigClassesInEvent(InputEvent());

  Double_t centrality = fMuonEventCuts->GetCentrality(InputEvent());
  Int_t centralityBin = GetCentralityClasses()->FindBin(centrality);
  TString centralityBinLabel = GetCentralityClasses()->GetBinLabel(centralityBin);

  TString selKey = ( physSel == kPhysSelPass ) ? "yes" : "no";
  for ( Int_t itrig=0; itrig<selectTrigClasses->GetEntries(); ++itrig ) {
    TString trigName = selectTrigClasses->At(itrig)->GetName();
    fEventCounters->Count(Form("trigger:%s/selected:%s/centrality:%s/run:%i", trigName.Data(), selKey.Data(), centralityBinLabel.Data(),fCurrentRunNumber));
  }

  ProcessEvent(fPhysSelKeys->At(physSel)->GetName(), *selectTrigClasses, centralityBinLabel);

  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliVAnalysisMuon::Terminate(Option_t *)
{
  //
  /// Draw some histogram at the end.
  //
  
  if ( ! fTerminateOptions ) SetTerminateOptions();
  
  if ( gROOT->IsBatch() ) return;
    
  fOutputList = dynamic_cast<TObjArray*>(GetOutputData(1));
  if ( ! fOutputList ) return;
  fEventCounters = static_cast<AliCounterCollection*>(fOutputList->FindObject("eventCounters"));
  fMergeableCollection = static_cast<AliMergeableCollection*>(fOutputList->FindObject("outputObjects"));
  
  if ( ! fMergeableCollection ) return;
  AliInfo(Form("Mergeable collection size %g MB", fMergeableCollection->EstimateSize()/1024.0/1024.0));
  if ( fTerminateOptions->At(3) ) {
    TString sopt = fTerminateOptions->At(3)->GetName();
    if ( sopt.Contains("verbose") ) fMergeableCollection->Print("*"); 
  }
  SetCentralityClassesFromOutput();
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
    AliVParticle* matchedMCTrack = MCEvent()->GetTrack(trackLabel);
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
  
  Int_t imother = AliAnalysisMuonUtility::GetMotherIndex(mcParticle);
  
  Int_t den[3] = {100, 1000, 1};
  
  Int_t motherType = kDecayMu;
  while ( imother >= 0 ) {
    AliVParticle* part = MCEvent()->GetTrack(imother);
    //if ( ! part ) return motherType;
    
    Int_t absPdg = TMath::Abs(part->PdgCode());
    
    Bool_t isPrimary = AliAnalysisMuonUtility::IsPrimary(part, MCEvent());
    
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
    
    imother = AliAnalysisMuonUtility::GetMotherIndex(part);
    
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
  /// - objectPattern must be in the form match1,match2
  ///   meaning that the object name must contain match1 or match2
  ///   wildcard * is allowed
  
  if ( ! fMergeableCollection ) return 0x0;
  
  // Get centrality range
  Int_t firstCentrality = 1;
  Int_t lastCentrality = GetCentralityClasses()->GetNbins();
  
  TObjArray* centralityRange = centrality.Tokenize("_");
  Float_t range[2] = {0., 100.};
  if ( centralityRange->GetEntries() >= 2 ) {
    for ( Int_t irange=0; irange<2; ++irange ) {
      range[irange] = ((TObjString*)centralityRange->At(irange))->GetString().Atof();
    }
    firstCentrality = GetCentralityClasses()->FindBin(range[0]+0.0001);
    lastCentrality = GetCentralityClasses()->FindBin(range[1]-0.0001);
  }
  delete centralityRange;
  
  TString sumCentralityString = "";
  for ( Int_t icent=firstCentrality; icent<=lastCentrality; ++icent ) {
    if ( ! sumCentralityString.IsNull() ) sumCentralityString += ",";
    sumCentralityString += GetCentralityClasses()->GetBinLabel(icent);
  }
  
//  objectPattern.ReplaceAll(" ","");
//  TObjArray* objPatternList = objectPattern.Tokenize("&");
//  TObjArray objPatternMatrix(objPatternList->GetEntries());
//  objPatternMatrix.SetOwner();
//  for ( Int_t ikey=0; ikey<objPatternList->GetEntries(); ikey++ ) {
//    TObjArray* subKeyList = ((TObjString*)objPatternList->At(ikey))->GetString().Tokenize(",");
//    objPatternMatrix.AddAt(subKeyList, ikey);
//  }
//  delete objPatternList;
  

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
  
  TObjArray* objPatternList = objectPattern.Tokenize(",");

  TString matchingObjectNames = "";
  for ( Int_t iobj=0; iobj<objectNameInCollection.GetEntries(); iobj++ ) {
    TString objName = objectNameInCollection.At(iobj)->GetName();
    for ( Int_t ipat=0; ipat<objPatternList->GetEntries(); ipat++ ) {
      TString currPattern = ((TObjString*)objPatternList->At(ipat))->GetString();
      if ( currPattern.Contains("*") ) {
        if ( ! objName.Contains(TRegexp(currPattern.Data(),kTRUE)) ) continue;
      }
      else if ( objName != currPattern ) continue;

      if ( ! matchingObjectNames.IsNull() ) matchingObjectNames.Append(",");
      matchingObjectNames += objName;
    }
  }
  delete objPatternList;
  
//  for ( Int_t iobj=0; iobj<objectNameInCollection.GetEntries(); iobj++ ) {
//    TString objName = objectNameInCollection.At(iobj)->GetName();
//    Bool_t matchAnd = kTRUE;
//    for ( Int_t ikey=0; ikey<objPatternMatrix.GetEntries(); ikey++ ) {
//      TObjArray*  subKeyList = (TObjArray*)objPatternMatrix.At(ikey);
//      Bool_t matchOr = kFALSE;
//      for ( Int_t isub=0; isub<subKeyList->GetEntries(); isub++ ) {
//        TString subKeyString = ((TObjString*)subKeyList->At(isub))->GetString();
//        if ( subKeyString.Contains("*") ) {
//          if ( objName.Contains(TRegexp(subKeyString.Data())) ) {
//            matchOr = kTRUE;
//            break;
//          }
//        }
//        else if ( objName == subKeyString ) {
//          matchOr = kTRUE;
//          break;
//        }
//      }
//      if ( ! matchOr ) {
//        matchAnd = kFALSE;
//        break;
//      }
//    }
//    if ( ! matchAnd ) continue;
//    if ( ! matchingObjectNames.IsNull() ) matchingObjectNames.Append(",");
//    matchingObjectNames += objName;
//  }

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
  
  // Keep for backward compatibility
  
  return AliAnalysisMuonUtility::SetSparseRange(gridSparse,ivar,labelName,varMin, varMax,option);
}

//________________________________________________________________________
TString AliVAnalysisMuon::GetDefaultTrigClassPatterns() const
{
  /// Get default trigger class patterns
  return fMuonEventCuts->GetDefaultTrigClassPatterns();
}

//________________________________________________________________________
void AliVAnalysisMuon::SetTrigClassPatterns(const TString pattern)
{
  /// Set trigger classes
  TString currPattern = pattern;
  if ( currPattern.IsNull() ) { 
    currPattern = GetDefaultTrigClassPatterns();
    currPattern.Append(",!CMUP*"); // by default do not account for UltraPeripheral events
  }
  fMuonEventCuts->SetTrigClassPatterns(currPattern);
}

//________________________________________________________________________
TList* AliVAnalysisMuon::GetAllSelectedTrigClasses() const
{
  /// Get trigger classes
  return fMuonEventCuts->GetAllSelectedTrigClasses();
}

//________________________________________________________________________
void AliVAnalysisMuon::SetCentralityClasses(Int_t nCentralityBins, Double_t* centralityBins)
{
  //
  /// Set centrality classes
  //
  fMuonEventCuts->SetCentralityClasses(nCentralityBins, centralityBins);
}

//________________________________________________________________________
TAxis* AliVAnalysisMuon::GetCentralityClasses() const
{
  //
  /// Set centrality classes
  //
  return fMuonEventCuts->GetCentralityClasses();
}

//________________________________________________________________________
Bool_t AliVAnalysisMuon::SetCentralityClassesFromOutput()
{
  //
  /// Get axis of centrality classes from output key
  //
  if ( ! fMergeableCollection ) return kFALSE;
  TList* centrKeyList = fMergeableCollection->CreateListOfKeys(2);
  TObjArray centrLimitsList;
  centrLimitsList.SetOwner();
  if ( ! centrKeyList ) return kFALSE;
  for ( Int_t ikey=0; ikey<centrKeyList->GetEntries(); ikey++ ) {
    TString centr = static_cast<TObjString*>(centrKeyList->At(ikey))->GetString();
    TObjArray* array = centr.Tokenize("_");
    for ( Int_t ilim=0; ilim<array->GetEntries(); ilim++ ) {
      TString currLim = static_cast<TObjString*>(array->At(ilim))->GetString();
      if ( ! centrLimitsList.FindObject(currLim.Data()) ) centrLimitsList.Add(new TObjString(currLim));
    }
    delete array;
  }
  delete centrKeyList;
  
  // Get unsorted array
  TArrayD bins(centrLimitsList.GetEntries());
  for ( Int_t ibin=0; ibin<centrLimitsList.GetEntries(); ibin++ ) {
    bins[ibin] = static_cast<TObjString*>(centrLimitsList.At(ibin))->GetString().Atof();
  }
  
  // Sort it
  Int_t index[bins.GetSize()];
  TMath::Sort(bins.GetSize(),bins.GetArray(),index, kFALSE);
  
  TArrayD sortedBins(bins.GetSize());
  for ( Int_t ibin=0; ibin<centrLimitsList.GetEntries(); ibin++ ) {
    sortedBins[ibin] = bins[index[ibin]];
  }
  
  SetCentralityClasses(sortedBins.GetSize()-1, sortedBins.GetArray());
  return kTRUE;
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
  //
  /// Init keys
  //
  TString chargeKeys = "MuMinus MuPlus";
  fChargeKeys = chargeKeys.Tokenize(" ");
  
  TString srcKeys = "CharmMu BeautyMu QuarkoniumMu WbosonMu DecayMu SecondaryMu Hadron Unidentified";
  fSrcKeys = srcKeys.Tokenize(" ");
  
  TString physSelKeys = "PhysSelPass PhysSelReject";
  fPhysSelKeys = physSelKeys.Tokenize(" ");
}


//________________________________________________________________________
void AliVAnalysisMuon::SetWeight ( TObject* wgtObj )
{
  /// Set weight
  if ( fWeights->FindObject(wgtObj->GetName()) ) {
    AliWarning(Form("Weight object %s is already in the list",wgtObj->GetName()));
    return;
  }
  fWeights->Add(wgtObj);
}

//________________________________________________________________________
TObject* AliVAnalysisMuon::GetWeight ( const Char_t* wgtName )
{
  /// Get weight
  return fWeights->FindObject(wgtName);
}



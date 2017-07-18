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

#include "AliMuonEventCuts.h"

#include "THashList.h"
#include "TList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TAxis.h"
#include "TArrayI.h"

#include "AliLog.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"

#include "AliMuonTriggerCombo.h"
#include "AliAnalysisMuonUtility.h"

/// \cond CLASSIMP
ClassImp(AliMuonEventCuts) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliMuonEventCuts::AliMuonEventCuts() :
  AliAnalysisCuts(),
  fPhysicsSelectionMask(0),
  fVertexMinNContributors(0),
  fVertexVzMin(0.),
  fVertexVzMax(0.),
  fCheckMask(0),
  fSelectedTrigPattern(0x0),
  fRejectedTrigPattern(0x0),
  fSelectedTrigCombination(0x0),
  fAllSelectedTrigClasses(0x0),
  fCentralityClasses(0x0),
  fAnalysisUtils(0x0),
  fEventTriggerMask(0),
  fEventL0Inputs(0),
  fEventL1Inputs(0),
  fEventL2Inputs(0),
  fEventPS(0),
  fSelectedTrigClassesInEvent(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliMuonEventCuts::AliMuonEventCuts(const char* name, const char* title ) :
AliAnalysisCuts(name, title),
  fPhysicsSelectionMask(0),
  fVertexMinNContributors(0),
  fVertexVzMin(0.),
  fVertexVzMax(0.),
  fCheckMask(0xFFFF),
  fSelectedTrigPattern(0x0),
  fRejectedTrigPattern(0x0),
  fSelectedTrigCombination(0x0),
  fAllSelectedTrigClasses(0x0),
  fCentralityClasses(0x0),
  fAnalysisUtils(0x0),
  fEventTriggerMask(0),
  fEventL0Inputs(0),
  fEventL1Inputs(0),
  fEventL2Inputs(0),
  fEventPS(0),
  fSelectedTrigClassesInEvent(0x0)
{
  /// Constructor
  SetDefaultParameters();
  SetDefaultFilterMask();
  SetDefaultTrigClassPatterns();
  SetCentralityClasses();
  fAnalysisUtils = new AliAnalysisUtils();
}

//________________________________________________________________________
AliMuonEventCuts::AliMuonEventCuts(const AliMuonEventCuts& obj) :
  AliAnalysisCuts(obj),
  fPhysicsSelectionMask(obj.fPhysicsSelectionMask),
  fVertexMinNContributors(obj.fVertexMinNContributors),
  fVertexVzMin(obj.fVertexVzMin),
  fVertexVzMax(obj.fVertexVzMax),
  fCheckMask(obj.fCheckMask),
  fSelectedTrigPattern(( obj.fSelectedTrigPattern ) ? static_cast<TObjArray*>(obj.fSelectedTrigPattern->Clone() ) : 0x0),
  fRejectedTrigPattern(( obj.fRejectedTrigPattern ) ? static_cast<TObjArray*>(obj.fRejectedTrigPattern->Clone() ) : 0x0),
  fSelectedTrigCombination(( obj.fSelectedTrigCombination ) ? static_cast<TObjArray*>(obj.fSelectedTrigCombination->Clone() ) : 0x0),
  fAllSelectedTrigClasses(( obj.fAllSelectedTrigClasses ) ? static_cast<THashList*>(obj.fAllSelectedTrigClasses->Clone() ) : 0x0),
  fCentralityClasses(( obj.fCentralityClasses ) ? static_cast<TAxis*>(obj.fCentralityClasses->Clone() ) : 0x0),
  fAnalysisUtils(( obj.fAnalysisUtils ) ? static_cast<AliAnalysisUtils*>(obj.fAnalysisUtils->Clone() ) : 0x0),
  fEventTriggerMask(obj.fEventTriggerMask),
  fEventL0Inputs(obj.fEventL0Inputs),
  fEventL1Inputs(obj.fEventL1Inputs),
  fEventL2Inputs(obj.fEventL2Inputs),
  fEventPS(obj.fEventPS),
  fSelectedTrigClassesInEvent(( obj.fSelectedTrigClassesInEvent ) ? static_cast<TObjArray*>(obj.fSelectedTrigClassesInEvent->Clone() ) : 0x0)
{
  /// Copy constructor
}


//________________________________________________________________________
AliMuonEventCuts& AliMuonEventCuts::operator=(const AliMuonEventCuts& obj)
{
  /// Assignment operator
  if ( this != &obj ) {
    AliAnalysisCuts::operator=(obj);
    fPhysicsSelectionMask = obj.fPhysicsSelectionMask;
    fVertexMinNContributors = obj.fVertexMinNContributors,
    fVertexVzMin = obj.fVertexVzMin;
    fVertexVzMax = obj.fVertexVzMax;
    fCheckMask = obj.fCheckMask;
    delete fSelectedTrigPattern;
    fSelectedTrigPattern = ( obj.fSelectedTrigPattern ) ? static_cast<TObjArray*>(obj.fSelectedTrigPattern->Clone() ) : 0x0;
    delete fRejectedTrigPattern;
    fRejectedTrigPattern = ( obj.fRejectedTrigPattern ) ? static_cast<TObjArray*>(obj.fRejectedTrigPattern->Clone() ) : 0x0;
    delete fSelectedTrigCombination;
    fSelectedTrigCombination = ( obj.fSelectedTrigCombination ) ? static_cast<TObjArray*>(obj.fSelectedTrigCombination->Clone() ) : 0x0;
    delete fAllSelectedTrigClasses;
    fAllSelectedTrigClasses = ( obj.fAllSelectedTrigClasses ) ? static_cast<THashList*>(obj.fAllSelectedTrigClasses->Clone() ) : 0x0;
    delete fCentralityClasses;
    fCentralityClasses = ( obj.fCentralityClasses ) ? static_cast<TAxis*>(obj.fCentralityClasses->Clone() ) : 0x0;
    delete fAnalysisUtils;
    fAnalysisUtils = ( obj.fAnalysisUtils ) ? static_cast<AliAnalysisUtils*>(obj.fAnalysisUtils->Clone() ) : 0x0;
    fEventTriggerMask = obj.fEventTriggerMask;
    fEventL0Inputs = obj.fEventL0Inputs;
    fEventL1Inputs = obj.fEventL1Inputs;
    fEventL2Inputs = obj.fEventL2Inputs;
    fEventPS = obj.fEventPS;
    delete fSelectedTrigClassesInEvent;
    fSelectedTrigClassesInEvent = ( obj.fSelectedTrigClassesInEvent ) ? static_cast<TObjArray*>(obj.fSelectedTrigClassesInEvent->Clone() ) : 0x0;
  }
  return *this;
}


//________________________________________________________________________
AliMuonEventCuts::~AliMuonEventCuts()
{
  /// Destructor
  delete fSelectedTrigPattern;
  delete fRejectedTrigPattern;
  delete fSelectedTrigCombination;
  delete fAllSelectedTrigClasses;
  delete fSelectedTrigClassesInEvent;
  delete fCentralityClasses;
  delete fAnalysisUtils;
}

//________________________________________________________________________
Bool_t AliMuonEventCuts::IsSelected( TObject* obj )
{
  /// Track is selected
  UInt_t filterMask = GetFilterMask();
  UInt_t selectionMask = GetSelectionMask(obj);

  AliDebug(1, Form("Is event selected %i  mask 0x%x", ( selectionMask & filterMask ) == filterMask, selectionMask ));

  return ( ( selectionMask & filterMask ) == filterMask );
}


//________________________________________________________________________
UInt_t AliMuonEventCuts::GetSelectionMask( const TObject* obj )
{
  /// Get selection mask

  UInt_t selectionMask = 0;

  UInt_t checkMask = fCheckMask | GetFilterMask();

  const AliInputEventHandler* inputHandler = static_cast<const AliInputEventHandler*> ( obj );

  UInt_t physicsSelection = const_cast<AliInputEventHandler*>(inputHandler)->IsEventSelected();

  if ( checkMask & kPhysicsSelected ) {
    if ( physicsSelection & fPhysicsSelectionMask ) selectionMask |= kPhysicsSelected;
  }

  const AliVEvent* event = inputHandler->GetEvent();

  Double_t centrality = GetCentrality(event);
  if ( centrality >= fCentralityClasses->GetXmin() && centrality <= fCentralityClasses->GetXmax() ) selectionMask |= kSelectedCentrality;

  UpdateEvent(event,physicsSelection);

  if ( fSelectedTrigClassesInEvent->GetEntries() > 0 ) selectionMask |= kSelectedTrig;

  if ( checkMask & kGoodVertex ) {
    const AliVVertex* vertex = event->GetPrimaryVertexSPD();
    if ( vertex->GetNContributors() >= GetVertexMinNContributors() &&
      vertex->GetZ() >= GetVertexVzMin() && vertex->GetZ() <= GetVertexVzMax() ) selectionMask |= kGoodVertex;
  }

  if ( checkMask & kNoPileup ) {
    if ( ! fAnalysisUtils->IsPileUpEvent(const_cast<AliVEvent*>(event)) ) selectionMask |= kNoPileup;
    //  // Uncomment to use settings for pPb
    //  if ( fRejectPileup ) {
    //    Int_t nTracklets = ( event.IsA() == AliESDEvent::Class() ) ? static_cast<AliESDEvent*>(event)->GetMultiplicity()->GetNumberOfTracklets() : static_cast<AliAODEvent*>(event)->GetTracklets()->GetNumberOfTracklets();
    //    Int_t nContrib = ( nTracklets < 40 ) ? 3 : 5;
    //    Double_t dist = 0.8;
    //    Bool_t isPielup = ( event.IsA() == AliESDEvent::Class() ) ? static_cast<AliESDEvent*>(event)->IsPileupFromSPD(nContrib,dist) : static_cast<AliAODEvent*>(event)->IsPileupFromSPD(nContrib,dist);
    //    if ( isPielup ) return;
    //  }
  }

  AliDebug(1, Form("Selection mask 0x%x\n", selectionMask));
  return selectionMask;
}


//________________________________________________________________________
Bool_t AliMuonEventCuts::IsSelected( TList* /* list */)
{
  /// Not implemented
  AliError("Function not implemented: Use IsSelected(TObject*)");
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliMuonEventCuts::UpdateEvent ( const AliVEvent* event, UInt_t physicsSelection )
{
  /// Update the transient data member per event

  UInt_t l0Inputs = event->GetHeader()->GetL0TriggerInputs();
  UInt_t l1Inputs = event->GetHeader()->GetL1TriggerInputs();
  UInt_t l2Inputs = event->GetHeader()->GetL2TriggerInputs();

  if ( fSelectedTrigClassesInEvent && ( fEventTriggerMask == event->GetTriggerMask() ) &&
      ( fEventL0Inputs == l0Inputs ) && ( fEventL1Inputs == l1Inputs ) && ( fEventL2Inputs == l2Inputs ) && ( physicsSelection == fEventPS ) ) return kFALSE;

  BuildTriggerClasses(event->GetFiredTriggerClasses(), l0Inputs, l1Inputs, l2Inputs, physicsSelection);

  fEventTriggerMask = event->GetTriggerMask();
  fEventL0Inputs = l0Inputs;
  fEventL1Inputs = l1Inputs;
  fEventL2Inputs = l2Inputs;
  fEventPS = physicsSelection;

  return kTRUE;
}

//________________________________________________________________________
TString AliMuonEventCuts::GetDefaultTrigClassPatterns () const
{
  /// Get the default patterns
  /// (done in such a way to get the most commonly used muon triggers)
  return "CINT[78]-(B|S|SC)-NOPF-[A-Z]+,C0?M(SL|SH|UL|LL)[78]?-(B|S|SC)-NOPF-[A-Z]+";
}

//________________________________________________________________________
TString AliMuonEventCuts::GetDefaultTrigInputsMap () const
{
  /// Get the default trigger inputs
  ///

  TString trigInputsMap = "0VBA:0,";
  trigInputsMap += "0VBC:1,";
  trigInputsMap += "0SMB:2,";
  trigInputsMap += "0TVX:3,";
  trigInputsMap += "0VGC:4,";
  trigInputsMap += "0VGA:5,";
  trigInputsMap += "0SH1:6,";
  trigInputsMap += "0SH2:7,";
  trigInputsMap += "0HPT:8,";
  trigInputsMap += "0AMU:9,";
  trigInputsMap += "0OB0:10,";
  trigInputsMap += "0ASL:11,";
  trigInputsMap += "0MSL:12,";
  trigInputsMap += "0MSH:13,";
  trigInputsMap += "0MUL:14,";
  trigInputsMap += "0MLL:15,";
  trigInputsMap += "0EMC:16,";
  trigInputsMap += "0PH0:17,";
  trigInputsMap += "0HWU:18,";
  trigInputsMap += "0LSR:19,";
  trigInputsMap += "0T0A:20,";
  trigInputsMap += "0BPA:21,";
  trigInputsMap += "0BPC:22,";
  trigInputsMap += "0T0C:23,";

  trigInputsMap += "1EJE:0,";
  trigInputsMap += "1EGA:1,";
  trigInputsMap += "1EJ2:2,";
  trigInputsMap += "1EG2:3,";
  trigInputsMap += "1PHL:4,";
  trigInputsMap += "1PHM:5,";
  trigInputsMap += "1PHH:6,";
  trigInputsMap += "1HCO:8,";
  trigInputsMap += "1HJT:9,";
  trigInputsMap += "1HSE:10,";
  trigInputsMap += "1DUM:11,";
  trigInputsMap += "1HQU:12,";
  trigInputsMap += "1H14:13,";
  trigInputsMap += "1ZMD:14,";
  trigInputsMap += "1ZMB:16,";
  trigInputsMap += "1ZED:17,";
  trigInputsMap += "1ZAC:18,";
  trigInputsMap += "1EJE:19";

  return trigInputsMap;
}

//________________________________________________________________________
void AliMuonEventCuts::SetDefaultTrigClassPatterns ()
{
  /// Set the default patterns
  /// (done in such a way to get all muon triggers)
  SetTrigClassPatterns(GetDefaultTrigClassPatterns(),GetDefaultTrigInputsMap());
}


//________________________________________________________________________
void AliMuonEventCuts::SetTrigClassPatterns ( TString trigPattern, TString trigInputsMap )
{
  /// Set trigger class patterns or
  /// combinations of trigger classes, trigger inputs or physics selection bits
  /// \param trigPattern String containing a comma separated list of patterns
  /// \param trigInputsMap Map of the trigger inputs
  /// See AliMuonTriggerCombo::Init for details on the accepted patterns/combinations
  /// and on the syntax for the trigger inputs map

  if ( ! fSelectedTrigCombination ) {
    fSelectedTrigCombination = new TObjArray();
    fSelectedTrigCombination->SetOwner();
  }
  else if ( fSelectedTrigCombination->GetEntries() > 0 ) fSelectedTrigCombination->Delete();

  if ( ! fSelectedTrigPattern ) {
    fSelectedTrigPattern = new TObjArray();
    fSelectedTrigPattern->SetOwner();
  }
  else if ( fSelectedTrigPattern->GetEntries() > 0 ) fSelectedTrigPattern->Delete();

  if ( ! fRejectedTrigPattern ) {
    fRejectedTrigPattern = new TObjArray();
    fRejectedTrigPattern->SetOwner();
  }
  else if ( fRejectedTrigPattern->GetEntries() > 0 ) fRejectedTrigPattern->Delete();

  if ( ! fAllSelectedTrigClasses ) {
    fAllSelectedTrigClasses = new THashList();
    fAllSelectedTrigClasses->SetOwner();
  }

  TObjArray tmpCombinationList;

  AliMuonTriggerCombo* trigCombo = 0x0;
  TString pattern(trigPattern);
  pattern.ReplaceAll(" ","");
  TObjArray* fullList = pattern.Tokenize(",");
  TIter next(fullList);
  TObjString* objString = 0x0;
  TString badSyntax = "";
  while ( ( objString = static_cast<TObjString*>(next()) ) ) {
    TObjArray* arr = objString->String().Tokenize(":");
    TString matchPtLevel = ( arr->GetEntries() == 2 ) ? arr->At(1)->GetName() : "";
    trigCombo = new AliMuonTriggerCombo(arr->At(0)->GetName(), trigInputsMap.Data(), matchPtLevel.Data());
    delete arr;

    if ( trigCombo->GetType() == AliMuonTriggerCombo::kBadPattern ) {
      badSyntax += Form(" %s",trigCombo->GetName());
      delete trigCombo;
    }
    // else if ( trigCombo->GetType() == AliMuonTriggerCombo::kRejectPattern ) {
    //   fRejectedTrigPattern->Add(trigCombo);
    //   AliDebug(1,Form("Adding %s to reject pattern",trigCombo->GetName()));
    // }
    else if ( trigCombo->GetType() == AliMuonTriggerCombo::kRegex || ( trigCombo->GetType() == AliMuonTriggerCombo::kComboSimple && trigCombo->HasTriggerClasses() ) ) {
      fSelectedTrigPattern->Add(trigCombo);
      AliDebug(1,Form("Adding %s to match pattern",trigCombo->GetName()));
    }
    else {
      tmpCombinationList.Add(trigCombo);
    }
  }

  TString duplicated = "";
  TIter nextCombo(&tmpCombinationList);
  while ( ( trigCombo = static_cast<AliMuonTriggerCombo*>(nextCombo()) ) ) {
    // Check for combinations that are already accounted for in the expressions
    if ( trigCombo->GetType() == AliMuonTriggerCombo::kComboSimple ) {
      TString trigName = trigCombo->GetName();
      if ( CheckTriggerClassPattern(trigName) ) {
        duplicated += Form("%s ",trigCombo->GetName());
        continue;
      }
    }
    fSelectedTrigCombination->Add(trigCombo);
    AliDebug(1,Form("Adding %s to trigger combination (type %u)",trigCombo->GetName(),trigCombo->GetType()));
  }
  delete fullList;

  if ( ! duplicated.IsNull() )
    AliWarning(Form("Triggers %s already accounted in patterns",duplicated.Data()));
  if ( ! badSyntax.IsNull() )
    AliWarning(Form("%s : illegal expressions. Must be in the form:\n   pattern* => keep class if it contains pattern\n   !pattern* => reject class if it contains pattern\n   class&input = keep class if it satisfies the expression (exact matching required)",badSyntax.Data()));
}


//________________________________________________________________________
TArrayI AliMuonEventCuts::GetTrigClassPtCutLevel ( TString trigClassName ) const
{
  /// Get trigger class pt cut level for tracking/trigger matching

  TArrayI ptCutLevel(2);
  ptCutLevel.Reset();

  AliMuonTriggerCombo* trigCombo = static_cast<AliMuonTriggerCombo*>(fAllSelectedTrigClasses->FindObject(trigClassName.Data()));
  if ( ! trigCombo ) {
    AliWarning(Form("Class %s not in the list!", trigClassName.Data()));
    return -1;
  }

  ptCutLevel[0] = trigCombo->GetTrigMatchLevel();
  if ( trigCombo->IsDimuTrigger()
      ) ptCutLevel[1] = trigCombo->GetTrigMatchLevel();

  AliDebug(3,Form("Class %s ptCutLevel %i %i",trigClassName.Data(),ptCutLevel[0],ptCutLevel[1]));

  return ptCutLevel;
}

//________________________________________________________________________
const TObjArray*
AliMuonEventCuts::GetSelectedTrigClassesInEvent ( const TString& firedTriggerClasses,
                                                 UInt_t l0Inputs, UInt_t l1Inputs,
                                                 UInt_t l2Inputs, UInt_t physicsSelection )
{
  /// Return the selected trigger classes in the fired trigger classes
  /// give also the L0,L1,L2 input bit masks

  BuildTriggerClasses(firedTriggerClasses,l0Inputs,l1Inputs,l2Inputs,physicsSelection);

  return fSelectedTrigClassesInEvent;
}

//________________________________________________________________________
const TObjArray*
AliMuonEventCuts::GetSelectedTrigClassesInEvent ( const TString& firedTriggerClasses,
                                                  UInt_t l0Inputs, UInt_t l1Inputs,
                                                  UInt_t l2Inputs )
{
  /// Return the selected trigger classes in the fired trigger classes
  /// give also the L0,L1,L2 input bit masks

  BuildTriggerClasses(firedTriggerClasses,l0Inputs,l1Inputs,l2Inputs,0xFFFFFFFF);

  return fSelectedTrigClassesInEvent;
}

//________________________________________________________________________
const TObjArray* AliMuonEventCuts::GetSelectedTrigClassesInEvent ( const AliVEvent* event )
{
  /// Return the selected trigger classes in the current event
  UpdateEvent(event,0xFFFFFFFF);
  return fSelectedTrigClassesInEvent;
}

//________________________________________________________________________
const TObjArray* AliMuonEventCuts::GetSelectedTrigClassesInEvent ( const AliInputEventHandler* eventHandler )
{
  /// Return the selected trigger classes in the current event
  const AliVEvent* event = eventHandler->GetEvent();
  UpdateEvent(event,const_cast<AliInputEventHandler*>(eventHandler)->IsEventSelected());
  return fSelectedTrigClassesInEvent;
}


//________________________________________________________________________
void AliMuonEventCuts::BuildTriggerClasses ( TString firedTrigClasses,
                                             UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs,
                                             UInt_t physicsSelection )
{
  //
  /// Return the list of trigger classes to be considered
  /// for current event. Update the global list if necessary
  //

  if ( fSelectedTrigClassesInEvent ) fSelectedTrigClassesInEvent->Delete();
  else {
    fSelectedTrigClassesInEvent = new TObjArray(0);
    fSelectedTrigClassesInEvent->SetOwner();
  }

  TString firedTrigClassesAny = "ANY " + firedTrigClasses;

  if ( fSelectedTrigPattern->GetEntriesFast() > 0 ) {
    TObjArray* firedTrigClassesList = firedTrigClassesAny.Tokenize(" ");

    for ( Int_t itrig=0; itrig<firedTrigClassesList->GetEntriesFast(); itrig++ ) {
      TString trigName = static_cast<TObjString*>(firedTrigClassesList->At(itrig))->GetString();

      AliMuonTriggerCombo* foundCombo = static_cast<AliMuonTriggerCombo*>(fAllSelectedTrigClasses->FindObject(trigName.Data()));
      AliMuonTriggerCombo* matchCombo = 0x0;
      if ( ! foundCombo ) {
        matchCombo = CheckTriggerClassPattern(trigName);
        if ( ! matchCombo ) continue;
      }

      AddToEventSelectedClass ( trigName, foundCombo, matchCombo );
    } // loop on trigger classes

    delete firedTrigClassesList;
  }

  for ( Int_t icomb=0; icomb<fSelectedTrigCombination->GetEntriesFast(); icomb++ ) {
    AliMuonTriggerCombo* trigCombo = static_cast<AliMuonTriggerCombo*>(fSelectedTrigCombination->At(icomb));
    if ( trigCombo->MatchEvent(firedTrigClassesAny, l0Inputs, l1Inputs, l2Inputs,physicsSelection) ) {
      AliMuonTriggerCombo* foundCombo = static_cast<AliMuonTriggerCombo*>(fAllSelectedTrigClasses->FindObject(trigCombo->GetName()));
      AddToEventSelectedClass ( trigCombo->GetName(), foundCombo, trigCombo );
    }
  }
}

//_____________________________________________________________________________
AliMuonTriggerCombo*
AliMuonEventCuts::CheckTriggerClassPattern ( const TString& toCheck ) const
{
  // Check if the "toCheck" class matches the user pattern

  TIter nextMatch(fSelectedTrigPattern);
  AliMuonTriggerCombo* matchCombo = 0x0;
  while ( (matchCombo = static_cast<AliMuonTriggerCombo*>(nextMatch())) ) {
    if ( matchCombo->MatchEvent(toCheck,0,0,0,0) ) {
      AliMuonTriggerCombo* rejectCombo = 0x0;
      TIter nextReject(fRejectedTrigPattern);
      while ( (rejectCombo = static_cast<AliMuonTriggerCombo*>(nextReject())) ) {
        if ( ! rejectCombo->MatchEvent(toCheck,0,0,0,0) ) return 0x0;
      }
      return matchCombo;
    }
  }

  return 0x0;
}


//_____________________________________________________________________________
void
AliMuonEventCuts::AddToEventSelectedClass ( const TString& toCheck, const AliMuonTriggerCombo* foundCombo, const AliMuonTriggerCombo* matchCombo )
{
  /// Add current trigger to the selected class for the event

  fSelectedTrigClassesInEvent->Add(new TObjString(toCheck));

  if ( foundCombo ) return;

  AliMuonTriggerCombo* addCombo = 0x0;
  if ( matchCombo->GetType() == AliMuonTriggerCombo::kRegex ) {
    addCombo = new AliMuonTriggerCombo(toCheck.Data(),"",matchCombo->GetTitle());
  }
  else addCombo = static_cast<AliMuonTriggerCombo*>(matchCombo->Clone());
  fAllSelectedTrigClasses->Add(addCombo);
  TString trigLevelInfo = Form("trig level %i ", addCombo->GetTrigMatchLevel());
  trigLevelInfo += addCombo->IsDimuTrigger() ? "di-muon" : "single-muon";
  AliInfo(Form("Adding %s (%s) to considered trigger classes",toCheck.Data(),trigLevelInfo.Data()));
}

//________________________________________________________________________
void AliMuonEventCuts::SetCentralityClasses(Int_t nCentralityBins, Double_t* centralityBins)
{
  //
  /// Set centrality classes
  //
  Double_t* bins = centralityBins;
  Int_t nbins = nCentralityBins;

  Double_t defaultCentralityBins[] = {-5., 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 105.};
  if ( ! centralityBins ) {
    bins = defaultCentralityBins;
    nbins = sizeof(defaultCentralityBins)/sizeof(defaultCentralityBins[0])-1;
  }

  TString centralityEstimator = "V0M";
  if ( fCentralityClasses ) {
    centralityEstimator = GetCentralityEstimator();
    delete fCentralityClasses;
  }
  fCentralityClasses = new TAxis(nbins, bins);
  TString currClass = "";
  for ( Int_t ibin=1; ibin<=fCentralityClasses->GetNbins(); ++ibin ){
    currClass = Form("%.0f_%.0f",fCentralityClasses->GetBinLowEdge(ibin),fCentralityClasses->GetBinUpEdge(ibin));
    fCentralityClasses->SetBinLabel(ibin, currClass.Data());
  }

  SetCentralityEstimator(centralityEstimator);
}

//________________________________________________________________________
void AliMuonEventCuts::SetCentralityEstimator ( TString centralityEstimator )
{
  /// Set centrality estimator
  fCentralityClasses->SetName(centralityEstimator.Data());
}


//________________________________________________________________________
TString AliMuonEventCuts::GetCentralityEstimator () const
{
  /// Get centrality estimator
  return fCentralityClasses->GetName();
}

//________________________________________________________________________
Double_t AliMuonEventCuts::GetCentrality ( const AliVEvent* event ) const
{
  /// Get centrality
  AliVEvent* evt = const_cast<AliVEvent*>(event);
  AliMultSelection* multSelection = static_cast<AliMultSelection*>(evt->FindListObject("MultSelection"));
  return ( multSelection ) ? multSelection->GetMultiplicityPercentile(GetCentralityEstimator()) : evt->GetCentrality()->GetCentralityPercentile(GetCentralityEstimator());
}

//________________________________________________________________________
void AliMuonEventCuts::SetDefaultParameters ()
{
  /// Standard parameters for muon event
  SetPhysicsSelectionMask(AliVEvent::kAny);
  SetVertexMinNContributors(1);
  SetVertexVzLimits();
}


//________________________________________________________________________
void AliMuonEventCuts::SetDefaultFilterMask ()
{
  /// Standard cuts for muon event
  SetFilterMask ( kPhysicsSelected | kSelectedTrig | kGoodVertex );
}

//________________________________________________________________________
void AliMuonEventCuts::Print(Option_t* option) const
{
  //
  /// Print info
  //
  TString sopt(option);
  sopt.ToLower();
  if ( sopt.IsNull() || sopt.Contains("*") || sopt.Contains("all") ) sopt = "mask param";
  UInt_t filterMask = GetFilterMask();
  if ( sopt.Contains("mask") ) {
    printf(" *** Muon event filter mask: *** \n");
    printf("  0x%x\n", filterMask);
    if ( filterMask & kPhysicsSelected ) printf("  Pass physics selection 0x%x\n", fPhysicsSelectionMask);
    if ( filterMask & kSelectedCentrality ) printf(  "%g < centrality (%s) < %g\n", fCentralityClasses->GetXmin(), GetCentralityEstimator().Data(), fCentralityClasses->GetXmax() );
    if ( filterMask & kSelectedTrig ) printf("  Has selected trigger classes\n");
    if ( filterMask & kGoodVertex ) printf("  SPD vertex with %i contributors && %g < Vz < %g\n", GetVertexMinNContributors(), GetVertexVzMin(), GetVertexVzMax());
    if ( filterMask & kNoPileup ) printf("  Reject pileup with SPD\n");
    printf(" ******************** \n");
  }
  if ( sopt.Contains("param") ) {
    printf(" *** Muon event parameters: *** \n");
    printf("  Centrality estimator: %s\n", GetCentralityEstimator().Data());
    printf(" ******************** \n");
  }
}

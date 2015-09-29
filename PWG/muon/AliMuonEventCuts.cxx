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

#include "TMath.h"
#include "TFormula.h"
#include "THashList.h"
#include "TList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "TParameter.h"
#include "TKey.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TArrayI.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TDataMember.h"

#include "AliLog.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliCentrality.h"
#include "AliAnalysisUtils.h"

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
  fSelectedTrigLevel(0x0),
  fSelectedTrigCombination(0x0),
  fTrigInputsMap(0x0),
  fPhysSelBits(0x0),
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
  fSelectedTrigPattern(new TObjArray()),
  fRejectedTrigPattern(new TObjArray()),
  fSelectedTrigLevel(new TObjArray()),
  fSelectedTrigCombination(new TObjArray()),
  fTrigInputsMap(new THashList()),
  fPhysSelBits(new THashList()),
  fAllSelectedTrigClasses(new THashList()),
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
  SetTrigClassLevels();
  SetCentralityClasses();
  fAnalysisUtils = new AliAnalysisUtils();
  fAllSelectedTrigClasses->SetOwner();
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
  fSelectedTrigLevel(( obj.fSelectedTrigLevel ) ? static_cast<TObjArray*>(obj.fSelectedTrigLevel->Clone() ) : 0x0),
  fSelectedTrigCombination(( obj.fSelectedTrigCombination ) ? static_cast<TObjArray*>(obj.fSelectedTrigCombination->Clone() ) : 0x0),
  fTrigInputsMap(( obj.fTrigInputsMap ) ? static_cast<THashList*>(obj.fTrigInputsMap->Clone() ) : 0x0),
  fPhysSelBits(( obj.fPhysSelBits ) ? static_cast<THashList*>(obj.fPhysSelBits->Clone()) : 0x0 ),
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
    delete fSelectedTrigLevel;
    fSelectedTrigLevel = ( obj.fSelectedTrigLevel ) ? static_cast<TObjArray*>(obj.fSelectedTrigLevel->Clone() ) : 0x0;
    delete fSelectedTrigCombination;
    fSelectedTrigCombination = ( obj.fSelectedTrigCombination ) ? static_cast<TObjArray*>(obj.fSelectedTrigCombination->Clone() ) : 0x0;
    delete fTrigInputsMap;
    fTrigInputsMap = ( obj.fTrigInputsMap ) ? static_cast<THashList*>(obj.fTrigInputsMap->Clone() ) : 0x0;
    delete fPhysSelBits;
    fPhysSelBits = ( obj.fPhysSelBits ) ? static_cast<THashList*>(obj.fPhysSelBits->Clone()) : 0x0;
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
  delete fSelectedTrigLevel;
  delete fSelectedTrigCombination;
  delete fTrigInputsMap;
  delete fPhysSelBits;
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
  /// (done in such a way to get all muon triggers)
  return "CM*,C0M*,CINT*,CPBI*,CCENT*,CV*,!*ABCE*,!*-ACE-*,!*-AC-*,!*-E-*,!*WU*,!*EGA*,!*EJE*,!*PHS*";
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
  /// Set trigger classes
  ///
  /// 1) specify trigger class pattern and reject pattern (wildcard * accepted)
  /// Classes will be filled dynamycally according to the pattern
  /// - if name contains ! (without spaces): reject it
  /// - otherwise keep it
  /// e.g. CMBAC*,!*ALLNOTRD*
  /// keeps classes beginning with CMBAC, and not containing ALLNOTRD.
  ///
  /// CAVEATs:
  ///   a ) if a wildcard is not specified, exact match is required
  ///   b ) if you use an fCFContainer and you want an axis to contain the trigger classes,
  ///       please be sure that each pattern matches only 1 trigger class, or triggers will be mixed up
  ///       when merging different chuncks.
  ///
  ///
  /// 2) specify a combination of triggers (or physics selection bits)
  /// combined through a logical AND "&" or a logical OR "|" (wildcard * NOT accepted)
  /// It is also possible to ask for a trigger class containing a specific trigger input:
  /// e.g. CMSL7-B-NOPF-MUON&0MSH,CMSL7-B-NOPF-MUON,CMSL7-B-NOPF-MUON|CMSL8-B-NOPF-MUON
  /// will give the events with:
  /// - the single low trigger class fired and containing a single high trigger input
  /// - the single low trigger class fired
  /// - the single low trigger class 7 or 8 fired
  /// Also, kMUU7|kMUL7 will give events firing either AliVEvent::kMUU7 or AliVEvent::kMUL7
  /// By default, when specific trigger combinations are provided, the most general case
  /// based on trigger pattern is disabled...but it can be activated with the disableTrigPattern flag
  ///
  /// example:
  /// SetTrigClassPatterns("CINT7*,!*-ACE-*,CMSL7-B-NOPF-MUON,CMSL7-B-NOPF-MUON&0MSH")
  
  
  fSelectedTrigCombination->SetOwner();
  if ( fSelectedTrigCombination->GetEntries() > 0 ) fSelectedTrigCombination->Delete();
  fSelectedTrigPattern->SetOwner();
  if ( fSelectedTrigPattern->GetEntries() > 0 ) fSelectedTrigPattern->Delete();
  fRejectedTrigPattern->SetOwner();
  if ( fRejectedTrigPattern->GetEntries() > 0 ) fRejectedTrigPattern->Delete();

  SetTrigInputsMap(trigInputsMap);
  SetPhysSelBits();

  TString badSyntax = "", duplicated = "";
  TString listName[kNtypes];
  listName[kL0Input] = "L0";
  listName[kL1Input] = "L1";
  listName[kL2Input] = "L2";
  listName[kPhysSelBit] = "PhysSelBit";
  listName[kTrigClass] = "trigClass";
  
  TString pattern(trigPattern);
  pattern.ReplaceAll(" ","");
  TObjArray* fullList = pattern.Tokenize(",");
  TIter next(fullList);
  TObjString* objString = 0x0;
  
  TObjArray combinationList;
  // First search for patterns
  while ( ( objString = static_cast<TObjString*>(next()) ) ) {
    TString currPattern = objString->String();
    Bool_t isCombination = ( currPattern.Contains("&") || currPattern.Contains("|") );
    Bool_t isSingleTrigger = ( ! isCombination && ! currPattern.BeginsWith("0") && ! currPattern.BeginsWith("1") && ! currPattern.BeginsWith("2") & ! currPattern.BeginsWith("k") );
    Bool_t isMatchPattern = ( currPattern.Contains("*") || isSingleTrigger );
    Bool_t isRejectPattern = kFALSE;
    if ( isMatchPattern && currPattern.Contains("!") ) {
      currPattern.ReplaceAll("!","");
      isRejectPattern = kTRUE;
    }
    if ( isCombination && ( isMatchPattern || isRejectPattern ) ) {
      badSyntax += Form(" %s;", currPattern.Data());
      continue;
    }
    if ( isRejectPattern ) {
      fRejectedTrigPattern->AddLast(new TObjString(currPattern));
      AliDebug(1,Form("Adding %s to reject pattern",currPattern.Data()));
    }
    else if ( isMatchPattern ) {
      fSelectedTrigPattern->AddLast(new TObjString(currPattern));
      AliDebug(1,Form("Adding %s to match pattern",currPattern.Data()));
    }
    else combinationList.Add(objString);
  }
  
  // Then check for combinations
  TIter nextCombo(&combinationList);
  while ( ( objString = static_cast<TObjString*>(nextCombo()) ) ) {
    TString currPattern = objString->String();
    
    TString tn (currPattern);
    Bool_t hasAND = kFALSE, hasOR = kFALSE, hasNOT = kFALSE;
    if ( tn.Contains("&") ) {
      tn.ReplaceAll("&",":");
      hasAND = kTRUE;
    }
    if ( tn.Contains("|") ) {
      tn.ReplaceAll("|",":");
      hasOR = kTRUE;
    }
    if ( tn.Contains("!") ) {
      tn.ReplaceAll("!","");
      hasNOT = kTRUE;
    }
    if ( tn.Contains("(") || tn.Contains(")") ) {
      tn.ReplaceAll("(","");
      tn.ReplaceAll(")","");
    }
    
    if ( ! hasAND && ! hasOR && ! hasNOT ) {
      if ( CheckTriggerClassPattern(currPattern) ) {
        duplicated += Form("%s ", currPattern.Data());
        continue;
      }
    }
    
    TObjArray* trigCombo = new TObjArray();
    trigCombo->SetOwner();
    trigCombo->SetName(currPattern.Data());
    
    UInt_t uniqueID = kComboSimple;
    if ( ( hasAND && hasOR ) || hasNOT ) uniqueID = kComboFormula;
    else if ( hasAND ) uniqueID = kComboAND;
    else if ( hasOR ) uniqueID = kComboOR;
    
    trigCombo->SetUniqueID(uniqueID);
    
    TObjArray* arr = tn.Tokenize(":");
    
    TIter nextA(arr);
    TObjString* an = 0x0;
    while ( ( an = static_cast<TObjString*>(nextA()) ) )
    {
      Int_t listIdx = kTrigClass;
      if ( an->String().BeginsWith("0") ) listIdx = kL0Input;
      else if ( an->String().BeginsWith("1") ) listIdx = kL1Input;
      else if ( an->String().BeginsWith("2") ) listIdx = kL2Input;
      else if ( an->String().BeginsWith("k") ) listIdx = kPhysSelBit;
      
      TObjArray* currList = static_cast<TObjArray*>(trigCombo->FindObject(listName[listIdx].Data()));
      if ( ! currList ) {
        currList = new TObjArray();
        currList->SetOwner();
        currList->SetName(listName[listIdx].Data());
        currList->SetUniqueID(listIdx);
        trigCombo->Add(currList);
      }
      TObjString* currStr = new TObjString(an->String());
      
      TObject* toBeAdded = currStr;

      Bool_t isOk = kTRUE;

      if ( listIdx <= kL2Input ) {
        // that's an input
        TObject* trigInput = fTrigInputsMap->FindObject(an->String().Data());
        if ( trigInput ) currStr->SetUniqueID(trigInput->GetUniqueID());
        else {
          AliError(Form("Uknown input %s in formula %s", an->String().Data(), currPattern.Data()));
          isOk = kFALSE;
        }
      }
      else if ( listIdx == kPhysSelBit ) {
        // That's a physics selection bit
        TObject* physSelBit = fPhysSelBits->FindObject(an->String().Data());
        // FIXME: When AliBits will be in place, change this with:
        // toBeAdded = physSelBit
        if ( physSelBit ) currStr->SetUniqueID(physSelBit->GetUniqueID());
        else {
          AliError(Form("Uknown physSelBit %s in formula %s", an->String().Data(), currPattern.Data()));
          isOk = kFALSE;
        }
      }
      if ( ! isOk ) {
        delete trigCombo;
        trigCombo = 0x0;
        break;
      }
      currList->AddLast(toBeAdded);
    }
    delete arr;
    if ( trigCombo ) {
      fSelectedTrigCombination->AddLast(trigCombo);
      AliDebug(1,Form("Adding %s to trigger combination (type %u)",currPattern.Data(),trigCombo->GetUniqueID()));
    }
  }
  delete fullList;
  
  if ( ! duplicated.IsNull() )
    AliWarning(Form("Triggers %s already accounted in patterns",duplicated.Data()));
  if ( ! badSyntax.IsNull() )
    AliWarning(Form("%s : illegal expressions. Must be in the form:\n   pattern* => keep class if it contains pattern\n   !pattern* => reject class if it contains pattern\n   class&input = keep class if it satisfies the expression (exact matching required)",badSyntax.Data()));
}


//________________________________________________________________________
void AliMuonEventCuts::SetTrigClassLevels ( TString pattern )
{
  /// Set trigger cut level associated to the trigger class
  ///
  /// example:
  /// SetTrigClassLevels("MSL:Lpt,MSH:Hpt,MUL:LptLpt")
  ///
  /// For the trigger classes defined in SetTrigClassPatterns
  /// it check if they contains the keywords MSL or MSH
  /// Hence, in the analysis, the function
  /// TrackPtCutMatchTrigClass(track, "CPBIMSL") returns true if track match Lpt
  /// TrackPtCutMatchTrigClass(track, "CPBIMSH") returns true if track match Hpt
  /// TrackPtCutMatchTrigClass(track, "CMBAC") always returns true
  
  fSelectedTrigLevel->SetOwner();
  if ( fSelectedTrigLevel->GetEntries() > 0 ) fSelectedTrigLevel->Delete();
  
  pattern.ReplaceAll(" ","");
  
  TObjArray* fullList = pattern.Tokenize(",");
  UInt_t offset = 2;
  for ( Int_t ipat=0; ipat<fullList->GetEntries(); ++ipat ) {
    TString currPattern = fullList->At(ipat)->GetName();
    TObjArray* arr = currPattern.Tokenize(":");
    TObjString* trigClassPattern = new TObjString(arr->At(0)->GetName());
    TString selTrigLevel = arr->At(1)->GetName();
    selTrigLevel.ToUpper();
    UInt_t trigLevel = 0;
    if ( selTrigLevel.Contains("LPT") ) {
      trigLevel = 2;
      if ( selTrigLevel.Contains("LPTLPT") ) trigLevel += 2<<offset;
    }
    else if ( selTrigLevel.Contains("HPT") ) {
      trigLevel = 3;
      if ( selTrigLevel.Contains("HPTHPT") ) trigLevel += 3<<offset;
    }
    trigClassPattern->SetUniqueID(trigLevel);
    fSelectedTrigLevel->AddLast(trigClassPattern);
    delete arr;
  }
  
  delete fullList;
}

//________________________________________________________________________
UInt_t AliMuonEventCuts::GetTriggerInputBitMaskFromInputName(const char* inputName) const
{
  // Get trigger input bit from its name
  
  if (!fTrigInputsMap)
  {
    AliError("No Inputs Map available");
    return TMath::Limits<UInt_t>::Max();
  }
  
  TObjString* s = static_cast<TObjString*>(fTrigInputsMap->FindObject(inputName));
  if (!s)
  {
    AliError(Form("Did not find input %s",inputName));
    return TMath::Limits<UInt_t>::Max();    
  }
  return s->GetUniqueID();
}

//________________________________________________________________________
TArrayI AliMuonEventCuts::GetTrigClassPtCutLevel ( TString trigClassName ) const
{
  /// Get trigger class pt cut level for tracking/trigger matching
  ///
  /// CAVEAT: this functionality fully works with trigger class names,
  ///   but it can have a problem to extract the correct information in
  ///   combinations of trigger classes/inputs. In this case it provides:
  ///   - the highest pt level among classes/inputs combined through a logical AND "&"
  ///   - the lowest pt level among classes/inputs combined through a logical OR "|"
  ///   The first should be fine, but the second could not be the proper matching.
  TObject* obj = fAllSelectedTrigClasses->FindObject(trigClassName.Data());
  if ( ! obj ) {
    AliWarning(Form("Class %s not in the list!", trigClassName.Data()));
    return -1;
  }
  
  TArrayI ptCutLevel(2);
  ptCutLevel.Reset();
  ptCutLevel[0] = obj->GetUniqueID() & 0x3;
  ptCutLevel[1] = ( obj->GetUniqueID() >> 2 ) & 0x3;
  
  AliDebug(3,Form("Class %s ptCutLevel %i %i",trigClassName.Data(),ptCutLevel[0],ptCutLevel[1]));
  
  return ptCutLevel;
}

//________________________________________________________________________
void AliMuonEventCuts::SetPhysSelBits ()
{
  /// Set the correspondence between the physics selection
  /// bit name and its actual value

  if ( ! fPhysSelBits ) fPhysSelBits = new THashList();
  fPhysSelBits->SetOwner();

  TList* dmList = AliVEvent::Class()->GetListOfDataMembers();
  TDataMember* dm = 0x0;
  TIter next(dmList);
  TString typeName = "";
  while (( dm = static_cast<TDataMember*>(next()))) {
    typeName = dm->GetTypeName();
    if ( typeName != "AliVEvent::EOfflineTriggerTypes" ) continue;
    UInt_t bitVal = (UInt_t)gROOT->ProcessLineFast(Form("AliVEvent::%s",dm->GetName()));
    TObjString* physSelBit = new TObjString(dm->GetName());
    physSelBit->SetUniqueID(bitVal);
    fPhysSelBits->Add(physSelBit);
    AliDebug(3,Form("PhysSelBit %s 0x%x",dm->GetName(),bitVal));
  }
}


//________________________________________________________________________
void AliMuonEventCuts::SetTrigInputsMap ( TString trigInputsMap )
{
  /// Set trigger input mask
  /// The inputs must be in the form:
  /// input1:ID1,input2:ID2,...
  /// CAVEAT: the input ID is ID_aliceLogbook - 1
  /// since this is the ID in the OCDB
  
  fTrigInputsMap->SetOwner();
  if ( fTrigInputsMap->GetEntries() > 0 ) fTrigInputsMap->Delete();

  if ( trigInputsMap.IsNull() ) {
    AliWarning("Trigger input map not specified: using default");
    trigInputsMap = GetDefaultTrigInputsMap();
  }

  trigInputsMap.ReplaceAll(" ","");

  TObjArray* fullList = trigInputsMap.Tokenize(",");
  for ( Int_t ipat=0; ipat<fullList->GetEntries(); ++ipat ) {
    TString currPattern = fullList->At(ipat)->GetName();
    TObjArray* arr = currPattern.Tokenize(":");
    TObjString* trigInput = new TObjString(arr->At(0)->GetName());
    UInt_t trigID = (UInt_t)static_cast<TObjString*>(arr->At(1))->GetString().Atoi();
    trigInput->SetUniqueID(1<<trigID);
    fTrigInputsMap->Add(trigInput);
    delete arr;
  }
  delete fullList;
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
  
  AliDebug(2,Form("Fired classes: %s  Inputs 0x%x 0x%x 0x%x",firedTrigClasses.Data(),l0Inputs,l1Inputs,l2Inputs));
  
  if ( fSelectedTrigClassesInEvent ) fSelectedTrigClassesInEvent->Delete();
  else {
    fSelectedTrigClassesInEvent = new TObjArray(0);
    fSelectedTrigClassesInEvent->SetOwner();
  }
  

  TString firedTrigClassesAny = "ANY " + firedTrigClasses;
  
  if ( fSelectedTrigPattern->GetEntries() > 0 ) {
    TObjArray* firedTrigClassesList = firedTrigClassesAny.Tokenize(" ");
  
    for ( Int_t itrig=0; itrig<firedTrigClassesList->GetEntries(); ++itrig ) {
      TString trigName = ((TObjString*)firedTrigClassesList->At(itrig))->GetString();

      TObjString* foundTrig = static_cast<TObjString*>(fAllSelectedTrigClasses->FindObject(trigName.Data()));
      if ( ! foundTrig ) {
        if ( ! CheckTriggerClassPattern(trigName) ) continue;
      }
      
      AddToEventSelectedClass ( trigName, foundTrig );
    } // loop on trigger classes
  
    delete firedTrigClassesList;
  }
  
  for ( Int_t icomb=0; icomb<fSelectedTrigCombination->GetEntries(); icomb++ ) {
    TObjArray* currComb = static_cast<TObjArray*>(fSelectedTrigCombination->At(icomb));
    if ( CheckTriggerClassCombination(currComb, firedTrigClassesAny, l0Inputs, l1Inputs, l2Inputs,physicsSelection) ) {
      TObjString* foundTrig = static_cast<TObjString*>(fAllSelectedTrigClasses->FindObject(currComb->GetName()));
      AddToEventSelectedClass ( currComb->GetName(), foundTrig, currComb->GetUniqueID() );
    }
  }
}

//_____________________________________________________________________________
Bool_t
AliMuonEventCuts::CheckTriggerClassPattern ( const TString& toCheck ) const
{
  // Check if the "toCheck" class matches the user pattern
  
  for ( Int_t ipat=0; ipat<fRejectedTrigPattern->GetEntries(); ++ipat ) {
    if ( toCheck.Contains(TRegexp(fRejectedTrigPattern->At(ipat)->GetName(),kTRUE) ) ) return kFALSE;
  } // loop on reject pattern

  for ( Int_t ipat=0; ipat<fSelectedTrigPattern->GetEntries(); ++ipat ) {
    if ( toCheck.Contains(TRegexp(fSelectedTrigPattern->At(ipat)->GetName(),kTRUE) ) ) return kTRUE;
  } // loop on keep pattern
  
  return kFALSE;
}


//_____________________________________________________________________________
Bool_t
AliMuonEventCuts::CheckTriggerClassCombination ( const TObjArray* combo,
                                                 const TString& firedTriggerClasses,
                                                 UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs,
                                                 UInt_t physicsSelection ) const
{
  // Check if the "toCheck" class (or logical combination of classes and L0 inputs)
  // are within the "firedTriggerClasses"
  
  Bool_t ok(kFALSE);
  
  TString comp(combo->GetName());
  UInt_t trigInputs[3] = {l0Inputs, l1Inputs, l2Inputs};
  
  Bool_t exitLoop = kFALSE;
  
  TIter nextObj(combo);
  TObjArray* currList = 0x0;
  while ( ( currList = static_cast<TObjArray*>(nextObj()) ) ) {
    Int_t listIdx = currList->GetUniqueID();
    TIter nextA(currList);
    TObject* an = 0x0;
    while ( ( an = static_cast<TObjString*>(nextA()) ) )
    {
      if ( listIdx <= kL2Input ) {
        UInt_t bit = an->GetUniqueID();
        ok = ( (trigInputs[listIdx] & bit) == bit );
      }
      else if ( listIdx == kPhysSelBit ) {
        UInt_t bit = an->GetUniqueID();
        ok = ( physicsSelection & bit );
      }
      else {
        TPRegexp re(Form("(^|[ ])%s([ ]|$)",an->GetName()));
        ok = firedTriggerClasses.Contains(re);
      }
      if ( combo->GetUniqueID() == kComboFormula ) comp.ReplaceAll(an->GetName(),Form("%d",ok));
      else if ( ( combo->GetUniqueID() == kComboAND && ! ok ) || ( combo->GetUniqueID() == kComboOR && ok ) ) {
        exitLoop = kTRUE;
        break;
      }
    }
    if ( exitLoop ) break;
  }
  
  if ( combo->GetUniqueID() == kComboFormula ) {
    TFormula formula("TriggerClassFormulaCheck", comp.Data());
#if ROOT_VERSION_CODE < ROOT_VERSION(6,3,0)
    if ( formula.Compile() > 0 ) {
      AliError(Form("Could not evaluate formula %s",comp.Data()));
      ok = kFALSE;
    }
    else
#endif
      ok = formula.Eval(0);
  }
  
  AliDebug(2,Form("tname %s => %d comp=%s  inputs 0x%x 0x%x 0x%x",combo->GetName(),ok,comp.Data(),l0Inputs, l1Inputs, l2Inputs));
  
  return ok;
}

//_____________________________________________________________________________
void
AliMuonEventCuts::AddToEventSelectedClass ( const TString& toCheck, const TObjString* foundTrig, UInt_t comboType )
{
  /// Add current trigger to the selected class for the event
  
  // Compute the trigger pt cut level of the current class
  UInt_t trigLevel = 0;
  if ( foundTrig ) trigLevel = foundTrig->GetUniqueID();
  else {
    // The assigned trigger pt cut level is:
    // - the correct one if "toCheck" is a single trigger class
    // - the highest pt cut among the matching ones in case "toCheck" is a trigger class/input
    //   combined through (at least one) logical AND "&"
    // - the lowest pt cut among the macthing ones in case "toCheck" is a trigger class/input
    //   combined through (only) logical OR "|"
    // This may lead to errors in case of complex combinations of trigger/inputs

    // First eliminate trigger classes which are negated in combinations
    TString checkStr(toCheck);
    while ( checkStr.Contains("!") ) {
      Int_t startNot = checkStr.Index("!");
      Int_t endNot = startNot;
      Int_t npars = 0;
      for ( endNot = startNot; endNot<checkStr.Length(); endNot++ ) {
        if ( checkStr[endNot] == '(' ) npars++;
        else if ( checkStr[endNot] == ')' ) npars--;

        if ( npars == 0 ) {
          if ( checkStr[endNot] == '&' || checkStr[endNot] == '|' ) break;
        }
      }
      checkStr.Remove(startNot,endNot-startNot);
    }

    // Then check if they match the Lpt or Hpt
    Bool_t isFirst = kTRUE;
    for ( Int_t ipat=0; ipat<fSelectedTrigLevel->GetEntries(); ++ipat ) {
      if ( checkStr.Contains(fSelectedTrigLevel->At(ipat)->GetName() ) ) {
        UInt_t currLevel = fSelectedTrigLevel->At(ipat)->GetUniqueID();
        if ( comboType == kComboAND ) trigLevel = TMath::Max(trigLevel, currLevel);
        else if ( comboType == kComboOR || comboType == kComboFormula ) {
          if ( isFirst ) {
            trigLevel = currLevel;
            isFirst = kFALSE;
          }
          else trigLevel = TMath::Min(trigLevel, currLevel);
        }
        else {
          trigLevel = currLevel;
          break;
        }
      }
    } // loop on trig level patterns
  }
  
  TObjString* currTrig = new TObjString(toCheck);
  currTrig->SetUniqueID(trigLevel);
  fSelectedTrigClassesInEvent->AddLast(currTrig);
  
  if ( foundTrig ) return;
  TObjString* addTrig = new TObjString(toCheck);
  addTrig->SetUniqueID(trigLevel);
  fAllSelectedTrigClasses->Add(addTrig);
  TString trigLevelInfo = Form("trig level %i ", trigLevel & 0x3);
  trigLevelInfo += ( trigLevel > 3 ) ? "di-muon" : "single-muon";
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
  return evt->GetCentrality()->GetCentralityPercentile(GetCentralityEstimator());
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

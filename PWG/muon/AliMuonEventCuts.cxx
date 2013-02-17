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

#include "AliLog.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliCentrality.h"

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
  fDefaultTrigClassPatterns(""),
  fSelectedTrigPattern(0x0),
  fRejectedTrigPattern(0x0),
  fSelectedTrigLevel(0x0),
  fSelectedTrigCombination(0x0),
  fTrigInputsMap(0x0),
  fAllSelectedTrigClasses(0x0),
  fCentralityClasses(0x0),
  fEventTriggerMask(0),
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
  fDefaultTrigClassPatterns(""),
  fSelectedTrigPattern(new TObjArray()),
  fRejectedTrigPattern(new TObjArray()),
  fSelectedTrigLevel(new TObjArray()),
  fSelectedTrigCombination(new TObjArray()),
  fTrigInputsMap(new THashList()),
  fAllSelectedTrigClasses(new THashList()),
  fCentralityClasses(0x0),
  fEventTriggerMask(0),
  fSelectedTrigClassesInEvent(new TObjArray())
{
  /// Constructor
  SetDefaultParameters();
  SetDefaultFilterMask();
  SetDefaultTrigClassPatterns();
  SetTrigClassLevels();
  SetDefaultTrigInputsMap();
  SetCentralityClasses();
  fAllSelectedTrigClasses->SetOwner();
  fSelectedTrigClassesInEvent->SetOwner();
}

//________________________________________________________________________
AliMuonEventCuts::AliMuonEventCuts(const AliMuonEventCuts& obj) :
  AliAnalysisCuts(obj),
  fPhysicsSelectionMask(obj.fPhysicsSelectionMask),
  fVertexMinNContributors(obj.fVertexMinNContributors),
  fVertexVzMin(obj.fVertexVzMin),
  fVertexVzMax(obj.fVertexVzMax),
  fDefaultTrigClassPatterns(obj.fDefaultTrigClassPatterns),
  fSelectedTrigPattern(obj.fSelectedTrigPattern),
  fRejectedTrigPattern(obj.fRejectedTrigPattern),
  fSelectedTrigLevel(obj.fSelectedTrigLevel),
  fSelectedTrigCombination(obj.fSelectedTrigCombination),
  fTrigInputsMap(obj.fTrigInputsMap),
  fAllSelectedTrigClasses(obj.fAllSelectedTrigClasses),
  fCentralityClasses(obj.fCentralityClasses),
  fEventTriggerMask(obj.fEventTriggerMask),
  fSelectedTrigClassesInEvent(obj.fSelectedTrigClassesInEvent)
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
    fDefaultTrigClassPatterns = obj.fDefaultTrigClassPatterns;
    fSelectedTrigPattern = obj.fSelectedTrigPattern;
    fRejectedTrigPattern = obj.fRejectedTrigPattern;
    fSelectedTrigLevel = obj.fSelectedTrigLevel;
    fSelectedTrigCombination = obj.fSelectedTrigCombination;
    fTrigInputsMap = obj.fTrigInputsMap;
    fAllSelectedTrigClasses = obj.fAllSelectedTrigClasses;
    fCentralityClasses = obj.fCentralityClasses;
    fEventTriggerMask = obj.fEventTriggerMask;
    fSelectedTrigClassesInEvent = obj.fSelectedTrigClassesInEvent;
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
  delete fAllSelectedTrigClasses;
  delete fSelectedTrigClassesInEvent;
  delete fCentralityClasses;
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
  
  const AliInputEventHandler* inputHandler = static_cast<const AliInputEventHandler*> ( obj );
  
  if ( const_cast<AliInputEventHandler*>(inputHandler)->IsEventSelected() & fPhysicsSelectionMask ) selectionMask |= kPhysicsSelected;
  
  const AliVEvent* event = inputHandler->GetEvent();
  
  Double_t centrality = GetCentrality(event);
  if ( centrality >= fCentralityClasses->GetXmin() && centrality <= fCentralityClasses->GetXmax() ) selectionMask |= kSelectedCentrality;
  
  UpdateEvent(event);
  
  if ( fSelectedTrigClassesInEvent->GetEntries() > 0 ) selectionMask |= kSelectedTrig;
  
  AliVVertex* vertex = AliAnalysisMuonUtility::GetVertexSPD(event);
  if ( vertex->GetNContributors() >= GetVertexMinNContributors() && 
      vertex->GetZ() >= GetVertexVzMin() && vertex->GetZ() <= GetVertexVzMax() ) selectionMask |= kGoodVertex;
  
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
Bool_t AliMuonEventCuts::UpdateEvent ( const AliVEvent* event )
{
  /// Update the transient data member per event
  
  if ( fSelectedTrigClassesInEvent && ( fEventTriggerMask == event->GetTriggerMask() ) ) return kFALSE;
  
  BuildTriggerClasses(AliAnalysisMuonUtility::GetFiredTriggerClasses(event), AliAnalysisMuonUtility::GetL0TriggerInputs(event), AliAnalysisMuonUtility::GetL1TriggerInputs(event), AliAnalysisMuonUtility::GetL2TriggerInputs(event));
  
  fEventTriggerMask = event->GetTriggerMask();
  
  return kTRUE;
}

//________________________________________________________________________
void AliMuonEventCuts::SetDefaultTrigClassPatterns ()
{
  /// Set the default patterns
  /// (done in such a way to get all muon triggers)
  fDefaultTrigClassPatterns = "CM*,C0M*,CINT*,CPBI*,CCENT*,CV*,!*ABCE*,!*-ACE-*,!*-AC-*,!*-E-*,!*WU*,!*EGA*,!*EJE*,!*PHS*";
  SetTrigClassPatterns(fDefaultTrigClassPatterns);
}


//________________________________________________________________________
void AliMuonEventCuts::SetTrigClassPatterns ( const TString trigPattern )
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
  /// 2) specify a combination of triggers
  /// combined through a logical AND "&" or a logical OR "|" (wildcard * NOT accepted)
  /// It is also possible to ask for a trigger class containing a specific trigger input:
  /// e.g. CMSL7-B-NOPF-MUON&0MSH,CMSL7-B-NOPF-MUON,CMSL7-B-NOPF-MUON|CMSL8-B-NOPF-MUON
  /// will give the events with:
  /// - the single low trigger class fired and containing a single high trigger input
  /// - the single low trigger class fired
  /// - the single low trigger class 7 or 8 fired
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
  
  TString badSyntax = "", duplicated = "";
  TString listName[4] = {"L0","L1","L2","trigClass"};
  
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
    Bool_t isMatchPattern = currPattern.Contains("*");
    Bool_t isRejectPattern = kFALSE;
    if ( currPattern.Contains("!") ) {
      currPattern.ReplaceAll("!","");
      isRejectPattern = kTRUE;
    }
    if ( isCombination && ( isMatchPattern || isRejectPattern ) ) {
      badSyntax += Form(" %s;", currPattern.Data());
      continue;
    }
    if ( isRejectPattern ) {
      fRejectedTrigPattern->AddLast(new TObjString(currPattern));
      AliDebug(2,Form("Adding %s to reject pattern",currPattern.Data()));
    }
    else if ( isMatchPattern ) {
      fSelectedTrigPattern->AddLast(new TObjString(currPattern));
      AliDebug(2,Form("Adding %s to match pattern",currPattern.Data()));
    }
    else combinationList.Add(objString);
  }
  
  // Then check for combinations
  TIter nextCombo(&combinationList);
  while ( ( objString = static_cast<TObjString*>(nextCombo()) ) ) {
    TString currPattern = objString->String();
    
    TString tn (currPattern);
    Bool_t isCombination = kFALSE;
    Bool_t requiresFromula = kFALSE;
    if ( tn.Contains("&") ) {
      tn.ReplaceAll("&",":");
      isCombination = kTRUE;
    }
    if ( tn.Contains("|") ) {
      tn.ReplaceAll("|",":");
      isCombination = kTRUE;
      requiresFromula = kTRUE;
    }
    if ( tn.Contains("(") || tn.Contains(")") ) {
      tn.ReplaceAll("(","");
      tn.ReplaceAll(")","");
    }
    
    if ( ! isCombination ) {
      if ( CheckTriggerClassPattern(currPattern) ) {
        duplicated += Form("%s ", currPattern.Data());
        continue;
      }
    }
    
    TObjArray* trigCombo = new TObjArray();
    trigCombo->SetOwner();
    trigCombo->SetName(currPattern.Data());
    
    
    if ( requiresFromula ) trigCombo->SetUniqueID(1);
    
    TObjArray* arr = tn.Tokenize(":");
    
    TIter nextA(arr);
    TObjString* an = 0x0;
    while ( ( an = static_cast<TObjString*>(nextA()) ) )
    {
      Int_t listIdx = 3;
      if ( an->String().BeginsWith("0") ) listIdx = 0;
      else if ( an->String().BeginsWith("1") ) listIdx = 1;
      else if ( an->String().BeginsWith("2") ) listIdx = 2;
      
      TObjArray* currList = static_cast<TObjArray*>(trigCombo->FindObject(listName[listIdx].Data()));
      if ( ! currList ) {
        currList = new TObjArray();
        currList->SetOwner();
        currList->SetName(listName[listIdx].Data());
        currList->SetUniqueID(listIdx);
        trigCombo->AddAt(currList,listIdx);
      }
      TObjString* currStr = new TObjString(an->String());
      
      if ( listIdx < 3 ) {
        // that's an input
        TObject* trigInput = fTrigInputsMap->FindObject(an->String().Data());
        if ( trigInput ) currStr->SetUniqueID(trigInput->GetUniqueID());
        else {
          AliError(Form("Uknown input %s in formula %s", an->String().Data(), currPattern.Data()));
          delete trigCombo;
          trigCombo = 0x0;
          break;
        }
      }
      currList->AddLast(currStr);
    }
    delete arr;
    if ( trigCombo) {
      fSelectedTrigCombination->AddLast(trigCombo);
      AliDebug(2,Form("Adding %s to trigger combination",currPattern.Data()));
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
TArrayI AliMuonEventCuts::GetTrigClassPtCutLevel ( const TString trigClassName ) const
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
  
  AliDebug(1,Form("Class %s ptCutLevel %i %i",trigClassName.Data(),ptCutLevel[0],ptCutLevel[1]));
  
  return ptCutLevel;
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
void AliMuonEventCuts::SetDefaultTrigInputsMap ()
{
  /// Set default trigger input mask
  
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
  
  SetTrigInputsMap(trigInputsMap);
}

//________________________________________________________________________
TObjArray* AliMuonEventCuts::GetSelectedTrigClassesInEvent( const AliVEvent* event )
{
  /// Return the selected trigger classes in the current event
  UpdateEvent(event);
  return fSelectedTrigClassesInEvent;
}


//________________________________________________________________________
void AliMuonEventCuts::BuildTriggerClasses ( const TString firedTrigClasses,
                                             UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs )
{
  //
  /// Return the list of trigger classes to be considered
  /// for current event. Update the global list if necessary
  //
  
  AliDebug(2,Form("Fired classes: %s  Inputs 0x%x 0x%x 0x%x",firedTrigClasses.Data(),l0Inputs,l1Inputs,l2Inputs));
  
  if ( fSelectedTrigClassesInEvent) fSelectedTrigClassesInEvent->Delete();
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
    if ( CheckTriggerClassCombination(currComb, firedTrigClassesAny, l0Inputs, l1Inputs, l2Inputs) ) {
      TObjString* foundTrig = static_cast<TObjString*>(fAllSelectedTrigClasses->FindObject(currComb->GetName()));
      AddToEventSelectedClass ( currComb->GetName(), foundTrig );
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
                                                 UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs ) const
{
  // Check if the "toCheck" class (or logical combination of classes and L0 inputs)
  // are within the "firedTriggerClasses"
  
  Bool_t ok(kFALSE);
  
  TString comp(combo->GetName());
  UInt_t trigInputs[3] = {l0Inputs, l1Inputs, l2Inputs};
  Bool_t requiresFromula = ( combo->GetUniqueID() == 1 );
  
  Bool_t exitLoop = kFALSE;
  
  TIter nextObj(combo);
  TObjArray* currList = 0x0;
  while ( ( currList = static_cast<TObjArray*>(nextObj()) ) ) {
    Int_t listIdx = currList->GetUniqueID();
    TIter nextA(currList);
    TObjString* an = 0x0;
    while ( ( an = static_cast<TObjString*>(nextA()) ) )
    {
      if ( listIdx < 3 ) {
        UInt_t bit = an->GetUniqueID();
        Bool_t matchInput = ( (trigInputs[listIdx] & bit) == bit );
        if ( requiresFromula ) comp.ReplaceAll(an->String().Data(),( matchInput ) ? "1" : "0");
        else ok = matchInput;
      }
      else {
        TPRegexp re(Form("(^|[ ])%s([ ]|$)",an->String().Data()));
        Bool_t matchTrig = firedTriggerClasses.Contains(re);
        if ( requiresFromula ) comp.ReplaceAll(an->String().Data(),Form("%d",matchTrig));
        else ok = matchTrig;
      }
      if ( ! requiresFromula && ! ok ) {
        exitLoop = kTRUE;
        break;
      }
    }
    if ( exitLoop ) break;
  }
  
  if ( requiresFromula ) {
    TFormula formula("TriggerClassFormulaCheck", comp.Data());
    if ( formula.Compile() > 0 ) AliError(Form("Could not evaluate formula %s",comp.Data()));
    else ok = formula.Eval(0);
  }
  
  AliDebug(2,Form("tname %s => %d comp=%s  inputs 0x%x 0x%x 0x%x",combo->GetName(),ok,comp.Data(),l0Inputs, l1Inputs, l2Inputs));
  
  return ok;
}

//_____________________________________________________________________________
void
AliMuonEventCuts::AddToEventSelectedClass ( const TString& toCheck, const TObjString* foundTrig )
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
    for ( Int_t ipat=0; ipat<fSelectedTrigLevel->GetEntries(); ++ipat ) {
      if ( toCheck.Contains(fSelectedTrigLevel->At(ipat)->GetName() ) ) {
        UInt_t currLevel = fSelectedTrigLevel->At(ipat)->GetUniqueID();
        if ( toCheck.Contains("&") ) trigLevel = TMath::Max(trigLevel, currLevel);
        else if ( toCheck.Contains("|") ) trigLevel = TMath::Min(trigLevel, currLevel);
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
void AliMuonEventCuts::SetCentralityEstimator ( const TString centralityEstimator )
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
    if ( filterMask & kSelectedCentrality ) printf(  "%g < centrality < %g\n", fCentralityClasses->GetXmin(), fCentralityClasses->GetXmax() );
    if ( filterMask & kSelectedTrig )    printf("  Has selected trigger classes\n");
    if ( filterMask & kGoodVertex )      printf("  SPD vertex with %i contributors && %g < Vz < %g\n", GetVertexMinNContributors(), GetVertexVzMin(), GetVertexVzMax());
    printf(" ******************** \n");
  }
}

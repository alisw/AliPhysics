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


#include "AliMuonTriggerCombo.h"

#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TList.h"
#include "THashList.h"
#include "TFormula.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TDataMember.h"
#include "TROOT.h"

#include "AliLog.h"
#include "AliVEvent.h"

/// \cond CLASSIMP
ClassImp(AliMuonTriggerCombo) // Class implementation in ROOT context
/// \endcond

//________________________________________________________________________
AliMuonTriggerCombo::AliMuonTriggerCombo () :
TNamed(),
fL0Inputs(0x0),
fL1Inputs(0x0),
fL2Inputs(0x0),
fPhysSelBits(0x0),
fTriggerClasses(0x0),
fComboType(0),
fTrigPtMatchLevel(0),
fIsDimuTrig(0)
{
  /// Defult ctr.
}

//________________________________________________________________________
AliMuonTriggerCombo::AliMuonTriggerCombo ( const char* name, const char* trigInputsString, const char* trigPtMatchLevel ) :
TNamed(name,trigPtMatchLevel),
fL0Inputs(0x0),
fL1Inputs(0x0),
fL2Inputs(0x0),
fPhysSelBits(0x0),
fTriggerClasses(0x0),
fComboType(0),
fTrigPtMatchLevel(0),
fIsDimuTrig(0)
{
  /// Ctr

  Init(name,trigInputsString,trigPtMatchLevel);
}


//________________________________________________________________________
Bool_t AliMuonTriggerCombo::Init ( const char* name, const char* trigInputsString, const char* trigPtMatchLevel )
{
  /// Initialise trigger pattern or trigger combination
  ///
  /// 1) A trigger pattern is a perl regular expression.
  /// The fired trigger will be kept only if it matches the regex
  /// Examples:
  /// - CINT[78]-B-NOPF-MUFAST* : matches CINT7-B-NOPF-MUFAST and CINT8-B-NOPF-MUFAST
  ///
  /// 2) specify a combination of trigger classes, trigger inputs or physics selection bits
  /// combined through a logical AND "&" or a logical OR "|" (regex * NOT accepted)
  /// Trigger input sbegins with 0, 1 or 2
  /// Physics selction bits begins with k, and correspond to the enum of AliVEvent
  /// Examples:
  /// - CMSL7-B-NOPF-MUON&0MSH : matches the single low trigger class fired and containing a single high trigger input
  /// - CMSL7-B-NOPF-MUON : matches the single low trigger class fired
  /// - CMSL7-B-NOPF-MUON|CMSL8-B-NOPF-MUON : maches the single low trigger class 7 or 8 fired
  /// - kMUU7|kMUL7 : matches events firing either AliVEvent::kMUU7 or AliVEvent::kMUL7
  ///
  /// The trigger inputs must be specified if one wants to use them.
  /// The inputs must be in the form:
  /// input1:ID1,input2:ID2,...
  /// CAVEAT: the input ID is ID_aliceLogbook - 1
  /// since this is the ID in the OCDB
  ///
  /// It is also possible to provide a trigger pt level for tracking/trigger matching
  /// associated to the trigger combination.
  /// This information can be used by AliMuonTrackCuts::TrackPtCutMatchTrigClass(AliVTrack* track, AliMuonEventCuts::GetTrigClassPtCutLevel);
  /// and by:
  /// AliMuonPairCuts::TrackPtCutMatchTrigClass(AliVTrack* track1, AliVTrack* track2, AliMuonEventCuts::GetTrigClassPtCutLevel);
  /// Valid levels are:
  /// - Lpt : single muon Lpt
  /// - Lpt2 : di-muon Lpt
  /// - Hpt : single muon Hpt
  /// - Hpt2 : di-muon Hpt
  /// Example:
  /// - CMSL7*:Lpt : matches a single Lpt trigger classes and assign a Lpt level for matching
  /// - kMUU7|kMUL7:Lpt2 : matches events fired by either AliVEvent::kMUU7 or AliVEvent::kMUL7 and assigns a di0muon Lpt trigger level

  // Set trigger level if any
  TString matchPt(trigPtMatchLevel);
  matchPt.ToUpper();
  if ( matchPt.Contains("LPT") ) fTrigPtMatchLevel = 2;
  else if ( matchPt.Contains("HPT") ) fTrigPtMatchLevel = 3;

  if ( matchPt.Contains("2") || matchPt.Contains("DI") ) fIsDimuTrig = kTRUE;

  TString tn (name);
  Bool_t hasAND = tn.Contains("&");
  Bool_t hasOR = tn.Contains("|");
  Bool_t hasNOT = tn.Contains("!");

  Bool_t isRegex = kFALSE;
  if ( tn.Contains("[") || tn.Contains(".") || tn.Contains("+") || tn.Contains("?") || tn.Contains("*") || tn.Contains("^") || tn.Contains("$") ) {
    // This is a regular expression
    if ( hasAND ||
        ( hasOR && ! tn.Contains("(") ) )
    {
       AliError(Form("%s : illegal expressions. Must be in the form:\n   regexp => keep class if it matches the regexp\n   class&input = keep class if it satisfies the expression (exact matching required)",GetName()));
       fComboType = kBadPattern;
       return kFALSE;
     }
     fComboType = kRegex;
  }
  else {
    tn.ReplaceAll("&",";");
    tn.ReplaceAll("|",";");
    tn.ReplaceAll("!","");
    tn.ReplaceAll("(","");
    tn.ReplaceAll(")","");

    if ( ( hasAND && hasOR ) || hasNOT ) fComboType = kComboFormula;
    else if ( hasAND ) fComboType = kComboAND;
    else if ( hasOR ) fComboType = kComboOR;
    else fComboType = kComboSimple;
  }

  THashList* trigInputsMap = GetTrigInputsMap(trigInputsString);
  THashList* physSelBits = GetPhysSelBits();

  TObjArray* arr = tn.Tokenize(";");

  TIter nextA(arr);
  TObjString* an = 0x0;
  TString currString = "";
  while ( ( an = static_cast<TObjString*>(nextA()) ) )
  {
    currString = an->GetString();
    if ( currString.BeginsWith("0") || currString.BeginsWith("1") || currString.BeginsWith("2") ) {
      TObject* foundInput = ( trigInputsMap ) ? trigInputsMap->FindObject(currString.Data()) : 0x0;
      if ( ! foundInput ) {
        AliError(Form("Uknown input %s in formula %s", currString.Data(), name));
        fComboType = kBadPattern;
        break;
      }
      TObjString* trigInput = new TObjString(currString);
      trigInput->SetUniqueID(foundInput->GetUniqueID());
      TString strLevel = an->String()[0];
      GetInputs(strLevel.Atoi())->Add(trigInput);
    }
    else if ( an->String().BeginsWith("k") ) {
      // That's a physics selection bit
      TObject* foundPhysSelBit = physSelBits->FindObject(currString.Data());
      // FIXME: When AliBits will be in place, change this with:
      // toBeAdded = physSelBit
      if ( ! foundPhysSelBit ) {
        AliError(Form("Uknown physSelBit %s in formula %s", currString.Data(), name));
        fComboType = kBadPattern;
        break;
      }
      TObjString* physSelBit = new TObjString(currString);
      physSelBit->SetUniqueID(foundPhysSelBit->GetUniqueID());
      if ( ! fPhysSelBits ) {
        fPhysSelBits = new TObjArray();
        fPhysSelBits->SetOwner();
      }
      fPhysSelBits->Add(physSelBit);
    }
    else {
      if ( ! fTriggerClasses ) {
        fTriggerClasses = new TObjArray();
        fTriggerClasses->SetOwner();
      }
      fTriggerClasses->Add(new TObjString(currString));
    }
  }

  delete trigInputsMap;
  delete physSelBits;
  delete arr;

  return ( fComboType != kBadPattern );
}

//________________________________________________________________________
AliMuonTriggerCombo::~AliMuonTriggerCombo ()
{
  /// Destructor
  delete fL0Inputs;
  delete fL1Inputs;
  delete fL2Inputs;
  delete fPhysSelBits;
  delete fTriggerClasses;
}

//________________________________________________________________________
TObjArray* AliMuonTriggerCombo::GetInputs ( Int_t level )
{
  /// Get trigger inputs

  if ( level == 0 ) {
    if ( ! fL0Inputs ) {
      fL0Inputs = new TObjArray();
      fL0Inputs->SetOwner();
    }
    return fL0Inputs;
  }
  else if ( level == 1 ) {
    if ( ! fL1Inputs ) {
      fL1Inputs = new TObjArray();
      fL1Inputs->SetOwner();
    }
    return fL1Inputs;
  }
  else if ( level == 2 ) {
    if ( ! fL2Inputs ) {
      fL2Inputs = new TObjArray();
      fL2Inputs->SetOwner();
    }
    return fL2Inputs;
  }
  return 0x0;
}


//________________________________________________________________________
THashList* AliMuonTriggerCombo::GetTrigInputsMap ( const char *trigInputsString ) const
{
  /// Get trigger inputs ID

  TString sTrigInputs(trigInputsString);
  sTrigInputs.ReplaceAll(" ","");

  THashList* trigInputsMap = new THashList();
  trigInputsMap->SetOwner();

  TObjArray* fullList = sTrigInputs.Tokenize(",");
  for ( Int_t ipat=0; ipat<fullList->GetEntries(); ipat++ ) {
    TString currPattern = fullList->At(ipat)->GetName();
    TObjArray* arr = currPattern.Tokenize(":");
    TObjString* trigInput = new TObjString(arr->At(0)->GetName());
    UInt_t trigID = (UInt_t)static_cast<TObjString*>(arr->At(1))->GetString().Atoi();
    trigInput->SetUniqueID(1<<trigID);
    trigInputsMap->Add(trigInput);
    delete arr;
  }
  delete fullList;

  return trigInputsMap;
}

//________________________________________________________________________
THashList* AliMuonTriggerCombo::GetPhysSelBits () const
{
  /// Set the correspondence between the physics selection
  /// bit name and its actual value

  THashList* physSelBits = new THashList();
  physSelBits->SetOwner();

  TList* dmList = AliVEvent::Class()->GetListOfDataMembers(); // Not owner
  TDataMember* dm = 0x0;
  TIter next(dmList);
  TString typeName = "";
  while (( dm = static_cast<TDataMember*>(next()))) {
    typeName = dm->GetTypeName();
    if ( typeName != "AliVEvent::EOfflineTriggerTypes" ) continue;
    UInt_t bitVal = (UInt_t)gROOT->ProcessLineFast(Form("AliVEvent::%s",dm->GetName()));
    TObjString* physSelBit = new TObjString(dm->GetName());
    physSelBit->SetUniqueID(bitVal);
    physSelBits->Add(physSelBit);
    AliDebug(3,Form("PhysSelBit %s 0x%x",dm->GetName(),bitVal));
  }
  return physSelBits;
}

//________________________________________________________________________
Bool_t AliMuonTriggerCombo::CheckElement ( TString& formula, const TObjString* element, Bool_t ok ) const
{
  /// Check if element is ok
  Bool_t stopCheck = kFALSE;
  if ( fComboType == kComboFormula ) formula.ReplaceAll(element->GetName(),Form("%d",ok));
  else if ( ( fComboType == kComboAND && ! ok ) || ( fComboType == kComboOR && ok ) ) stopCheck = kTRUE;
  return stopCheck;
}

//________________________________________________________________________
Bool_t AliMuonTriggerCombo::MatchEvent ( const TString& firedTriggerClasses, UInt_t l0Inputs, UInt_t l1Inputs, UInt_t l2Inputs, UInt_t physicsSelection )
{
  /// Check if current event satisfies the trigger combination requirements
  Bool_t ok(kFALSE);

  TString debugString = Form("Classes: %s  inputs: 0x%x 0x%x 0x%x  PhsySel: 0x%x  Match %s =>",firedTriggerClasses.Data(),l0Inputs,l1Inputs,l2Inputs,physicsSelection,GetName());

  TObjString* an = 0x0;
  if ( fComboType == kRegex ) {
    an = static_cast<TObjString*>(fTriggerClasses->At(0));
    TPRegexp re(an->GetName());
    ok = firedTriggerClasses.Contains(re);
    AliDebug(2,Form("%s %i",debugString.Data(),ok));
    return  ok;
  }

  // Obsolete
  // if ( fComboType == kRejectPattern ) {
  //   an = static_cast<TObjString*>(fTriggerClasses->At(0));
  //   ok = ( ! firedTriggerClasses.Contains(TRegexp(an->GetName(),kTRUE)) );
  //   AliDebug(2,Form("%s %i",debugString.Data(),ok));
  //   return ok;
  // }

  TString comp(GetName());

  TIter nextTrig(fTriggerClasses);
  while ( ( an = static_cast<TObjString*>(nextTrig()) ) ) {
    TPRegexp re(Form("(^|[ ])%s([ ]|$)",an->GetName()));
    ok = firedTriggerClasses.Contains(re);
    if ( CheckElement(comp, an, ok) ) {
      AliDebug(2,Form("%s %i",debugString.Data(),ok));
      return ok;
    }
  }

  TIter nextPhysSel(fPhysSelBits);
  while ( ( an = static_cast<TObjString*>(nextPhysSel()) ) ) {
    UInt_t bit = an->GetUniqueID();
    ok = ( physicsSelection & bit );
    if ( CheckElement(comp, an, ok) ) {
      AliDebug(2,Form("%s %i",debugString.Data(),ok));
      return ok;
    }
  }

  UInt_t trigInputs[3] = {l0Inputs, l1Inputs, l2Inputs};
  for ( Int_t ilevel=0; ilevel<3; ilevel++ ) {
    TIter nextInput(GetInputs(ilevel));
    while ( ( an = static_cast<TObjString*>(nextInput()) ) ) {
      UInt_t bit = an->GetUniqueID();
      ok = ( (trigInputs[ilevel] & bit) == bit );
      if ( CheckElement(comp, an, ok) ) {
        AliDebug(2,Form("%s %i",debugString.Data(),ok));
        return ok;
      }
    }
  }

  if ( fComboType == kComboFormula ) {
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

  AliDebug(2,Form("%s %i",debugString.Data(),ok));

  return ok;
}

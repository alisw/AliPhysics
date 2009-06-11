/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////
//
//  Finds events according to specified criteria
//
//  origin: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
//////////////////////////////////////////////////////////////////////////

#include "AliEveEventSelector.h"
#include "AliEveEventManager.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDRun.h"
#include <TEntryList.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TRegexp.h>
#include <TROOT.h>

ClassImp(AliEveEventSelector)

//_____________________________________________________________________________
AliEveEventSelector::AliEveEventSelector(AliEveEventManager* evman):
  fPEventManager(evman),
  fWrapAround(kTRUE),
  fMaxEventId(fPEventManager->GetMaxEventId()),
  fSelectOnString(kFALSE),
  fString(""),
  fPEntryList(NULL),
  fEntryListId(0),
  fSelectOnTriggerType(kFALSE),
  fTriggerType(""),
  fSelectOnTriggerString(kFALSE),
  fTriggerSelectionString(""),
  fTriggerMaskPatternString("trgmask"),
  fSelectOnMultiplicity(kFALSE),
  fMultiplicityLow(0),
  fMultiplicityHigh(0)
{
  //ctor
}

//_____________________________________________________________________________
void AliEveEventSelector::SetSelectionString( const TString& str )
{
  //selection string setter

  TTree* pESDTree = fPEventManager->GetESDTree();
  if (!pESDTree || fString==str) return;//don't waist time
  fString = str;
  if (str == "" ) return;//on reset don't recalculate
  fMaxEventId = fPEventManager->GetMaxEventId();
  pESDTree->Draw( ">>listofentries", fString, "entrylist");
  fPEntryList = dynamic_cast<TEntryList*>(gDirectory->Get("listofentries"));
}

//_____________________________________________________________________________
void AliEveEventSelector::SetSelectionString( const char* str )
{
  TString ts = str;
  SetSelectionString(ts);
}

//_____________________________________________________________________________
void AliEveEventSelector::SetTriggerType( const TString& type )
{
  fTriggerType = type;
}

//_____________________________________________________________________________
void AliEveEventSelector::SetTriggerType( const char* type)
{
  TString ts = type;
  SetTriggerType(ts);
}

//_____________________________________________________________________________
void AliEveEventSelector::UpdateEntryList()
{
  //update the entrylist from file if file changed
  TTree* pESDTree = fPEventManager->GetESDTree();
  if (!pESDTree) return;
  
  Long64_t treesize = fPEventManager->GetMaxEventId();
  if (treesize<=fMaxEventId) return; //nothing changed, do nothing
  pESDTree->Draw(">>+fPEntryList", fString, "entrylist",
      fMaxEventId, treesize-fMaxEventId );
  fMaxEventId = treesize;
}

//_____________________________________________________________________________
void AliEveEventSelector::Update()
{
  //refresh stuff

  UpdateEntryList();
}

//_____________________________________________________________________________
Bool_t AliEveEventSelector::FindNext( Int_t& output )
{
  //does the magick of selecting the next event

  static const TEveException kEH("AliEveEventSelector::GetNext ");
  
  TTree* pESDTree = fPEventManager->GetESDTree();
  if (!pESDTree) return kFALSE;
  AliESDEvent* pESDEvent = fPEventManager->GetESD();
  Int_t eventId = fPEventManager->GetEventId();
  fMaxEventId = fPEventManager->GetMaxEventId();

  if (!fSelectOnString)
  {
    // pure non-string
    for (Int_t i = eventId+1; i<fMaxEventId+1; i++)
    {
      if (pESDTree->GetEntry(i) <= 0)
        throw (kEH + "failed getting required event from ESD.");
      if (CheckOtherSelection(pESDEvent))
      {
        output = i;
        return kTRUE;
      }
    }
    if (!fWrapAround) return kFALSE;
    for (Int_t i = 0; i<eventId+1; i++)
    {
      if (pESDTree->GetEntry(i) <= 0)
        throw (kEH + "failed getting required event from ESD.");
      if (CheckOtherSelection(pESDEvent))
      {
        output = i;
        return kTRUE;
      }
    }
    return kFALSE;
  }
  else
  {
    //select using the entrylist
    for (Long64_t i=fEntryListId+1; i<fPEntryList->GetN(); i++ )
    {
      Long64_t entry = fPEntryList->GetEntry(i);
      if (pESDTree->GetEntry(entry) <= 0)
        throw (kEH + "failed getting required event from ESD.");
      if (CheckOtherSelection(pESDEvent))
      {
        output = entry;
        fEntryListId = i;
        return kTRUE;
      }
    }
    if (!fWrapAround) return kFALSE;
    for (Long64_t i=0; i<fEntryListId+1; i++ )
    {
      Long64_t entry = fPEntryList->GetEntry(i);
      if (pESDTree->GetEntry(entry) <= 0)
        throw (kEH + "failed getting required event from ESD.");
      if (CheckOtherSelection(pESDEvent))
      {
        output = entry;
        fEntryListId=i;
        return kTRUE;
      }
    }
    return kFALSE;
  }
  return kFALSE;

}

//_____________________________________________________________________________
Bool_t AliEveEventSelector::FindPrev( Int_t& output )
{
  //does the magick of selecting the previous event

  static const TEveException kEH("AliEveEventSelector::GetNext ");
  
  TTree* pESDTree = fPEventManager->GetESDTree();
  if (!pESDTree) return kFALSE;
  AliESDEvent* pESDEvent = fPEventManager->GetESD();
  Int_t eventId = fPEventManager->GetEventId();
  fMaxEventId = fPEventManager->GetMaxEventId();

  if (!fSelectOnString)
  {
    // pure non-string
    for (Int_t i = eventId-1; i>-1; i--)
    {
      if (pESDTree->GetEntry(i) <= 0)
        throw (kEH + "failed getting required event from ESD.");
      if (CheckOtherSelection(pESDEvent))
      {
        output = i;
        return kTRUE;
      }
    }
    if (!fWrapAround) return kFALSE;
    for (Int_t i = fMaxEventId; i>eventId-1; i--)
    {
      if (pESDTree->GetEntry(i) <= 0)
        throw (kEH + "failed getting required event from ESD.");
      if (CheckOtherSelection(pESDEvent))
      {
        output = i;
        return kTRUE;
      }
    }
    return kFALSE;
  }
  else
  {
    //select using the entrylist
    for (Long64_t i=fEntryListId-1; i>-1; i--)
    {
      Long64_t entry = fPEntryList->GetEntry(i);
      if (pESDTree->GetEntry(entry) <= 0)
        throw (kEH + "failed getting required event from ESD.");
      if (CheckOtherSelection(pESDEvent))
      {
        output = entry;
        fEntryListId = i;
        return kTRUE;
      }
    }
    if (!fWrapAround) return kFALSE;
    for (Long64_t i=fPEntryList->GetN()-1; i>fEntryListId-1; i--)
    {
      Long64_t entry = fPEntryList->GetEntry(i);
      if (pESDTree->GetEntry(entry) <= 0)
        throw (kEH + "failed getting required event from ESD.");
      if (CheckOtherSelection(pESDEvent))
      {
        output = entry;
        fEntryListId=i;
        return kTRUE;
      }
    }
    return kFALSE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliEveEventSelector::CheckOtherSelection( AliESDEvent* pESDEvent )
{
  //checks the event for any other hardcoded selection criteria
  Bool_t ret=kTRUE;

  //trigger stuff
  if (fSelectOnTriggerType)
  {
    TString firedtrclasses = pESDEvent->GetFiredTriggerClasses();
    printf(firedtrclasses.Data());
    printf("\n");
    printf(fTriggerType.Data());
    printf("\n");
    if (!(firedtrclasses.Contains(fTriggerType))) return kFALSE;
    //if (!pESDEvent->IsTriggerClassFired(fTriggerType.Data())) return kFALSE;
  }

  if (fSelectOnMultiplicity)
  {
    Int_t mult = pESDEvent->GetNumberOfTracks();
    Int_t mhigh = (fMultiplicityHigh==0)?100000000:fMultiplicityHigh;
    if (mult<fMultiplicityLow || mult>mhigh) return kFALSE;
  }

  if (fSelectOnTriggerString)
  {
    ULong64_t triggermask = pESDEvent->GetTriggerMask();
    TString triggermaskstr;
    triggermaskstr += triggermask;
    TString selstr(fTriggerSelectionString); //make copy
    selstr.ReplaceAll(fTriggerMaskPatternString,triggermaskstr);
    Int_t returncode;
    Bool_t result = static_cast<Bool_t>(gROOT->ProcessLine(selstr,&returncode));
    //if (!returncode) return kFALSE;
    if (!result) return kFALSE;
  }

  return ret;
}

//_____________________________________________________________________________
void AliEveEventSelector::SetTriggerSelectionString( TString str )
{
  const AliESDRun* run = fPEventManager->GetESD()->GetESDRun();
  for (Int_t i=0; i<run->kNTriggerClasses; i++)
  {
    TString name(run->GetTriggerClass(i));
    if (name.IsNull()) continue;
    TString valuestr("(");
    valuestr += fTriggerMaskPatternString;
    valuestr += "&";
    valuestr += static_cast<ULong64_t>(1<<i);
    valuestr += ")";
    str.ReplaceAll(name,valuestr);
  }//for i
  fTriggerSelectionString = str;
}


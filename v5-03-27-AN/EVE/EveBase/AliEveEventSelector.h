/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//
//  AliEveEventSelector class
//  selects events according to given criteria
//
//  origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
////////////////////////////////////////////////////////////////////////////

#ifndef AliEveEventSelector_H
#define AliEveEventSelector_H

#include "AliEveEventManager.h"
#include <TObject.h>
#include <TString.h>

class TTree;
class TList;
class AliESDEvent;
class TEntryList;
class AliEveEventManager;
class TString;
class TObjArray;

class AliEveEventSelector : public TObject
{
public:
  AliEveEventSelector(AliEveEventManager* evman);
  void SetEventManager(AliEveEventManager* evman) { fPEventManager = evman; }
  AliEveEventManager* GetEventManager() { return fPEventManager; }
  TTree* GetESDTree() const { return fPEventManager->GetESDTree(); }

  Bool_t FindNext(Int_t& i);
  Bool_t FindPrev(Int_t& i);
  
  void SetSelectionString( const TString& str );
  void SetSelectionString( const char* str );
  TString GetSelectionString() {return fString;}
  Bool_t GetSelectOnString() {return fSelectOnString;}
  void SetSelectOnString( const Bool_t b = kTRUE ) {fSelectOnString = b;}
  void SetWrapAround( const Bool_t b = kTRUE ) {fWrapAround = b;}
  Bool_t GetWrapAround() const { return fWrapAround; }
  void Update();

  void SetSelectOnTriggerType( const Bool_t b = kTRUE ) {fSelectOnTriggerType = b;}
  void SetTriggerType( const TString& type );
  void SetTriggerType( const char* type );
  TString GetTriggerType() const {return fTriggerType;}
  Bool_t GetSelectOnTriggerType() {return fSelectOnTriggerType;}

  void SetSelectOnTriggerString( const Bool_t b=kTRUE ) { fSelectOnTriggerString = b; }
  void SetTriggerSelectionString( const TString m );
  TString GetTriggerSelectionString() const { return fTriggerSelectionString; }
  Bool_t GetSelectOnTriggerString() const { return fSelectOnTriggerString; }
  void SetTriggerMaskPatternString( TString str) { fTriggerMaskPatternString = str; }
  TString GetTriggerMaskPatternString() const { return fTriggerMaskPatternString; }

  void SetSelectOnMultiplicity( const Bool_t b = kTRUE) { fSelectOnMultiplicity = b; }
  Bool_t GetSelectOnMultiplicity() { return fSelectOnMultiplicity; }
  void SetMultiplicityLow( const Int_t low ) { fMultiplicityLow = low; }
  void SetMultiplicityHigh( const Int_t high ) { fMultiplicityHigh = high; }
  Int_t GetMultiplicityLow() { return fMultiplicityLow; }
  Int_t GetMultiplicityHigh() { return fMultiplicityHigh; }
  
protected:
  void UpdateEntryList();
  Bool_t CheckOtherSelection(AliESDEvent* ESD);
  
private:
  AliEveEventSelector(const AliEveEventSelector&);
  AliEveEventSelector& operator=(const AliEveEventSelector&);

  AliEveEventManager *fPEventManager; //event manager
  Bool_t fWrapAround;                 //go back to the first event after the last?
  
  Bool_t  fSelectOnString;            //whether to select on string expression
  TString fString;                    //expression to select on
  TEntryList* fPEntryList;            //list of enry numbers after selection
  Int_t  fEntryListId;                //current entry no. in entry list
  Int_t  fLastTreeSize;               //stores the last treesize

  Bool_t  fSelectOnTriggerType;       //whether to select on trigger type
  TString fTriggerType;               //trigger type

  Bool_t  fSelectOnTriggerString;     //whether to select on trigger type with a formula
  TString fTriggerSelectionString;    //trigger selection formula
  TString fTriggerMaskPatternString;  //trigger mask placeholder for evaluation

  Bool_t fSelectOnMultiplicity;       //whether to select on multiplicity
  Int_t  fMultiplicityLow;            //lower bound on multiplicity
  Int_t  fMultiplicityHigh;           //higher bound on multiplicity


  ClassDef(AliEveEventSelector, 2); 
};

#endif

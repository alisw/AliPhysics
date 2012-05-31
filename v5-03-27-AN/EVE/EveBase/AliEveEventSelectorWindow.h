/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//
//  AliEveEventSelectorWindow class
//  GUI for setting event and trigger selections
//
//  origin: Mikolaj Krzewicki, Nikhef, Mikolaj.Krzewicki@cern.ch
//
////////////////////////////////////////////////////////////////////////////


#ifndef AliEveEventSelectorWindow_H
#define AliEveEventSelectorWindow_H

#include "TGedFrame.h"
#include "TGComboBox.h"
#include <TObjString.h>

class AliEveEventManager;
class AliEveEventSelector;
class TGTextEntry;
class TGNumberEntry;
class TGCheckButton;
class TGComboBox;
class TRootEmbeddedCanvas;

//==============================================================================
// AliEveEventSelectorWindow
//==============================================================================

//______________________________________________________________________________
// Short description of AliEveEventSelectorWindow
//

class AliEveEventSelectorWindow : public TGMainFrame
{
public:
  AliEveEventSelectorWindow(const TGWindow *p, UInt_t w, UInt_t h, AliEveEventSelector* sel);
  virtual ~AliEveEventSelectorWindow();
  void SetEventSelector(AliEveEventSelector* sel) {fPSelector = sel;}
  void DoSetSelectionString();
  void DoSetTriggerSelectionString();
  void DoHandleTriggerFromComboBox(const char* str);
  void DoSetMultiplicityRange();
  void DoDrawHistogram();
  void SetupTriggerSelect();

protected:

private:
  AliEveEventSelector* fPSelector; //event selector
  TRootEmbeddedCanvas* fPCanvas;   //the canvas for histograms

  TGTextEntry*   fPDrawFormula;           //test draw input field
  TGTextEntry*   fPEntryFormula;          //selectin formula field
  TGNumberEntry* fPEntryLowerBound;       //lower boung for the formula
  TGNumberEntry* fPEntryHigherBound;      //higher bound for the formula
  TGTextButton*  fPButtonTextDone;        //done button for selection formula
  
  TGComboBox*    fPComboBoxTrigger;       //trigger selection box
  TGTextEntry*   fPEntryTriggerSelection; //trigger selection formula entry field
  TGCheckButton* fPCheckTriggerSimple;    //use simple trigger select
  TGCheckButton* fPCheckTriggerString;    //use trigger select formula
  TGNumberEntry* fPEntryMultHigh;         //lowest allowed multiplicity field
  TGNumberEntry* fPEntryMultLow;          //higest allowed multiplicity field

  AliEveEventSelectorWindow(const AliEveEventSelectorWindow&);
  AliEveEventSelectorWindow& operator=(const AliEveEventSelectorWindow&);
  
  ClassDef(AliEveEventSelectorWindow, 1); // GUI window for AliEveEventSelector
};

#endif

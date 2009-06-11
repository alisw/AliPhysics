/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

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
  AliEveEventSelector* fPSelector;
  TRootEmbeddedCanvas* fPCanvas;

  TGTextEntry* fPDrawFormula;
  TGTextEntry* fPEntryFormula;
  TGNumberEntry* fPEntryLowerBound;
  TGNumberEntry* fPEntryHigherBound;
  TGTextButton* fPButtonTextDone;
  
  TGComboBox* fPComboBoxTrigger;
  TGTextEntry* fPEntryTriggerSelection;
  TGCheckButton* fPCheckTriggerSimple;
  TGCheckButton* fPCheckTriggerString;
  TGNumberEntry* fPEntryMultHigh;
  TGNumberEntry* fPEntryMultLow;

  AliEveEventSelectorWindow(const AliEveEventSelectorWindow&);
  AliEveEventSelectorWindow& operator=(const AliEveEventSelectorWindow&);
  
  ClassDef(AliEveEventSelectorWindow, 1); // GUI window for AliEveEventSelector
};

#endif

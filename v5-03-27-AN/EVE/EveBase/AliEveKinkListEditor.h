// $Id$
// Author: Paraskevi Ganoti: 2009

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveKinkListEditor_H
#define AliEveKinkListEditor_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TEveGDoubleValuator;
class TGComboBox;

class AliEveKinkList;

//______________________________________________________________________________
// Short description of AliEveKinkListEditor
//

class AliEveKinkListEditor : public TGedFrame
{
public:
  AliEveKinkListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                     UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveKinkListEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  void DoMinMaxRCut();
  void DoMinMaxKinkAngleCut();
  void DoMinMaxPt();
  void DoMinMaxInvariantMass();
  void DoSelectDaugPid(Int_t rDaugPid);
  void DoCheckDaugPid();
  void DoSelectDaugProb();

protected:
  AliEveKinkList            *fM; // Model object.

  // Declare widgets
  // TGSomeWidget*   fXYZZ;
  TEveGDoubleValuator* fMinMaxRCut;
  TEveGDoubleValuator* fMinMaxKinkAngleCut; 
  TEveGDoubleValuator* fMinMaxPt;
  TEveGDoubleValuator* fMinMaxInvariantMass;
  TGComboBox*          fDaughterSpecies; 
  TGCheckButton*       fDaughterCheckMaxPidProbability;
  TGNumberEntry*       fDaughterLevelPidProbability;  

private:
  AliEveKinkListEditor(const AliEveKinkListEditor&);            // Not implemented
  AliEveKinkListEditor& operator=(const AliEveKinkListEditor&); // Not implemented

  ClassDef(AliEveKinkListEditor, 0); // GUI editor for AliEveKinkList.
};

#endif

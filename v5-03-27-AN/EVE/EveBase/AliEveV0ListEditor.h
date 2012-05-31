// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveV0ListEditor_H
#define AliEveV0ListEditor_H

#include "TGedFrame.h"

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TEveGDoubleValuator;
class TGComboBox;

class AliEveV0List;

//______________________________________________________________________________
// Short description of AliEveV0ListEditor
//

class AliEveV0ListEditor : public TGedFrame
{
public:
  AliEveV0ListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                     UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveV0ListEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  void DoMinMaxRCut();
  void DoMinMaxDaughterDCA();
  void DoMinMaxPt();
  void DoSelectNegPid(Int_t rNegPid);
  void DoCheckNegPid();
  void DoSelectNegProb();
  void DoSelectPosPid(Int_t rPosPid);
  void DoCheckPosPid();
  void DoSelectPosProb();
  void DoMinMaxInvariantMass();

protected:
  AliEveV0List            *fM; // Model object.

  // Declare widgets
  // TGSomeWidget*   fXYZZ;
  TEveGDoubleValuator* fMinMaxRCut;
  TEveGDoubleValuator* fMinMaxDaughterDCA;
  TEveGDoubleValuator* fMinMaxPt;
  TGComboBox*          fNegativeSpecies;
  TGComboBox*          fPositiveSpecies;
  TGCheckButton*       fNegativeCheckMaxPidProbability;
  TGCheckButton*       fPositiveCheckMaxPidProbability;
  TGNumberEntry*       fNegativeLevelPidProbability;
  TGNumberEntry*       fPositiveLevelPidProbability;
  TEveGDoubleValuator* fMinMaxInvariantMass;

private:
  AliEveV0ListEditor(const AliEveV0ListEditor&);            // Not implemented
  AliEveV0ListEditor& operator=(const AliEveV0ListEditor&); // Not implemented

  ClassDef(AliEveV0ListEditor, 0); // GUI editor for AliEveV0List.
};

#endif

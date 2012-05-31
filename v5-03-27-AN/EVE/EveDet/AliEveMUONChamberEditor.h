// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveMUONChamberEditor_H
#define AliEveMUONChamberEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;
class TGHSlider;

class TEveGValuator;


class AliEveMUONChamber;

class AliEveMUONChamberEditor : public TGedFrame
{
public:
  AliEveMUONChamberEditor(const TGWindow* p = 0,
                          Int_t width = 170, Int_t height = 30,
                          UInt_t options = kChildFrame,
                          Pixel_t back = GetDefaultFrameBackground());

  virtual ~AliEveMUONChamberEditor();

  virtual void SetModel(TObject* obj);

  void DoThreshold();
  void DoMaxVal();
  void DoClusterSize();
  void DoHitSize();

protected:
  AliEveMUONChamber* fM; // fModel dynamic-casted to AliEveMUONChamberEditor

  TEveGValuator *fThreshold;   // digit ADC min
  TEveGValuator *fMaxVal;      // digit ADC max
  TEveGValuator *fClusterSize; // cluster point size
  TEveGValuator *fHitSize;     // hit point size

private:
  AliEveMUONChamberEditor(const AliEveMUONChamberEditor&);            // Not implemented
  AliEveMUONChamberEditor& operator=(const AliEveMUONChamberEditor&); // Not implemented

  ClassDef(AliEveMUONChamberEditor, 0); // Editor for AliEveMUONChamber
};

#endif

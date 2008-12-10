// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTOFSectorEditor_H
#define AliEveTOFSectorEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;

class TGHSlider;

class TEveGValuator;
class TEveGDoubleValuator;

class AliEveTOFSector;

class AliEveTOFSectorEditor : public TGedFrame
{
public:
  AliEveTOFSectorEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
			UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveTOFSectorEditor();

  virtual void SetModel(TObject* obj);
  void DoSectorID();
  void DoAutoTrans();
  void DoPlate0();
  void DoPlate1();
  void DoPlate2();
  void DoPlate3();
  void DoPlate4();

  void DoPlate(Int_t nPlate);
  void DoThreshold();
  void DoMaxVal();

protected:
  AliEveTOFSector*  fM;         // Model object.

  TEveGValuator*    fSectorID;  // Valuator for sector id.

  TGCheckButton*    fAutoTrans; // Check-button for automatic translation.

  TGCheckButton**   fPlate;     // Check-buttons for plates.

  TGCheckButton*    fPlate0;    // Check-button for plate 0.
  TGCheckButton*    fPlate1;    // Check-button for plate 1.
  TGCheckButton*    fPlate2;    // Check-button for plate 2.
  TGCheckButton*    fPlate3;    // Check-button for plate 3.
  TGCheckButton*    fPlate4;    // Check-button for plate 4.

  TEveGValuator*    fThreshold; // Valuator for threshold.
  TEveGValuator*    fMaxVal;    // Valuator for maximum value.

private:
  AliEveTOFSectorEditor(const AliEveTOFSectorEditor&);            // Not implemented
  AliEveTOFSectorEditor& operator=(const AliEveTOFSectorEditor&); // Not implemented

  ClassDef(AliEveTOFSectorEditor, 0); // Editor for AliEveTOFSector
};

#endif

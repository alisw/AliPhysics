// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_TOFSectorEditor_H
#define ALIEVE_TOFSectorEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGDoubleHSlider;

class TGHSlider;

class TEveGValuator;
class TEveGDoubleValuator;
class TEveTransSubEditor;


class AliEveTOFSector;

class AliEveTOFSectorEditor : public TGedFrame
{
private:
  AliEveTOFSectorEditor(const AliEveTOFSectorEditor&);            // Not implemented
  AliEveTOFSectorEditor& operator=(const AliEveTOFSectorEditor&); // Not implemented

protected:
  AliEveTOFSector*  fM; // Model object.

  TEveGValuator*    fSectorID;

  TGCheckButton*    fAutoTrans;

  TGCheckButton**   fPlate;

  TGCheckButton*    fPlate0;
  TGCheckButton*    fPlate1;
  TGCheckButton*    fPlate2;
  TGCheckButton*    fPlate3;
  TGCheckButton*    fPlate4;

  TEveGValuator*    fThreshold;
  TEveGValuator*    fMaxVal;

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

  ClassDef(AliEveTOFSectorEditor, 0); // Editor for AliEveTOFSector
}; // endclass AliEveTOFSectorEditor

#endif

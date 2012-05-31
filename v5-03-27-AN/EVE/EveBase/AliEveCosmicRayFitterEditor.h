// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_CosmicRayFitterEditor_H
#define ALIEVE_CosmicRayFitterEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class AliEveCosmicRayFitter;

class AliEveCosmicRayFitterEditor : public TGedFrame
{
public:
  AliEveCosmicRayFitterEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveCosmicRayFitterEditor() {}

  virtual void SetModel(TObject* obj);

  void DoStart();
  void DoFit();
  void DoReset();
  void DoStop();
  void DoGraph();

protected:
  AliEveCosmicRayFitter* fM; // fModel dynamic-casted to AliEveCosmicRayFitterEditor

  TGTextButton* fFit;    // button to fit selection
  TGTextButton* fReset;  // button to reset selection
  TGTextButton* fStart;  // button to connect to signal
  TGTextButton* fStop;   // button to disconnect from signal
  TGTextButton* fGraph;  // button to draw graph

private:
  AliEveCosmicRayFitterEditor(const AliEveCosmicRayFitterEditor&);            // Not implemented
  AliEveCosmicRayFitterEditor& operator=(const AliEveCosmicRayFitterEditor&); // Not implemented

  ClassDef(AliEveCosmicRayFitterEditor, 0); // Editor for AliEveCosmicRayFitter class.
};

#endif

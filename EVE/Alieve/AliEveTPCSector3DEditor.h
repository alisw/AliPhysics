// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

#ifndef ALIEVE_TPCSector3DEditor_H
#define ALIEVE_TPCSector3DEditor_H

#include <Alieve/AliEveTPCSector2DEditor.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class TEveGValuator;
class TEveGDoubleValuator;


class AliEveTPCSector3D;

class AliEveTPCSector3DEditor : public TGedFrame
{
  AliEveTPCSector3DEditor(const AliEveTPCSector3DEditor&);            // Not implemented
  AliEveTPCSector3DEditor& operator=(const AliEveTPCSector3DEditor&); // Not implemented

protected:
  AliEveTPCSector3D*      fM; // fModel dynamic-casted to AliEveTPCSector3DEditor

  TGCheckButton*    fRnrFrame;
  TEveGValuator* fDriftVel;

  TEveGValuator* fPointFrac;
  TEveGValuator* fPointSize;

public:
  AliEveTPCSector3DEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		    UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());  
  virtual ~AliEveTPCSector3DEditor();

  virtual void SetModel(TObject* obj);

  void DoRnrFrame();
  void DoDriftVel();

  void DoPointFrac();
  void DoPointSize();

  ClassDef(AliEveTPCSector3DEditor, 0); // Editor for AliEveTPCSector3D
}; // endclass AliEveTPCSector3DEditor

#endif

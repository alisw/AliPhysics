// $Header$

#ifndef ALIEVE_TPCSector3DEditor_H
#define ALIEVE_TPCSector3DEditor_H

#include <Alieve/TPCSector2DEditor.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {
class RGValuator;
class RGDoubleValuator;
}

namespace Alieve {

class TPCSector3D;

class TPCSector3DEditor : public TGedFrame
{
protected:
  TPCSector3D*      fM; // fModel dynamic-casted to TPCSector3DEditor

  TGCheckButton*    fRnrFrame;
  Reve::RGValuator* fDriftVel;

  Reve::RGValuator* fPointFrac;

public:
  TPCSector3DEditor(const TGWindow* p, Int_t id, Int_t width=170, Int_t height=30,
		    UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());  
  virtual ~TPCSector3DEditor();

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);

  void DoRnrFrame();
  void DoDriftVel();

  void DoPointFrac();

  ClassDef(TPCSector3DEditor, 0); // Editor for TPCSector3D
}; // endclass TPCSector3DEditor

}

#endif

// $Header$

#ifndef ALIEVE_AliEVEHOMERManagerEditor_H
#define ALIEVE_AliEVEHOMERManagerEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGTextButton;
class TGNumberEntry;
class TGColorSelect;

class AliEVEHOMERManager;

class AliEVEHOMERManagerEditor : public TGedFrame
{
private:
  AliEVEHOMERManagerEditor(const AliEVEHOMERManagerEditor&);            // Not implemented
  AliEVEHOMERManagerEditor& operator=(const AliEVEHOMERManagerEditor&); // Not implemented

protected:
  AliEVEHOMERManager* fM; // fModel dynamic-casted to AliEVEHOMERManagerEditor

  // Declare widgets
  // TGSomeWidget*   fXYZZ;
  TGTextButton  *fButt;

public:
  AliEVEHOMERManagerEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEVEHOMERManagerEditor();

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();
  void DoButt();

  ClassDef(AliEVEHOMERManagerEditor, 0); // Editor for AliEVEHOMERManager
}; // endclass AliEVEHOMERManagerEditor

#endif

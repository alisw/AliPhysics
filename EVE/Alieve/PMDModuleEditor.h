// $Header$

#ifndef ALIEVE_PMDModuleEditor_H
#define ALIEVE_PMDModuleEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Alieve {

class PMDModule;

class PMDModuleEditor : public TGedFrame
{
private:
  PMDModuleEditor(const PMDModuleEditor&);            // Not implemented
  PMDModuleEditor& operator=(const PMDModuleEditor&); // Not implemented

  void CreateInfoFrame();

protected:
  PMDModule* fM; // fModel dynamic-casted to PMDModuleEditor

  TGVerticalFrame*  fInfoFrame;

  TGLabel*   fInfoLabel0;
  TGLabel*   fInfoLabel1;
  TGLabel*   fInfoLabel2;
  TGLabel*   fInfoLabel3;
  TGLabel*   fInfoLabel4;
  TGLabel*   fInfoLabel5;

public:
  PMDModuleEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~PMDModuleEditor();

  virtual void SetModel(TObject* obj);
  void DisplayHistos();
  //  void PrintADC();



  // Declare callback/slot methods
  // void DoXYZZ();

  ClassDef(PMDModuleEditor, 0); // Editor for PMDModule
}; // endclass PMDModuleEditor

}

#endif

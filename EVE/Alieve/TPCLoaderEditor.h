// $Header$

#ifndef ALIEVE_TPCLoaderEditor_H
#define ALIEVE_TPCLoaderEditor_H

#include <TGedFrame.h>

class TGTextButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGTextEntry;

namespace Reve {
class RGValuator;
}

namespace Alieve {

class TPCLoader;

class TPCLoaderEditor : public TGedFrame
{
  TPCLoaderEditor(const TPCLoaderEditor&);            // Not implemented
  TPCLoaderEditor& operator=(const TPCLoaderEditor&); // Not implemented

protected:
  TPCLoader* fM; // fModel dynamic-casted to TPCLoaderEditor

  TGTextEntry*  fFile;
  TGTextButton* fOpenFile;

  Reve::RGValuator* fEvent;
  TGCheckButton*    fDoubleSR;

  TGTextButton* fUpdateSectors;
  TGTextButton* fCreateSectors3D;
  TGTextButton* fDeleteSectors3D;

public:
  TPCLoaderEditor(const TGWindow* p=0, Int_t id=-1, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  ~TPCLoaderEditor();

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);

  void FileSelect();
  void FileChanged();
  void DoOpen();

  void DoEvent();
  void DoDoubleSR();

  void DoUpdateSectors();
  void DoCreateSectors3D();
  void DoDeleteSectors3D();

  ClassDef(TPCLoaderEditor, 0); // Editor for TPCLoader
}; // endclass TPCLoaderEditor

}

#endif

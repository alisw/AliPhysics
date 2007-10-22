// $Header$

#ifndef ALIEVE_ITSScaledModuleEditor_H
#define ALIEVE_ITSScaledModuleEditor_H

#include <TGedFrame.h>
#include <Reve/RGBAPaletteEditor.h>

class TGNumberEntry;
class TGColorSelect;
class TGComboBox;

namespace Reve 
{
class RGValuator;
class RGDoubleValuator;
class RGBAPalette;
}

namespace Alieve {

class DigitScaleInfo;
class ITSScaledModule;
class AliITSsegmentation;

/**************************************************************************/

class ITSScaledModuleEditor : public TGedFrame
{
private:
  ITSScaledModuleEditor(const ITSScaledModuleEditor&);            // Not implemented
  ITSScaledModuleEditor& operator=(const ITSScaledModuleEditor&); // Not implemented

  void CreateInfoFrame();

protected:
  TGVerticalFrame*  fInfoFrame;

  ITSScaledModule*  fModule; // fModel dynamic-casted to ITSScaledModuleEditor

  TGNumberEntry*    fScale;
  TGComboBox*       fStatistic;  

  TGLabel*          fInfoLabel0;
  TGLabel*          fInfoLabel1;

public:
  ITSScaledModuleEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~ITSScaledModuleEditor();

  virtual void   SetModel(TObject* obj);

  void DoScale();
  void DoStatType(Int_t t);

  ClassDef(ITSScaledModuleEditor, 0); // Editor for ITSScaledModule
}; // endclass ITSScaledModuleEditor

}

#endif

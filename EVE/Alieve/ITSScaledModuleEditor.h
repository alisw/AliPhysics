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
class ZTransSubEditor;
}

namespace Alieve {

class DigitScaleInfo;
class ITSScaledModule;
class AliITSsegmentation;

class ITSSDSubEditor : public Reve::RGBAPaletteSubEditor
{
private:
  ITSSDSubEditor(const ITSSDSubEditor&);            // Not implemented
  ITSSDSubEditor& operator=(const ITSSDSubEditor&); // Not implemented
  
  void GetSubDetScaleData(Int_t& cnx, Int_t& cnz, Int_t& total, Float_t& maxoc);
  void SetPaletteFromDigitInfo();
protected:
  ITSScaledModule*  fModule;

  TGNumberEntry*    fScale;
  TGComboBox*       fStatistic;  

  TGLabel*          fInfoLabel0;
  TGLabel*          fInfoLabel1;
  TGLabel*          fInfoLabel2;

public:
  ITSSDSubEditor(const TGWindow* p);
  virtual ~ITSSDSubEditor() {}

  void SetModel(ITSScaledModule* sm);

  void DoScale();
  void DoStatType(Int_t t);

  ClassDef(ITSSDSubEditor, 1); // SubEditor for RGBAPalelet and scaled digits
};  // endclass ITSSDPaletteSubEditor

/**************************************************************************/

class ITSScaledModuleEditor : public TGedFrame
{
private:
  ITSScaledModuleEditor(const ITSScaledModuleEditor&);            // Not implemented
  ITSScaledModuleEditor& operator=(const ITSScaledModuleEditor&); // Not implemented

protected:
  ITSScaledModule*       fM; // fModel dynamic-casted to ITSScaledModuleEditor

  Reve::ZTransSubEditor* fHMTrans;
  ITSSDSubEditor*        fSDPalette;
  
public:
  ITSScaledModuleEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~ITSScaledModuleEditor();

  virtual void   ActivateBaseClassEditors(TClass* cl);

  virtual void   SetModel(TObject* obj);

  ClassDef(ITSScaledModuleEditor, 0); // Editor for ITSScaledModule
}; // endclass ITSScaledModuleEditor

}

#endif

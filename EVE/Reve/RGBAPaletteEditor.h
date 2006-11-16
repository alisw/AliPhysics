// $Header$

#ifndef REVE_RGBAPaletteEditor_H
#define REVE_RGBAPaletteEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGColorSelect;
class TGComboBox;


namespace Reve {

class RGBAPalette;
class RGValuator;
class RGDoubleValuator;

class RGBAPaletteSubEditor : public TGVerticalFrame
{
private:
  RGBAPaletteSubEditor(const RGBAPaletteSubEditor&);            // Not implemented
  RGBAPaletteSubEditor& operator=(const RGBAPaletteSubEditor&); // Not implemented

protected:
  RGBAPalette*      fM;

  TGComboBox*       fUndershootAction;
  TGColorSelect*    fUnderColor;
  TGComboBox*       fOvershootAction;
  TGColorSelect*    fOverColor;

  RGDoubleValuator* fMinMax;

  TGCheckButton*    fInterpolate;
  TGCheckButton*    fShowDefValue;
  TGColorSelect*    fDefaultColor;

public:
  RGBAPaletteSubEditor(const TGWindow* p);
  virtual ~RGBAPaletteSubEditor() {}

  void SetModel(RGBAPalette* p);

  void Changed(); //*SIGNAL*

  void DoMinMax();

  void DoInterpolate();
  void DoShowDefValue();
  void DoDefaultColor(Pixel_t color);
  void DoUnderColor(Pixel_t color);
  void DoOverColor(Pixel_t color);
  void DoUndershootAction(Int_t mode);
  void DoOvershootAction(Int_t mode);

  ClassDef(RGBAPaletteSubEditor, 1); // SubEditor for RGBAPalette
}; // endclass RGBAPaletteSubEditor


/**************************************************************************/
/**************************************************************************/

class RGBAPaletteEditor : public TGedFrame
{
private:
  RGBAPaletteEditor(const RGBAPaletteEditor&);            // Not implemented
  RGBAPaletteEditor& operator=(const RGBAPaletteEditor&); // Not implemented

protected:
  RGBAPalette*          fM;
  RGBAPaletteSubEditor* fSE;

public:
  RGBAPaletteEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~RGBAPaletteEditor();

  virtual void SetModel(TObject* obj);

  ClassDef(RGBAPaletteEditor, 1); // Editor for RGBAPalette
}; // endclass RGBAPaletteEditor

}

#endif

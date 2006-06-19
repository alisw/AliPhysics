// $Header$

#ifndef REVE_RGValuators_H
#define REVE_RGValuators_H

#include <TGNumberEntry.h>

class TGLabel;
class TGHSlider;
class TGDoubleHSlider;


namespace Reve {

class RGValuatorBase: public TGCompositeFrame
{
protected:
  UInt_t      fLabelWidth;
  Bool_t      fAlignRight;
  Bool_t      fShowSlider;

  Int_t       fNELength;
  Int_t       fNEHeight;

  TGLabel*    fLabel;

public:
  RGValuatorBase(const TGWindow *p, const char* title, UInt_t w, UInt_t h);
  virtual ~RGValuatorBase() {}

  virtual void Build() = 0;

  void SetLabelWidth(Int_t w)        { fLabelWidth = w; }
  void SetAlignRight(Bool_t a)       { fAlignRight = a; }
  void SetShowSlider(Bool_t s=kTRUE) { fShowSlider = s; }

  void SetNELength(Int_t l)          { fNELength = l; }
  void SetNEGeight(Int_t h)          { fNEHeight = h; }

  ClassDef(RGValuatorBase, 0);
}; // endclass RGValuatorBase

/**************************************************************************/

class RGValuator: public RGValuatorBase
{
protected:
  Float_t        fValue;
  Float_t        fMin;
  Float_t        fMax;

  Bool_t         fSliderNewLine;
  Int_t          fSliderDivs;
  TGNumberEntry* fEntry;
  TGHSlider*     fSlider;

  Int_t CalcSliderPos(Float_t v);
  
public:
  RGValuator(const TGWindow *p, const char* title, UInt_t w, UInt_t h);
  virtual ~RGValuator() {}

  virtual void Build();
  
  Float_t GetValue() const { return fValue; }
  virtual void SetValue(Float_t v, Bool_t emit=kFALSE);

  void SliderCallback();
  void EntryCallback();
  void ValueSet(Double_t); //*SIGNAL* 

  TGHSlider*     GetSlider() { return fSlider; }
  TGNumberEntry* GetEntry()  { return fEntry; }

  void SetSliderNewLine(Bool_t nl) { fSliderNewLine = nl; }

  void SetLimits(Int_t min, Int_t max);
  void SetLimits(Float_t min, Float_t max, Int_t npos,
		 TGNumberFormat::EStyle nef=TGNumberFormat::kNESRealTwo);

  void SetToolTip(const Text_t* tip);
  void SetEnabled(Bool_t state);

  ClassDef(RGValuator, 0);
}; // endclass RGValuator

/**************************************************************************/

class RGDoubleValuator: public RGValuatorBase
{
protected:
  TGNumberEntry*    fMinEntry;
  TGNumberEntry*    fMaxEntry;
  TGDoubleHSlider*  fSlider;
  
public:
  RGDoubleValuator(const TGWindow *p, const char* title, UInt_t w, UInt_t h);
  virtual ~RGDoubleValuator() {}

  virtual void Build();
 
  void MinEntryCallback();
  void MaxEntryCallback();
  void SliderCallback();
  void ValueSet(); //*SIGNAL* 

  TGDoubleHSlider* GetSlider()   { return fSlider; }
  TGNumberEntry*   GetMinEntry() { return fMinEntry; }
  TGNumberEntry*   GetMaxEntry() { return fMaxEntry; }

  void SetLimits(Float_t min, Float_t max, TGNumberFormat::EStyle nef=TGNumberFormat::kNESRealTwo);
  void SetValues(Float_t min, Float_t max, Bool_t emit=kFALSE);

  void GetValues(Float_t& min, Float_t& max) const
  { min = fMinEntry->GetNumber(); max = fMaxEntry->GetNumber(); }
  Float_t GetMin() const { return fMinEntry->GetNumber(); }
  Float_t GetMax() const { return fMaxEntry->GetNumber(); }

  ClassDef(RGDoubleValuator, 0);
}; // endclass RGDoubleValuator

}

#endif

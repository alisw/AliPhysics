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
  RGValuatorBase(const RGValuatorBase&);            // Not implemented
  RGValuatorBase& operator=(const RGValuatorBase&); // Not implemented

protected:
  UInt_t      fLabelWidth;
  Bool_t      fAlignRight;
  Bool_t      fShowSlider;

  Int_t       fNELength; // Number-entry length (in characters)
  Int_t       fNEHeight; // Number-entry height (in pixels)

  TGLabel*    fLabel;

public:
  RGValuatorBase(const TGWindow *p, const char* title, UInt_t w, UInt_t h);
  virtual ~RGValuatorBase() {}

  virtual void Build(Bool_t connect=kTRUE) = 0;

  void SetLabelWidth(Int_t w)        { fLabelWidth = w; }
  void SetAlignRight(Bool_t a)       { fAlignRight = a; }
  void SetShowSlider(Bool_t s=kTRUE) { fShowSlider = s; }

  void SetNELength(Int_t l)          { fNELength = l; }
  void SetNEHeight(Int_t h)          { fNEHeight = h; }

  ClassDef(RGValuatorBase, 0);
}; // endclass RGValuatorBase

/**************************************************************************/

class RGValuator: public RGValuatorBase
{
  RGValuator(const RGValuator&);            // Not implemented
  RGValuator& operator=(const RGValuator&); // Not implemented

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

  virtual void Build(Bool_t connect=kTRUE);
  
  Float_t GetValue() const { return fValue; }
  virtual void SetValue(Float_t v, Bool_t emit=kFALSE);

  void SliderCallback();
  void EntryCallback();
  void ValueSet(Double_t); //*SIGNAL* 

  TGHSlider*     GetSlider() { return fSlider; }
  TGNumberEntry* GetEntry()  { return fEntry; }

  void SetSliderNewLine(Bool_t nl) { fSliderNewLine = nl; }

  void GetLimits(Float_t& min, Float_t& max) const { min = fMin; max = fMax; }
  Float_t GetLimitMin() const { return fMin; }
  Float_t GetLimitMax() const { return fMax; }
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
  RGDoubleValuator(const RGDoubleValuator&);            // Not implemented
  RGDoubleValuator& operator=(const RGDoubleValuator&); // Not implemented

protected:
  TGNumberEntry*    fMinEntry;
  TGNumberEntry*    fMaxEntry;
  TGDoubleHSlider*  fSlider;
  
public:
  RGDoubleValuator(const TGWindow *p, const char* title, UInt_t w, UInt_t h);
  virtual ~RGDoubleValuator() {}

  virtual void Build(Bool_t connect=kTRUE);
 
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

/**************************************************************************/

class RGTriVecValuator : public TGCompositeFrame
{
  RGTriVecValuator(const RGTriVecValuator&);            // Not implemented
  RGTriVecValuator& operator=(const RGTriVecValuator&); // Not implemented

protected:
  RGValuator* fVal[3];

  // Weed-size vars from RGValuator; copied.
  UInt_t      fLabelWidth;
  Int_t       fNELength; // Number-entry length (in characters)
  Int_t       fNEHeight; // Number-entry height (in pixels)

public:
  RGTriVecValuator(const TGWindow *p, const char* name, UInt_t w, UInt_t h);
  virtual ~RGTriVecValuator();

  void Build(Bool_t vertical, const char* lab0, const char* lab1, const char* lab2);

  RGValuator* GetValuator(Int_t i) const { return fVal[i]; }

  Float_t GetValue(Int_t i) const   { return fVal[i]->GetValue(); }
  void SetValue(Int_t i, Float_t v) { fVal[i]->SetValue(v); }

  void GetValues(Float_t& v0, Float_t& v1, Float_t& v2) const
  { v0 = GetValue(0); v1 = GetValue(1); v2 = GetValue(2); }
  void GetValues(Float_t v[3]) const
  { v[0] = GetValue(0); v[1] = GetValue(1); v[2] = GetValue(2); }
  void GetValues(Double_t v[3]) const
  { v[0] = GetValue(0); v[1] = GetValue(1); v[2] = GetValue(2); }

  void SetValues(Float_t v0, Float_t v1, Float_t v2)
  { SetValue(0, v0); SetValue(1, v1); SetValue(2, v2); }
  void SetValues(Float_t v[3])
  { SetValue(0, v[0]); SetValue(1, v[1]); SetValue(2, v[2]); }
  void SetValues(Double_t v[3])
  { SetValue(0, v[0]); SetValue(1, v[1]); SetValue(2, v[2]); }

  void ValueSet(); //*SIGNAL*
  
  // Weed-size vars from RGValuator; copied.
  void SetLabelWidth(Int_t w) { fLabelWidth = w; }
  void SetNELength(Int_t l)   { fNELength   = l; }
  void SetNEHeight(Int_t h)   { fNEHeight   = h; }

  void SetLimits(Int_t min, Int_t max);
  void SetLimits(Float_t min, Float_t max,
		 TGNumberFormat::EStyle nef=TGNumberFormat::kNESRealTwo);

  ClassDef(RGTriVecValuator, 0)
};

}

#endif

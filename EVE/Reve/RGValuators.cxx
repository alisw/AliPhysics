// $Header$

#include "RGValuators.h"

#include <TGLabel.h>
#include <TGSlider.h>
#include <TGDoubleSlider.h>

using namespace Reve;

/**************************************************************************/
// RGValuatorBase
/**************************************************************************/

ClassImp(RGValuatorBase)

RGValuatorBase::RGValuatorBase(const TGWindow *p, const char* name,
			       UInt_t w, UInt_t h) :
  TGCompositeFrame(p, w, h),

  fLabelWidth (0),
  fAlignRight (kFALSE),
  fShowSlider (kTRUE),

  fNELength (5),
  fNEHeight (20),

  fLabel (0)
{
   SetName(name);
}

/**************************************************************************/
// RGValuator
/**************************************************************************/

ClassImp(RGValuator)

RGValuator::RGValuator(const TGWindow *p, const char* title,
		       UInt_t w, UInt_t h) :
  RGValuatorBase(p, title, w, h),

  fValue (0),
  fMin   (0),
  fMax   (0),

  fSliderNewLine (kFALSE),
  fSliderSteps   (-1),
  fEntry  (0),
  fSlider (0)
{}

void RGValuator::Build()
{
  TGCompositeFrame *hf1, *hfs;
  if(fShowSlider && fSliderNewLine) {
    SetLayoutManager(new TGVerticalLayout(this));
    hf1 = new TGHorizontalFrame(this);
    hf1->SetLayoutManager(new TGHorizontalLayout(hf1));
    AddFrame(hf1, new TGLayoutHints(kLHintsTop, 0,0,0,0));
    hfs = new TGHorizontalFrame(this);
    hfs->SetLayoutManager(new TGHorizontalLayout(hfs));
    AddFrame(hfs, new TGLayoutHints(kLHintsTop, 0,0,0,0));
  } else {
    hf1 = this;
    hfs = this;
    SetLayoutManager(new TGHorizontalLayout(this));
  }

  // label
  TGLayoutHints* lh;
  if(fAlignRight)
    lh = new TGLayoutHints(kLHintsRight | kLHintsBottom, 4,0,0,0);
  else
    lh = new TGLayoutHints(kLHintsLeft  | kLHintsBottom, 0,4,0,0);
  
  if (fLabelWidth > 0) {
    // printf("fLabelWidth > 0 \n");
    TGCompositeFrame *lf = new TGHorizontalFrame(hf1, fLabelWidth, fNEHeight, kFixedSize);
    fLabel = new TGLabel(lf, GetName());
    lf->AddFrame(fLabel, lh); 
    // add label frame to top horizontal frame
    TGLayoutHints* lfh = new TGLayoutHints(kLHintsLeft, 0,0,0,0);
    hf1->AddFrame(lf, lfh);
  } else {
    fLabel = new TGLabel(hf1, GetName());
    hf1->AddFrame(fLabel, lh);  
  }

  // entry
  TGLayoutHints*  elh =  new TGLayoutHints(kLHintsLeft, 0,0,0,0);
  fEntry = new TGNumberEntry(hf1, 0, fNELength);
  fEntry->SetHeight(fNEHeight);
  fEntry->GetNumberEntry()->SetToolTipText("Enter Slider Value");
  hf1->AddFrame(fEntry, elh);

  fEntry->Associate(this);   
  fEntry->Connect("ValueSet(Long_t)",
		  "Reve::RGValuator", this, "EntryCallback()");
  
  
  // slider
  if(fShowSlider) {
    fSlider = new TGHSlider(hfs, GetWidth(), kSlider1 | kScaleBoth);
    hfs->AddFrame(fSlider, new TGLayoutHints(kLHintsLeft, 1,1,0,0));
   
    fSlider->Associate(this);
    fSlider->Connect("PositionChanged(Int_t)",
		     "Reve::RGValuator", this, "SliderCallback()");
  }
}

void RGValuator::SetLimits(Float_t min, Float_t max, Int_t ndiv,
			   TGNumberFormat::EStyle nef) 
{
  fMin = Float_t(min);
  fMax = Float_t(max);
  fSliderSteps = ndiv;
  fEntry->SetFormat(nef);
  fEntry->SetLimits(TGNumberFormat::kNELLimitMinMax, min, max);
}

void RGValuator::SetLimits(Int_t min, Int_t max) 
{
  fMin = Float_t(min);
  fMax = Float_t(max);
  fEntry->SetFormat(TGNumberFormat::kNESInteger);
  fEntry->SetLimits(TGNumberFormat::kNELLimitMinMax, min, max);

  if(fSlider){
    fSlider->SetRange(min, max);
    fSliderSteps = max - min;
  }
}

void RGValuator::EntryCallback()
{
  fValue = fEntry->GetNumber();
  if(fSlider) {
    Int_t pos;
    if (fSliderSteps != -1) {
      pos = Int_t( fSliderSteps*(fValue - fMin)/(fMax - fMin));
    } else {
      pos = Int_t(fValue - fMin);
    }
    pos += Int_t(fMin);
    fSlider->SetPosition(pos);
    // printf( "RGValuator::EntryCallback() slider pos %d n", pos);
  }
  ValueSet(fValue);
}

void RGValuator::SliderCallback()
{
  Double_t val;
  if(fSliderSteps != -1)
    val = fMin + fSlider->GetPosition()*Double_t((fMax-fMin))/fSliderSteps;
  else 
    val = Double_t(fSlider->GetPosition());
  
  fValue = val;
  fEntry->SetNumber(fValue);
  ValueSet(fValue);
}


void RGValuator::ValueSet(Double_t val)
{
  Emit("ValueSet(Double_t)", val);
}

void RGValuator::SetValue(Float_t val, Bool_t emit)
{
  fValue = val;
  fEntry->SetNumber(fValue);

  if(fSlider){
    fSlider->SetPosition(Int_t((val-fMin)*fSliderSteps/(fMax-fMin)));
    // printf("RGValuator::ValueSet slider pos %d\n",fSlider->GetPosition() );
  }
  if(emit)
    ValueSet(val);
}

void RGValuator::SetToolTip(const Text_t* tip)
{
  fEntry->GetNumberEntry()->SetToolTipText(tip);
}

void RGValuator::SetEnabled(Bool_t state)
{
  fEntry->GetNumberEntry()->SetEnabled(state);
  fEntry->GetButtonUp()->SetEnabled(state);
  fEntry->GetButtonDown()->SetEnabled(state);
  if(fSlider) {
    if(state) fSlider->MapWindow();
    else      fSlider->UnmapWindow();
  }
}

/**************************************************************************/
// RGDoubleValuator
/**************************************************************************/

ClassImp(RGDoubleValuator)

RGDoubleValuator::RGDoubleValuator(const TGWindow *p, const char* title,
				   UInt_t w, UInt_t h) :
  RGValuatorBase(p, title, w, h),

  fMinEntry(0),
  fMaxEntry(0),
  fSlider(0)
{}

void RGDoubleValuator::Build()
{
  TGCompositeFrame *hf1, *hfs;
  if(fShowSlider) {
    SetLayoutManager(new TGVerticalLayout(this));
    hf1 = new TGHorizontalFrame(this);
    hf1->SetLayoutManager(new TGHorizontalLayout(hf1));
    AddFrame(hf1, new TGLayoutHints(kLHintsTop, 0,0,0,0));
    hfs = new TGHorizontalFrame(this);
    hfs->SetLayoutManager(new TGHorizontalLayout(hfs));
    AddFrame(hfs, new TGLayoutHints(kLHintsTop, 0,0,0,0));
  } else {
    hf1 = this;
    hfs = this;
    SetLayoutManager(new TGHorizontalLayout(this));
  }

  // label
  TGLayoutHints* lh;
  if(fAlignRight)
    lh = new TGLayoutHints(kLHintsRight | kLHintsBottom, 4,0,0,0);
  else
    lh = new TGLayoutHints(kLHintsLeft  | kLHintsBottom, 0,4,0,0);
  
  if(fLabelWidth > 0) {
    TGCompositeFrame *lf = new TGHorizontalFrame(hf1, fLabelWidth, fNEHeight, kFixedSize);
    fLabel = new TGLabel(lf, GetName());
    lf->AddFrame(fLabel, lh); 
    // add label frame to top horizontal frame
    TGLayoutHints* lfh = new TGLayoutHints(kLHintsLeft, 0,0,0,0);
    hf1->AddFrame(lf, lfh);
  } else {
    fLabel = new TGLabel(hf1, GetName());
    hf1->AddFrame(fLabel, lh);  
  }

  // entries
  fMinEntry = new TGNumberEntry(this, 0, fNELength);
  fMinEntry->SetHeight(fNEHeight);
  fMinEntry->GetNumberEntry()->SetToolTipText("Enter Slider Min Value");
  hf1->AddFrame(fMinEntry, new TGLayoutHints(kLHintsLeft, 0,0,0,0));
  fMinEntry->Connect("ValueSet(Long_t)",
		     "Reve::RGDoubleValuator", this, "MinEntryCallback()");
  fMinEntry->Associate(this);   
   
  fMaxEntry = new TGNumberEntry(this, 0, fNELength);
  fMaxEntry->SetHeight(fNEHeight);
  fMaxEntry->GetNumberEntry()->SetToolTipText("Enter Slider Max Value");
  hf1->AddFrame(fMaxEntry,  new TGLayoutHints(kLHintsLeft, 2,0,0,0));
  fMaxEntry->Connect("ValueSet(Long_t)",
		     "Reve::RGDoubleValuator", this, "MaxEntryCallback()");
  fMaxEntry->Associate(this);   
  
  // slider
  if(fShowSlider) {
    fSlider = new TGDoubleHSlider(hfs, GetWidth(), kDoubleScaleBoth);
    hfs->AddFrame(fSlider, new TGLayoutHints(kLHintsTop|kLHintsLeft, 0,0,1,0));
    fSlider->Associate(this);
    fSlider->Connect("PositionChanged()",
		     "Reve::RGDoubleValuator", this, "SliderCallback()");
  }
}

void RGDoubleValuator::SetLimits(Float_t min, Float_t max,
				 TGNumberFormat::EStyle nef) 
{
  //  printf("RGDoubleValuator::SetLimits(Float_t min, Float_t max, Int_ \n");
  fMinEntry->SetLimits(TGNumberFormat::kNELLimitMinMax, min, max);
  fMinEntry->SetFormat(nef);
  fMaxEntry->SetLimits(TGNumberFormat::kNELLimitMinMax, min, max);
  fMaxEntry->SetFormat(nef);
  
  if(fSlider) fSlider->SetRange(min, max);
}

void RGDoubleValuator::MinEntryCallback()
{
  if(GetMin() > GetMax())
    fMaxEntry->SetNumber(GetMin());
  if(fSlider) fSlider->SetPosition(GetMin(), GetMax());
  ValueSet();
}

void RGDoubleValuator::MaxEntryCallback()
{
  if(GetMax() < GetMin())
    fMinEntry->SetNumber(GetMax());
  if(fSlider) fSlider->SetPosition(GetMin(), GetMax());
  ValueSet();
}

void RGDoubleValuator::SliderCallback()
{
  Float_t minp, maxp;
  fSlider->GetPosition(minp, maxp);
  //printf("RGDoubleValuator::SliderCallback %f %f\n", minp, maxp);
  fMinEntry->SetNumber(minp);
  fMaxEntry->SetNumber(maxp); 
  ValueSet();
}

void RGDoubleValuator::SetValues(Float_t min, Float_t max, Bool_t emit)
{
  fMinEntry->SetNumber(min);
  fMaxEntry->SetNumber(max);

  if(fSlider) fSlider->SetPosition(min, max);
  if(emit)    ValueSet();
}

void RGDoubleValuator::ValueSet()
{
  Emit("ValueSet()");
}

/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/////////////////////////////////////////////////////////////////////////
// ALICE SLIDER FRAME CLASS                                            //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <TGFrame.h>
#include <TGLayout.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TEnv.h>

#include "AliSliderFrame.h"
#include "AliDisplay2.h"

ClassImp(AliSliderFrame)

//_____________________________________________________________
AliSliderFrame::AliSliderFrame(const TGWindow *p, UInt_t w, UInt_t h)
{
  // Constructor
  fMainFrame = new TGCompositeFrame(p, w, h, kVerticalFrame | kRaisedFrame);
  fLayout = new TGLayoutHints( kLHintsBottom | kLHintsRight | kLHintsExpandX,5,5,2,2);
  
  //Momentum Slider
  fMomentumFrame = new TGCompositeFrame(fMainFrame,0,0, kHorizontalFrame);
  fMomentumLayout = new TGLayoutHints( kLHintsLeft | kLHintsTop | kLHintsExpandX ,5,5,2,2);
  fMomentumSlider = new TGDoubleHSlider(fMomentumFrame, 400, kDoubleScaleBoth, kIdsMOMENTUM);
  fMomentumSlider->Connect("PositionChanged()","AliSliderFrame",this,"DoSlider()");
  fMomentumSlider->Connect("PositionChanged()","AliSliderFrame",this,"DoPositionChanged(Int_t)");
  fMomentumSlider->Connect("Released()","AliSliderFrame",this,"DoReleased()");
  SetMomentumRange(0,2);
  fMomentumSlider->SetPosition(0,2);
  fMomentumLabel = new TGLabel(fMomentumFrame,"Momentum");
  
  fMomentumMinValue = new TGNumberEntry(fMomentumFrame,GetMomentumMin(),7,kIdsMomentumMIN);
  fMomentumMinValue->Connect("ValueChanged(Long_t)","AliSliderFrame",this,"DoField(Long_t)");
  fMomentumMinValue->GetNumberEntry()->Connect("ReturnPressed()","AliSliderFrame",this,"DoField(Long_t)");
  fMomentumMinValue->SetButtonToNum(kFALSE);
  
  fMomentumMaxValue = new TGNumberEntry(fMomentumFrame,GetMomentumMax(),7,kIdsMomentumMAX);
  fMomentumMaxValue->Connect("ValueChanged(Long_t)","AliSliderFrame",this,"DoField(Long_t)");
  fMomentumMaxValue->GetNumberEntry()->Connect("ReturnPressed()","AliSliderFrame",this,"DoField(Long_t)");
  fMomentumMaxValue->SetButtonToNum(kFALSE);
  
  fMomentumFrame->AddFrame(fMomentumLabel,new TGLayoutHints( kLHintsLeft | kLHintsCenterY,5,5,2,2));
  fMomentumFrame->AddFrame(fMomentumMinValue,new TGLayoutHints( kLHintsLeft ,5,5,2,2));
  fMomentumFrame->AddFrame(fMomentumSlider,fMomentumLayout);
  fMomentumFrame->AddFrame(fMomentumMaxValue,new TGLayoutHints( kLHintsRight,5,5,2,2));
  fMomentumFrame->Resize(700,100);
  
  
  //Rapidity Slider
  fRapidityFrame = new TGCompositeFrame(fMainFrame,0,0, kHorizontalFrame);
  fRapidityLayout = new TGLayoutHints( kLHintsLeft | kLHintsTop | kLHintsExpandX,5,5,2,2);
  fRapiditySlider = new TGDoubleHSlider(fRapidityFrame, 400, kDoubleScaleBoth, kIdsRAPIDITY);
  fRapiditySlider->Connect("PositionChanged()","AliSliderFrame",this,"DoSlider()");
  fRapiditySlider->Connect("PositionChanged()","AliSliderFrame",this,"DoPositionChanged(Int_t)");
  fRapiditySlider->Connect("Released()","AliSliderFrame",this,"DoReleased()");
  SetRapidityRange(-1.5,1.5);
  fRapiditySlider->SetPosition(-1.5,1.5);
  fRapidityLabel = new TGLabel(fRapidityFrame,"Rapidity    ");
  fRapidityMinValue = new TGNumberEntry(fRapidityFrame,GetRapidityMin(),7,kIdsRapidityMIN);
  fRapidityMinValue->Connect("ValueChanged(Long_t)","AliSliderFrame",this,"DoField(Long_t)");
  fRapidityMinValue->SetButtonToNum(kFALSE);
  fRapidityMinValue->GetNumberEntry()->Connect("ReturnPressed()","AliSliderFrame",this,"DoField(Long_t)");
  
  fRapidityMaxValue = new TGNumberEntry(fRapidityFrame,GetRapidityMax(),7,kIdsRapidityMAX);
  fRapidityMaxValue->Connect("ValueChanged(Long_t)","AliSliderFrame",this,"DoField(Long_t)");
  fRapidityMaxValue->SetButtonToNum(kFALSE);
  fRapidityMaxValue->GetNumberEntry()->Connect("ReturnPressed()","AliSliderFrame",this,"DoField(Long_t)");
  
  fRapidityFrame->AddFrame(fRapidityLabel,new TGLayoutHints( kLHintsLeft | kLHintsCenterY,5,5,2,2));
  fRapidityFrame->AddFrame(fRapidityMinValue,new TGLayoutHints( kLHintsLeft ,5,5,2,2));
  fRapidityFrame->AddFrame(fRapiditySlider,fRapidityLayout);
  fRapidityFrame->AddFrame(fRapidityMaxValue,new TGLayoutHints( kLHintsRight ,5,5,2,2));
  fRapidityFrame->Resize(700,100);
  
  fMainFrame->AddFrame(fMomentumFrame,fLayout);
  fMainFrame->AddFrame(fRapidityFrame,fLayout);
  fMainFrame->MapSubwindows();
  fMainFrame->MapWindow();
  fMainFrame->Resize(700,100);
  LoadFromRC();
}

//_____________________________________________________________
AliSliderFrame::~AliSliderFrame(void)
{
  // Destructor
  delete fLayout;
  delete fRapidityLayout;
  delete fMomentumLayout;
  delete fMomentumMinValue;
  delete fMomentumMaxValue;
  delete fMainFrame;
  delete fRapidityFrame;
  delete fMomentumFrame;
  delete fMomentumSlider;
  delete fRapiditySlider;
  delete fMomentumLabel;
  delete fRapidityLabel;
  delete fRapidityMinValue;
  delete fRapidityMaxValue;
}

//_____________________________________________________________
void AliSliderFrame::CloseWindow(void)
{
  // Deletes this window
  delete this;
}

//_____________________________________________________________
void AliSliderFrame::DoSlider(Int_t /*pos*/)
{
  // Updates the values in case one moves a slider
  TGFrame *frame = (TGFrame *) gTQSender;
  TGDoubleSlider * ds = (TGDoubleSlider *) frame;
  int id = ds->WidgetId();
  char min[8];
  char max[8];
  switch(id){
  case kIdsMOMENTUM:{
    //sprintf(buf,"momentum min:%f max:%f",GetMomentumMin(),GetMomentumMax());
    //printf("\n%s",buf);
    sprintf(min,"%.4f",GetMomentumMin());
    sprintf(max,"%.4f",GetMomentumMax());
    fMomentumMaxValue->SetText(max);
    fMomentumMinValue->SetText(min);
  }
    break;
  case kIdsRAPIDITY:{
    sprintf(min,"%.4f",GetRapidityMin());
    sprintf(max,"%.4f",GetRapidityMax());
    fRapidityMaxValue->SetText(max);
    fRapidityMinValue->SetText(min);
  }
    break;
  default:break;
  }
}

//_____________________________________________________________
void AliSliderFrame::DoReleased(Int_t /*pos*/) const
{
  // Updates the display when the slider is released
  TGFrame *frame = (TGFrame *) gTQSender;
  TGDoubleSlider * ds = (TGDoubleSlider *) frame;
  int id = ds->WidgetId();
  switch(id){
  case kIdsRAPIDITY:
  case kIdsMOMENTUM:
    gAliDisplay2->Update(kmCUTS);
    break;
  default: break;
  }
}

//_____________________________________________________________
void AliSliderFrame::DoPositionChanged(Int_t /*pos*/) const
{
  // Updates the display when a slider has changed
  if(!gAliDisplay2->GetSliderUpdate()) return;
  TGFrame *frame = (TGFrame *) gTQSender;
  TGDoubleSlider * ds = (TGDoubleSlider *) frame;
  int id = ds->WidgetId();
  switch(id){
  case kIdsRapidityMIN:
  case kIdsRapidityMAX:
  case kIdsMomentumMIN:
  case kIdsMomentumMAX:{
    gAliDisplay2->Update(kmCUTS);
  }
    break;
  default: break;
  }
}

//_____________________________________________________________
void AliSliderFrame::DoField(Long_t pos)
{
  // Updates the display in case of changed min/max values for
  // momentum and/or rapidity
  TGNumberEntry *ne = (TGNumberEntry *) gTQSender;
  int id = ne->WidgetId();
  char max[8],min[8];
  Int_t sign=0;
  Float_t step=gAliDisplay2->GetSliderStep();
  if((pos/10000)==0){//Up button pressed
    sign = 1;
  }
  else sign = -1;
  
  switch(id){
  case kIdsMomentumMIN:{
    fMomentumMinValue->SetNumber(fMomentumMinValue->GetNumber()+step*sign);
    sprintf(min,"%.4f",fMomentumMinValue->GetNumber());
    fMomentumMinValue->SetText(min);
  }
    break;
  case kIdsMomentumMAX:{
    fMomentumMaxValue->SetNumber(fMomentumMaxValue->GetNumber()+step*sign);
    sprintf(max,"%.4f",fMomentumMaxValue->GetNumber());
    fMomentumMaxValue->SetText(max);
  }
    break;
  case kIdsRapidityMIN:{
    fRapidityMinValue->SetNumber(fRapidityMinValue->GetNumber()+step*sign);
    sprintf(min,"%.4f",fRapidityMinValue->GetNumber());
    fRapidityMinValue->SetText(min);
  }
    break;
  case kIdsRapidityMAX:{
    fRapidityMaxValue->SetNumber(fRapidityMaxValue->GetNumber()+step*sign);
    sprintf(max,"%.4f",fRapidityMaxValue->GetNumber());
    fRapidityMaxValue->SetText(max);
  }
    break;
  default:break;
  }
  
  fMomentumSlider->SetPosition(fMomentumMinValue->GetNumber(),fMomentumMaxValue->GetNumber());
  fRapiditySlider->SetPosition(fRapidityMinValue->GetNumber(),fRapidityMaxValue->GetNumber());
  //	gAliDisplay2->Draw();
  gAliDisplay2->Update(kmCUTS);
}

//_____________________________________________________________
void AliSliderFrame::SaveToRC() const
{
  // Saves settings in the .alidisplayrc file
  TEnv *rc=new TEnv(".alidisplayrc");
  rc->SetValue("AliDisplay.MomentumMin",GetMomentumMin());
  rc->SetValue("AliDisplay.MomentumMax",GetMomentumMax());
  rc->SetValue("AliDisplay.RapidityMin",GetRapidityMin());
  rc->SetValue("AliDisplay.RapidityMax",GetRapidityMax());
  rc->SaveLevel(kEnvLocal);
  rc->Save();
}

//_____________________________________________________________	
void AliSliderFrame::LoadFromRC()
{
  // Loads settings from the .alidisplayrc file
  // and uses the default momentum and rapidity ranges
  TEnv *rc=new TEnv(".alidisplayrc");
  Float_t a,b;
  a=rc->GetValue("AliDisplay.MomentumMin",0);
  b=rc->GetValue("AliDisplay.MomentumMax",2);
  fMomentumSlider->SetPosition(a,b);
  a=rc->GetValue("AliDisplay.RapidityMin",-1.5);
  b=rc->GetValue("AliDisplay.RapidityMax",1.5);
  fRapiditySlider->SetPosition(a,b);
}


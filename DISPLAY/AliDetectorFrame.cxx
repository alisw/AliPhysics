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
// ALICE DETECTOR FRAME CLASS                                          //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <TGWindow.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGLayout.h>
#include <TObjArray.h>

#include "AliDetectorFrame.h"
#include "AliDisplay2.h"
#include "AliModuleInfo.h"

#include "AliModule.h"

ClassImp(AliDetectorFrame)

int AliDetectorFrame::fgBaseId = 1000;

//_____________________________________________________________
AliDetectorFrame::AliDetectorFrame(const TGWindow *p, Int_t w, Int_t h,UInt_t bgc)
{
  // Constructor
  fMainFrame = new TGCompositeFrame(p,w,h,kVerticalFrame,bgc);
  TGLayoutHints *layout = new TGLayoutHints(kLHintsTop | kLHintsLeft,2,2,2,2);
  TGLayoutHints *layout2 = new TGLayoutHints(kLHintsTop | kLHintsRight,2,2,2,2);
  TGLayoutHints *layout3 = new TGLayoutHints(kLHintsTop | kLHintsExpandX,2,2,2,2);
  fCheckButton = new TGCheckButton*[gAliDisplay2->GetNbModules()];
  fCheckedButton = new Bool_t[gAliDisplay2->GetNbModules()];
  fCheckButtonId = new Int_t[gAliDisplay2->GetNbModules()];
  TGCompositeFrame *dframe;
  TGButton 		 *button;
  char			 text[32];
  AliModule 	*mod;
  fCheckedMode = kFALSE;
  for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
    mod = dynamic_cast<AliModule*> (gAliDisplay2->GetModules()->At(i));
    if(!mod) continue;
    dframe = new TGCompositeFrame(fMainFrame,150,20,kHorizontalFrame);
    fCheckButton[i] = new TGCheckButton(dframe,mod->GetName(),fgBaseId);
    fCheckButtonId[i]=fgBaseId;
    fCheckButton[i]->Connect("Clicked()","AliDetectorFrame",this,"DoCheckButton(Int_t)");		
    fCheckedButton[i]=kTRUE;
    dframe->AddFrame(fCheckButton[i],layout);
    fCheckButton[i]->SetState(kButtonDown);
    sprintf(text,"Specific %s view",mod->GetName());
    button = new TGTextButton(dframe,"Display",fgBaseId);
    button->SetToolTipText(text);
    button->Connect("Clicked()","AliDetectorFrame",this,"DoSpecific()");
    dframe->AddFrame(button,layout2);
    fMainFrame->AddFrame(dframe,layout3);
    gAliDisplay2->GetModuleInfo()->SetId((char*)mod->GetName(),fgBaseId);
    fgBaseId++;
  }	
  gAliDisplay2->GetModuleInfo()->Print();
  fButtonFrame = new TGCompositeFrame(fMainFrame,w,100,kHorizontalFrame,bgc);
  fButtonAll = new TGTextButton(fButtonFrame,"All",kIdbSelectALL);
  fButtonAll->Connect("Clicked()","AliDetectorFrame",this,"DoButton(Int_t)");
  fButtonFrame->AddFrame(fButtonAll,new TGLayoutHints(kLHintsBottom | kLHintsLeft,2,2,2,2));
  fButtonInvert = new TGTextButton(fButtonFrame,"Invert",kIdbSelectINVERT);
  fButtonInvert->Connect("Clicked()","AliDetectorFrame",this,"DoButton(Int_t)");
  fButtonFrame->AddFrame(fButtonInvert,new TGLayoutHints(kLHintsBottom | kLHintsRight,2,2,2,2));
  fMainFrame->AddFrame(fButtonFrame,new TGLayoutHints(kLHintsBottom | kLHintsLeft|kLHintsExpandX,0,0,2,2));		
  fCheckedMode = kTRUE;
}

//_____________________________________________________________
AliDetectorFrame::~AliDetectorFrame()
{
  // Destructor
  for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
    delete fCheckButton[i];
  }
  delete fCheckButton;
  delete fCheckedButton;
  delete fCheckButtonId;
  //delete [] fDetectorName;
  delete fButtonFrame;
  delete fMainFrame;
  delete fButtonAll;
  delete fButtonInvert;
}

//_____________________________________________________________
void AliDetectorFrame::DoButton(Int_t /*pos*/)
{
  // Update display if a button was used
  TGFrame *frame = (TGFrame *) gTQSender;
  TGButton *bu= (TGButton *) frame;
  int id = bu->WidgetId();	
  fCheckedMode = kFALSE;
  AliModule *mo;
  switch(id){
  case kIdbSelectALL:{
    for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
      mo = dynamic_cast<AliModule *> (gAliDisplay2->GetModules()->At(i));
      if(!mo) continue;
      fCheckButton[i]->SetState(kButtonDown);
      fCheckedButton[i]=kTRUE;
      gAliDisplay2->EnableDetector(mo->GetName());
    }
  }
    break;
  case kIdbSelectINVERT:{
    for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
      mo = dynamic_cast<AliModule *> (gAliDisplay2->GetModules()->At(i));
      if(!mo) continue;
      if(fCheckedButton[i]==kTRUE) {
	fCheckButton[i]->SetState(kButtonUp);
	fCheckedButton[i]=kFALSE;
	gAliDisplay2->DisableDetector(mo->GetName());
      }
      else if(fCheckedButton[i]==kFALSE)  {
	fCheckButton[i]->SetState(kButtonDown);
	fCheckedButton[i]=kTRUE;
	gAliDisplay2->EnableDetector(mo->GetName());
      }
    }
  }
    break;
  default:break;
  }
  gAliDisplay2->Update(kmMODULES);
  fCheckedMode = kTRUE;
}

//_____________________________________________________________
void AliDetectorFrame::DoCheckButton(Int_t /*pos*/)
{
  // Chech if any button was used
  if(fCheckedMode == kFALSE) return;
  TGFrame *frame = (TGFrame *) gTQSender;
  TGCheckButton *bu= (TGCheckButton *) frame;
  Int_t id = bu->WidgetId();
  AliModule *mo;
  for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
    mo = dynamic_cast<AliModule *> (gAliDisplay2->GetModules()->At(i));
    if(!mo) continue;
    if(id==fCheckButtonId[i]){
      if(fCheckedButton[i]==kTRUE) {
	fCheckedButton[i]=kFALSE;
	gAliDisplay2->DisableDetector(mo->GetName());
      }
      else {
	fCheckedButton[i]=kTRUE;
	gAliDisplay2->EnableDetector(mo->GetName());
      }
    }
  }
  gAliDisplay2->Update(kmMODULES);
}

//_____________________________________________________________
void AliDetectorFrame::DoSpecific() const
{
  // Draw detectors
  TGFrame *frame = (TGFrame *) gTQSender;
  TGCheckButton *bu= (TGCheckButton *) frame;
  Int_t id = bu->WidgetId();
  AliModule *mo;
  for(Int_t i=0;i<gAliDisplay2->GetNbModules();i++){
    mo = dynamic_cast<AliModule *> (gAliDisplay2->GetModules()->At(i));
    if(!mo) continue;
    if(id==fCheckButtonId[i]){
      gAliDisplay2->DrawDetector(mo->GetName());
    }
  }
}

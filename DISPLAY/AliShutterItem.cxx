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
// ALICE SHUTTER ITEM CLASS                                            //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <TGShutter.h>
#include <TGFrame.h>
#include <TGButton.h>

#include "AliShutterItem.h"
#include "AliDisplay2.h"

ClassImp(AliShutterItem)

//_____________________________________________________________
AliShutterItem::AliShutterItem(TGShutter *s, const char *text, UInt_t id)
{
  // Constructor
  fShutterItem = new TGShutterItem(s, new TGHotString(text), id);
  fMainFrame = (TGCompositeFrame *) fShutterItem->GetContainer();
  s->AddItem(fShutterItem);
}

//_____________________________________________________________
AliShutterItem::~AliShutterItem(void)
{
  // Destructor
  delete fButton;
  delete fShutterItem;
  delete fMainFrame;
}

//_____________________________________________________________
void AliShutterItem::AddTextButton(const char* text, const char *tiptext, UInt_t idb)
{
  //Add a TGTextButton in the TGShutterItem. This button will execute the fonction
  fButton = new TGTextButton(fMainFrame,new TGHotString(text),idb);
  fButton->Resize(100,fButton->GetDefaultHeight());
  fButton->Connect("Clicked()","AliShutterItem",this,"DoButton(Int_t)");
  fButton->SetToolTipText(tiptext);
  //fButton->Connect("Clicked()","AliDisplay2",gAliDisplay2,"DoViews(Int_t)");
  fMainFrame->AddFrame(fButton, new TGLayoutHints( kLHintsTop | kLHintsCenterX ,5,5,10,10));
}

//_____________________________________________________________
void AliShutterItem::AddPictureButton(const char* file, const char *tiptext, UInt_t idb)
{
  //Add a TGPictureButton in the TGShutterItem. The icon file must be in DISPLAY/icons
  TString filename=StrDup(gAliDisplay2->GetIconsPath());
  filename.Append(file);
  TGPicture *picture = (TGPicture *) gClient->GetPicture(filename);
  fButton = new TGPictureButton(fMainFrame,picture,idb);		
  fButton->SetToolTipText(tiptext);
  fButton->Connect("Clicked()","AliShutterItem",this,"DoButton(Int_t)");
  fMainFrame->AddFrame(fButton, new TGLayoutHints( kLHintsTop | kLHintsCenterX ,5,5,10,10));
}

//_____________________________________________________________
void AliShutterItem::AddCheckButton(const char *text,Int_t idb)
{
  // Add check button
  fButton = new TGCheckButton(fMainFrame,new TGHotString(text),idb);
  fButton->Resize(100,fButton->GetDefaultHeight());
  fButton->Connect("Clicked()","AliShutterItem",this,"DoButton(Int_t)");
  fMainFrame->AddFrame(fButton, new TGLayoutHints( kLHintsTop | kLHintsLeft ,5,5,10,10));
}

//_____________________________________________________________
void AliShutterItem::DoButton(Int_t /*pos*/) const
{
  // Show next/previous event if the buttom was used
  TGFrame *frame = (TGFrame *) gTQSender;
  TGButton *bu= (TGButton *) frame;
  int id = bu->WidgetId();
  switch(id){
  case kIdbNextEVENT:{
    gAliDisplay2->ShowNextEvent(1);
  }
    break;
  case kIdbPrevEVENT:{
    gAliDisplay2->ShowNextEvent(-1);
  }
    break;
  case kIdbCheckHITS:{
    if(gAliDisplay2->IsEnabled(kHits)) gAliDisplay2->Disable(kHits);
    else gAliDisplay2->Enable(kHits);	  
  }
    break;
  case kIdbCheckCLUSTERS:{
    if(gAliDisplay2->IsEnabled(kClusters)) gAliDisplay2->Disable(kClusters);
    else gAliDisplay2->Enable(kClusters);
  }
    break;
  case kIdbCheckHLT:{
    if(gAliDisplay2->IsEnabled(kHLT)) gAliDisplay2->Disable(kHLT);
    else gAliDisplay2->Enable(kHLT);
  }
    break;
  case kIdbCheckTRACKS:{
    if(gAliDisplay2->IsEnabled(kTracks)) gAliDisplay2->Disable(kTracks);
    else gAliDisplay2->Enable(kTracks);
  }
    break;
  case kIdbSIDEVIEW:{
    gAliDisplay2->DoView(kIdbSIDEVIEW);
  }
    break;
  case kIdbFRONTVIEW:{
    gAliDisplay2->DoView(kIdbFRONTVIEW);
  }
    break;
  case kIdbTOPVIEW:{
    gAliDisplay2->DoView(kIdbTOPVIEW);
  }
    break;
  case kIdbALLVIEW:{
    gAliDisplay2->DoView(kIdbALLVIEW);
  }
    break;
  default:break;
  }
}


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
// ALICE EVENT DISPLAY CLASS                                           //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#define do_mc

//standard modules
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

//ROOT

#include <TTree.h>
#include <TGLayout.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TGWindow.h>
#include <TEnv.h>
#include <TPad.h>

//AliRoot Module
#include "AliModule.h"
#include "AliDetector.h"
#include "AliRun.h"
#include "AliMC.h"

#include "AliModuleInfo.h"
#include "AliSliderFrame.h"
#include "AliShutterFrame.h"
#include "AliDisplayFrame.h"
#include "AliInfoFrame.h"
#include "AliDetectorFrame.h"
#include "AliMenu.h"

#include "AliDisplay2.h"


 AliDisplay2 *gAliDisplay2;


ClassImp(AliDisplay2)
//_____________________________________________________________
AliDisplay2::AliDisplay2(const TGWindow *p, UInt_t w, UInt_t h)
			:TObject()
{
  //  Constructor
  //gAlice->SetDisplay(this);
  gAliDisplay2=this;
  fSliderUpdate = kFALSE;
  fZoomMode = kFALSE;
  fZoomStep = 1.2;
  fZoomFactor = 1.5;
  fNbParticles = 0;
  fEventNumber = 0;
  fNbHits = 0;
  fSliderStep = 0.01;
  fClustersLoaded = kFALSE;
  fHitsLoaded = kFALSE;
  fHLTLoaded = kFALSE;
  fTracksLoaded = kFALSE;
  fMode =0;
  FindModules();
  
  fIconsPath = new char[32];
  strcpy(fIconsPath,gSystem->Getenv("ALICE_ROOT"));
  strcat(fIconsPath,"/DISPLAY/icons/");
  LoadFromRC();
  fMainFrame = new TGMainFrame(p,w,h,kVerticalFrame);
  fSubFrame = new TGCompositeFrame(fMainFrame,w,h,kHorizontalFrame);
  fLeftFrame = new TGCompositeFrame(fSubFrame,150,h,kVerticalFrame|kFixedWidth);
  fRightFrame = new TGCompositeFrame(fSubFrame,600,h,kVerticalFrame);
  //fMainFrame->Connect("ProcessedEvent(Event_t*)", "AliDisplay2", this,"HandleMouseWheel(Event_t*)");
  fMainFrame->Connect("ProcessedEvent(Event_t*)", "AliDisplay2",this,"HandleResize(Event_t*)");
  //MenuBar
  fMenu = new AliMenu(fMainFrame,1,1,kRaisedFrame|kHorizontalFrame);
  
  //Slider Frame
  fSliderFrameLayout = new TGLayoutHints( kLHintsBottom| kLHintsRight| kLHintsExpandX | kLHintsCenterX, 2, 2, 2, 2);
  fSliderFrame  = new AliSliderFrame(fRightFrame,600,150);
  fRightFrame->AddFrame(fSliderFrame->GetSliderFrame(),fSliderFrameLayout);
  
  //Info Frame
  fInfoFrameLayout = new TGLayoutHints( kLHintsTop | kLHintsLeft  | kLHintsExpandX  ,0,0,0,0);
  fInfoFrame = new AliInfoFrame(fLeftFrame,150,200);
  fLeftFrame->AddFrame(fInfoFrame->GetInfoFrame(),fInfoFrameLayout);
  
  
  //Shutter Frame
  fShutterFrameLayout = new TGLayoutHints( kLHintsTop | kLHintsLeft | kLHintsExpandY |kLHintsExpandX, 0, 0, 5, 0);
  fShutterFrame = new AliShutterFrame(fLeftFrame,150,300);
  fLeftFrame->AddFrame(fShutterFrame->GetShutterFrame(),fShutterFrameLayout);
  
  //Display Frame
  fDisplayFrameLayout = new TGLayoutHints( kLHintsTop | kLHintsRight | kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2); 
  fDisplayFrame = new AliDisplayFrame(fRightFrame, w-150,w-110);
  fRightFrame->AddFrame(fDisplayFrame->GetDisplayFrame(),fDisplayFrameLayout);
  fDisplayFrame->GetDisplayFrame()->Connect("ProcessedEvent(Event_t*)", "AliDisplay2", this,"HandleMouseWheel(Event_t*)");
  
  
  fLeftFrame->Layout();
  
  fSubFrame->AddFrame(fLeftFrame, new TGLayoutHints( kLHintsBottom | kLHintsLeft | kLHintsExpandY, 5, 5, 2, 2));
  fSubFrame->AddFrame(fRightFrame, new TGLayoutHints( kLHintsBottom | kLHintsRight | kLHintsExpandX | kLHintsExpandY, 5, 5, 2, 2));
  
  Int_t parts[] = {45,45,10};
  fStatusBar = new TGStatusBar(fMainFrame,50,10,kHorizontalFrame);
  fStatusBar->SetParts(parts,3);
  fStatusBar->SetText("AliDisplay v2.0",0);
  fMainFrame->AddFrame(fStatusBar,new TGLayoutHints(kLHintsBottom | kLHintsExpandX,0,0,0,0));
  
  fMainFrame->AddFrame(fSubFrame,new TGLayoutHints( kLHintsBottom | kLHintsLeft | kLHintsExpandX | kLHintsExpandY, 5, 5, 2, 2));
  fMainFrame->SetWindowName("Ali Display");
  
  fMainFrame->MapSubwindows();
  fMainFrame->MapWindow();
  LoadSettings();

	
  fMainFrame->Resize(w-10,h);
  fMainFrame->Resize(w,h);
  fMainFrame->SetWMSizeHints(500,500,1280, 1200,1,1);
}

//_____________________________________________________________
AliDisplay2::~AliDisplay2(void)
{
  //Destructor
  delete fModules;
  delete [] fEnabledModules;
  delete fModuleInfo;
  
  delete fSliderFrameLayout;
  delete fSliderFrame;	
  delete fDisplayFrameLayout;
  delete fDisplayFrame;
  //	delete fZoomFrameLayout;
  //	delete fZoomFrame;
  delete fShutterFrameLayout;
  delete fShutterFrame;
  delete fInfoFrameLayout;
  delete fInfoFrame;
  delete fDetectorFrameLayout;
  delete fDetectorFrame;
  
  delete fSubFrame;
  delete fLeftFrame;
  delete fRightFrame;
  delete fMainFrame;
  delete fAliDisplay2rc;
  
  delete fMenu;
  delete fStatusBar;
}

//_____________________________________________________________
void AliDisplay2::CloseWindow(void)
{
  // Deletes the current display
  delete this;
}

//_____________________________________________________________
void AliDisplay2::LoadFromRC()
{
  // Load the environment settings from .alidisplayrc file
  TEnv *rc=new TEnv(".alidisplayrc");
  SetSliderUpdate(rc->GetValue("AliDisplay.SliderUpdate",kFALSE));
  SetZoomStep(rc->GetValue("AliDisplay.ZoomStep",1.2));
  SetSliderStep(rc->GetValue("AliDisplay.SliderStep",0.01));
  char c[128];
  fRawDataPath = new char[128];
  strcpy(c,gSystem->Getenv("ALICE_ROOT"));
  sprintf(fRawDataPath,"%s%s",c,rc->GetValue("AliDisplay.RawDataPath","/raw"));
  printf("\nRaw data path %s",fRawDataPath);
}

//_____________________________________________________________
void AliDisplay2::SaveToRC() const
{
  // Saves the environment settings in .alidisplayrc file
  TEnv *rc=new TEnv(".alidisplayrc");
  rc->SetValue("AliDisplay.SliderUpdate",GetSliderUpdate());
  rc->SetValue("AliDisplay.ZoomStep",GetZoomStep());
  rc->SetValue("AliDisplay.SliderStep",GetSliderStep());
  rc->SetValue("AliDisplay.RawDataPath","/raw");
  rc->SaveLevel(kEnvLocal);
  rc->Save();
}

//_____________________________________________________________
void AliDisplay2::DoSaveSettings(void)
{
  // Saves the environment settings for the slider frame and display
  fSliderFrame->SaveToRC();
  SaveToRC();
}

//_____________________________________________________________
void AliDisplay2::LoadSettings()
{
  // Loads settings
  LoadFromRC();
}

//_____________________________________________________________
void AliDisplay2::Draw(Option_t */*options*/)
{
  // Draws display frame
  fDisplayFrame->DoView(fCurrentView);
}

//_____________________________________________________________
void AliDisplay2::DrawX3d()
{
  // Draws display frame using X3d
  fDisplayFrame->DrawX3d();
}

//_____________________________________________________________
void AliDisplay2::DrawGL()
{
  // Draws display frame using GL
  fDisplayFrame->DrawGL();
}

//_____________________________________________________________
void AliDisplay2::ShowNextEvent(Int_t delta)
{
  //Load the next event
  clock_t t1,t2;
  t1=clock();
  Int_t newEvent=0;
  if(delta!=0){
    gAlice->Clear();
    newEvent = fEventNumber + delta;
    if( newEvent < 0) return;
    gAlice->GetEvent(newEvent);
    fEventNumber += delta;
    //		if(!gAlice->TreeH()) return;
  }
  if(IsEnabled(kHits)) LoadHits();
  if(IsEnabled(kClusters)) LoadClusters(newEvent);
  if(IsEnabled(kHLT)) LoadHLTClusters(newEvent);
  t2=clock();
  //	printf("\nEvent loaded in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
  Update(kmMODULES);
}

//_____________________________________________________________
void AliDisplay2::FindModules()
{
  // Find the modules used for the simulation and assign 
  // these modules to the array fModules
  fModules = new TObjArray;
  TObjArray *modules = gAlice->Modules();
  AliModule *mod;
  Int_t nbm = 0;
  for(Int_t i=0;i<modules->GetEntriesFast();i++){
    mod = (AliModule *) modules->At(i);
    if(!mod) continue;
    const char *avoid = strstr("BODY MAG ABSO DIPO HALL FRAME SHIL PIPE",mod->GetName());
    if(avoid) continue;
    fModules->AddLast(mod);
    nbm++;
  }
  fEnabledModules = new Bool_t[nbm];
  fNbModules = nbm;
  fModuleInfo = new AliModuleInfo(nbm);
  for(Int_t j=0;j<fModules->GetEntriesFast();j++){
    fModuleInfo->Add(fModules->At(j)->GetName(),j);
    fEnabledModules[j]=kTRUE;
  }
}

//_____________________________________________________________
void AliDisplay2::LoadHits()
{
  //Load the detected hits from each detector to memory
  gAlice->ResetPoints();
  TIter next(gAlice->Modules());
  AliModule *module;
  Int_t ntracks = gAlice->GetMCApp()->GetNtrack();
  while((module = (AliModule*)next())) 
    {
      AliDetector* detector = dynamic_cast<AliDetector*>(module);
      if(detector) detector->SetTreeAddress();
    }
  next.Reset();
  for (Int_t track=0; track<ntracks;track++) {
    gAlice->ResetHits();
    while((module = (AliModule*)next())) {
      AliDetector* detector = dynamic_cast<AliDetector*>(module);
      if(detector)
	{
	  detector->TreeH()->GetEvent(track);
	  detector->LoadPoints(track);
	}
    }
    next.Reset();
  }
  fHitsLoaded = kTRUE;
}

//_____________________________________________________________
void AliDisplay2::LoadClusters(Int_t nevent)
{
  //clock_t t1,t2;
  fDisplayFrame->LoadClusters(nevent);
  fClustersLoaded = kTRUE;
  //	printf("\nClusters loaded in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
}

//_____________________________________________________________
void AliDisplay2::LoadHLTClusters(Int_t nevent)
{
  // Loads HLT clusters
  fDisplayFrame->LoadHLTClusters(nevent);
  fHLTLoaded = kTRUE;
}

//_____________________________________________________________
void AliDisplay2::Enable(Int_t m)
{
  // Enables the given mode m
  if(m==kHits){
    if((fMode&kHits)==kHits) return;
    fMode = kHits|fMode;
    if(!fHitsLoaded) LoadHits();
    Update(kmPOINTS);
  }
  if(m==kClusters){
    if((fMode&kClusters)==kClusters) return;
    fMode = kClusters|fMode;
    if(!fClustersLoaded) LoadClusters(fEventNumber);
    Update();
  }
  if(m==kHLT){
    if((fMode&kHLT)==kHLT) return;
    fMode = kHLT|fMode;
    if(!fHLTLoaded) {
      LoadHLTClusters(fEventNumber);
    }
    Update();
  }
  if(m==kTracks){
    if((fMode&kTracks)==kTracks) return;
    fMode = kTracks|fMode;
    Update();
  }
}

//_____________________________________________________________
void AliDisplay2::Disable(Int_t m)
{
  // Disables the given mode m
  if(m==kHits){
    fMode = fMode|kHits;
    fMode = fMode^kHits;
  }
  if(m==kClusters){
    fMode = fMode|kClusters;
    fMode = fMode^kClusters;
  }
  if(m==kHLT){
    fMode = fMode|kHLT;
    fMode = fMode^kHLT;
  }
  if(m==kTracks){
    fMode = fMode|kTracks;
    fMode = fMode^kTracks;
  }
  Update();
}

//_____________________________________________________________
Bool_t AliDisplay2::IsEnabled(Int_t m) const
{
  // Checks if the mode m is enabled
  if(m==kHits){
    if((fMode&kHits)==kHits) return kTRUE;
    return kFALSE;
  }
  if(m==kClusters){
    if((fMode&kClusters)==kClusters) return kTRUE;
    return kFALSE;
  }
  if(m==kHLT){
    if((fMode&kHLT)==kHLT) return kTRUE;
    return kFALSE;
  }
  if(m==kTracks){
    if((fMode&kTracks)==kTracks) return kTRUE;
    return kFALSE;
  }
  return kFALSE;
}

//_____________________________________________________________
void AliDisplay2::HandleMouseWheel(Event_t *event)
{
  //Handle mouve event, not working yet
  if(event->fType != kButtonPress && event->fType != kButtonRelease) return;
  
  if(event->fCode == kButton4){
    fZoomFactor *= fZoomStep;
    Draw();
  }
  
  if(event->fCode == kButton5){
    fZoomFactor /= fZoomStep;
    Draw();
  }
}

//_____________________________________________________________
void AliDisplay2::HandleResize(Event_t *event)
{
  // Handle resize event
  switch(event->fType){
  case kConfigureNotify:{
    Draw();
  }
    break;
  default:break;
  }
}	

//_____________________________________________________________
void AliDisplay2::Update(Int_t tag)
{
  // Update the view, if loading only the modified data from the previous
  // changes, the integer tag indicates the kind of modification
  if(tag==kmMODULES){
    LoadEnabledModules();
    if(((fMode)&kHits)==kHits){
      LoadEnabledHits();
      ApplyCuts();
    }
  }
  if(tag==kmCUTS){
    if(((fMode)&kHits)==kHits)ApplyCuts();
  }
  if(tag==kmPOINTS){
    if(((fMode)&kHits)==kHits){
      LoadEnabledHits();
      ApplyCuts();
    }
  }
  Draw();
  fInfoFrame->Update();
}



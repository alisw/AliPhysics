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
// ALICE DISPLAY FRAME CLASS                                           //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <time.h>

#include <TGFrame.h>
#include <TGTab.h>
#include <TObjArray.h>
#include <TStopwatch.h>
#include <TPad.h>
#include <TVirtualX.h>
#include <TCanvas.h>
#include <TView.h>
#include <TParticle.h>
#include <TGeometry.h>

#include "AliModuleInfo.h"
#include "AliDisplayHLT.h"
#include "AliDisplay2.h"
#include "AliModule.h"
#include "AliDetector.h"
#include "AliPoints.h"
#include "AliRun.h"


ClassImp(AliDisplayFrame);

//_____________________________________________________________
AliDisplayFrame::AliDisplayFrame(const TGWindow *p, UInt_t w, UInt_t h)
{
  // Constructor
  fClipMin=-20;
  fClipMax=20;
  fPreviousW=0;
  fPreviousH=0;
  fRange = 500;
  fPolyMarkers = new TObjArray(1000);
  
  fMainFrame = new TGCompositeFrame(p,w,h);
  fMainTab = new TGTab(fMainFrame, w, h);
  fFrame1 = fMainTab->AddTab("Main View");
  fMainEmbeddedCanvas = new TRootEmbeddedCanvas("Main12",fFrame1,w,h,kFixedWidth);
  fFrame1->AddFrame(fMainEmbeddedCanvas,new TGLayoutHints( kLHintsTop | kLHintsLeft|kLHintsExpandX| kLHintsExpandY, 0, 0, 0, 0));
  fMainCanvas = fMainEmbeddedCanvas->GetCanvas();
  fMainCanvas->SetFillColor(1);
  fMainCanvas->SetBorderMode(0);
  fMainCanvas->cd();
  fMainCanvas->SetFixedAspectRatio();
  fMainCanvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","AliDisplayFrame",this,"ExecuteEvent(Int_t,Int_t,Int_t,TObject*)");
  //fView = new TView(1);
  //DoView(kIdbFRONTVIEW);
  
  gAliDisplay2->SetCurrentView(kIdbFRONTVIEW);	
  
  fFrame2 = fMainTab->AddTab("No detector");
  fSelectionEmbeddedCanvas = new TRootEmbeddedCanvas("Selection",fFrame2,w,h);
  fSelectionCanvas = fSelectionEmbeddedCanvas->GetCanvas();
  fSelectionCanvas->SetFillColor(1);
  fSelectionCanvas->SetBorderMode(0);
  fSelectionCanvas->cd();
  fFrame2->AddFrame(fSelectionEmbeddedCanvas,new TGLayoutHints( kLHintsTop | kLHintsLeft|kLHintsExpandX| kLHintsExpandY, 0, 0, 0, 0));
  fMainFrame->AddFrame(fMainTab,new TGLayoutHints( kLHintsTop | kLHintsLeft|kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  fAllViews = kFALSE;
  fMainFrame->MapSubwindows();
  fMainFrame->MapWindow();
}

//_____________________________________________________________
AliDisplayFrame::~AliDisplayFrame(void)
{
  // Destructor
  delete fMainTab;
  delete fSelectionEmbeddedCanvas;
  delete fMainEmbeddedCanvas;
  delete fFrame1;
  delete fFrame2;
  delete fMainCanvas;
  delete fSelectionCanvas;
  delete fPoints2;
  delete fPoints;
  delete fModules;
  delete fMainFrame;
  delete [] fActivePoints;
  delete [] fClustersPos;
}

//_____________________________________________________________
void AliDisplayFrame::DoView(Int_t view)
{
  // Draws selected view
  Int_t x,y;
  char vname[16];
  y=fMainFrame->GetDefaultHeight();
  x=fMainFrame->GetDefaultWidth();
  gAliDisplay2->SetCurrentView(view);
  switch(view){
  case kIdbALLVIEW:{
    fAllViews=kTRUE;
    strcpy(vname,"All views");	
    fMainCanvas->cd();
    gPad->Clear();
    fMainCanvas->SetFillColor(15);
    fMainCanvas->Divide(2,2,0.005,0.005,1);
    
    fMainCanvas->cd(1);
    Draw(30,30,0);
    
    gAliDisplay2->SetCurrentView(kIdbTOPVIEW);
    fMainCanvas->cd(2);
    Draw(90,-90,90);
    
    gAliDisplay2->SetCurrentView(kIdbSIDEVIEW);
    fMainCanvas->cd(3);		
    Draw(90,0,-90);
    
    gAliDisplay2->SetCurrentView(kIdbFRONTVIEW);
    fMainCanvas->cd(4);
    Draw(0,-90,0);
    
    //fMainCanvas->cd();
    
  }
    break;
  case kIdbTOPVIEW:{
    strcpy(vname,"Top view  ");
    fAllViews=kFALSE;	
    fMainCanvas->cd();											
    gPad->SetFillColor(1);
    gPad->Clear();
    gPad->Draw();
    Draw(90,-90,90);
  }
    break;
  case kIdbSIDEVIEW:{
    strcpy(vname,"Side view");
    fAllViews=kFALSE;	
    fMainCanvas->cd();		
    gPad->SetFillColor(1);
    gPad->Clear();
    gPad->Draw();
    Draw(90,0,-90);
  }
    break;
  case kIdbFRONTVIEW:{
    strcpy(vname,"Front view");
    fAllViews=kFALSE;	
    fMainCanvas->cd();
    gPad->SetFillColor(1);
    gPad->Clear();
    gPad->Draw();
    
    Draw(0,-90,0);
  }
    break;
  default: break;
  }
  (fMainTab->GetTabTab(0))->SetText(new TGString(vname));
}

//_____________________________________________________________
void AliDisplayFrame::DrawDetector(const char *name)
{
  // Draws detector
  (fMainTab->GetTabTab(1))->SetText(new TGString(name));
}

//_____________________________________________________________
void AliDisplayFrame::EnableDetector(const char *name)
{
  // Enables detector
  AliModule *module = dynamic_cast<AliModule*>(gAlice->Modules()->FindObject(name));
  if(!module) return;
  gAliDisplay2->GetModuleInfo()->Enable((char*)name);
  module->Enable();
}

//_____________________________________________________________
void AliDisplayFrame::DisableDetector(const char *name)
{
  // Disables detector
  AliModule *module = dynamic_cast<AliModule*>(gAlice->Modules()->FindObject(name));
  if(!module) return;
  gAliDisplay2->GetModuleInfo()->Disable((char*)name);
  module->Disable();
}

//_____________________________________________________________
void AliDisplayFrame::Draw(Float_t theta, Float_t phi, Float_t psi)
{
  // Draws everything???
  //clock_t t1,t2;
  time_t t1,t2;
  //t1 = clock();
  TStopwatch timer;
  timer.Start();
  time(&t1);
  gPad->SetCursor(kWatch);
  gPad->SetEditable(kTRUE);
  gPad->SetFillColor(1);
  gPad->Clear();
  
  Int_t iret;
  
  TView *view = new TView(1);
  TGDimension dim=((TGCanvas*)fMainEmbeddedCanvas)->GetViewPort()->GetDefaultSize();
  Float_t aspectRatio = dim.fWidth/(Float_t) dim.fHeight;
  //printf("Dimension %d %d",dim.fWidth,dim.fHeight);
  if(gAliDisplay2->GetCurrentView()==kIdbFRONTVIEW){
    view->SetRange(-fRange*aspectRatio,-fRange,-fRange,fRange*aspectRatio,fRange,fRange);
  }
  if(gAliDisplay2->GetCurrentView()==kIdbTOPVIEW){
    view->SetRange(-fRange,-fRange,-fRange*aspectRatio,fRange,fRange,fRange*aspectRatio);
  }
  if(gAliDisplay2->GetCurrentView()==kIdbSIDEVIEW){
    view->SetRange(-fRange,-fRange,-fRange*aspectRatio,fRange,fRange,fRange*aspectRatio);
  }
  
  gAlice->GetGeometry()->Draw("same");
  if(gAliDisplay2->IsEnabled(kHits)) DrawHits();
  if(gAliDisplay2->IsEnabled(kClusters)) fClusters->Draw();
  if(gAliDisplay2->IsEnabled(kHLT)) fHLT->Draw();
  
  gAliDisplay2->AppendPad();
  view->SetView(phi,theta,psi,iret);
  
  view->ZoomView(gPad,gAliDisplay2->GetZoomFactor());
  //t2 = clock();
  time(&t2);
  //	printf("\nDrawn in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
  printf("\nDrawn in....%f sec", difftime(t2,t1));
  timer.Stop();
  timer.Print("m");
}

//_____________________________________________________________
void AliDisplayFrame::DrawHits()
{
  // Draws hits
  AliPoints *p;
  if(!fPoints2) return;
  for(Int_t i=0;i<fPoints2->GetEntries();i++){
    if(fActivePoints[i]){
      p=dynamic_cast<AliPoints *>(fPoints2->UncheckedAt(i));	
      if(!p) continue;
      p->Draw();
    }
  }
}

//_____________________________________________________________
void AliDisplayFrame::LoadEnabledModules()
{
  // Loads enabled modules
  clock_t t1,t2;
  t1=clock(); 
  TIter next(gAlice->Modules());
  AliModule *module;
  fModules = new TObjArray(0,32);
  while((module = dynamic_cast <AliModule*> (next()))){
    if(!module) continue;
    if(!module->IsActive()) continue;
    fModules->AddLast(module);
  }
  t2=clock();
  fNbModules = fModules->GetEntriesFast();
  //	printf("\nModules loaded in.....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
}

//_____________________________________________________________
void AliDisplayFrame::LoadClusters(Int_t nevent)
{
  // Loads clusters
  fClusters = new AliDisplayClusters();
  fClusters->LoadClusters("ITS TPC",nevent);
}

//_____________________________________________________________
void AliDisplayFrame::LoadHLTClusters(Int_t nevent)
{
  // Loads HLT clusters
  fHLT = new AliDisplayHLT();
  fHLT->LoadHLT("TPC",nevent);
}
	
//_____________________________________________________________
void AliDisplayFrame::LoadHits()
{
  // Loads hits
  clock_t t1,t2;

  t1=clock(); 
  fPoints2 = new TObjArray(0,1000);
  AliModule *module;
  TObjArray *points;
  for(Int_t i=0;i<fNbModules;i++){
    module = dynamic_cast<AliModule*>(fModules->UncheckedAt(i));
    if(!module) continue;
    points = module->Points();
    if(!points) {
      continue;
    }
    for(Int_t j=0;j<points->GetEntriesFast();j++){
      if(!points->UncheckedAt(j)) continue;
      fPoints2->AddLast((points->UncheckedAt(j)));
    }	
  }
  fActivePoints = new Bool_t[fPoints2->GetEntries()];
  for(Int_t k=0;k<fPoints2->GetEntriesFast();k++){
    fActivePoints[k]=kTRUE;
  }
  printf("\n nb hits %d",fPoints2->GetEntries());
  t2=clock();
  //	printf("\nPoints loaded in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
}

//_____________________________________________________________
void AliDisplayFrame::ApplyCuts()
{
  // Applies cuts
  clock_t t1,t2;
  t1=clock();
  
  Float_t		*pxyz;
  Float_t		r,theta,eta,cutmin,cutmax,etamin,etamax,pmom,smin,smax;
  Int_t		nbhits=0;
  AliPoints *pm;
  TParticle *particle;
  
  //Get momentum cut
  smin = gAliDisplay2->GetMomentumMin();
  smax = gAliDisplay2->GetMomentumMax();
  cutmin = 2.0*smin;
  if(smax<0.98) 	cutmax = 2.0*smax;
  else 			cutmax = 100000;
  
  //Get rapidity cut
  smax = gAliDisplay2->GetRapidityMax();
  smin = gAliDisplay2->GetRapidityMin();
  //etamin = 1.5*(2*smin-1);
  //etamax = 1.5*(2*smax-1);
  etamin = smin;
  etamax = smax;
  if(smin<-1.46) etamin = -1000;
  if(smax>1.46) etamax = 1000;
  
  
  if(!fPoints2) return;
  for(Int_t i=0;i<fPoints2->GetEntries();i++){
    pm = dynamic_cast<AliPoints*>(fPoints2->UncheckedAt(i));
    if(!pm) {
      fActivePoints[i]=kFALSE;
      continue;
    }
    particle = pm->GetParticle();
    if(!particle) {
      fActivePoints[i]=kFALSE;
      continue;
    }
    pmom = particle->P();
    if(pmom < cutmin) {
      fActivePoints[i]=kFALSE;
      continue;
    }
    if(pmom > cutmax) {
      fActivePoints[i]=kFALSE;
      continue;
    }
    pxyz = pm->GetP();
    r = TMath::Sqrt(pxyz[0]*pxyz[0]+pxyz[1]*pxyz[1]);
    theta = TMath::ATan2(r,TMath::Abs(pxyz[2]));
    if(theta) eta = -TMath::Log(TMath::Abs(TMath::Tan(0.5*theta)));
    else eta = 1e10;
    if(pxyz[2] < 0) eta = -eta;
    if((eta < etamin) || (eta > etamax)) {
      fActivePoints[i]=kFALSE;
      continue;
    }
    fActivePoints[i]=kTRUE;
    //pm->Draw();
    nbhits += pm->GetN();
  }
  gAliDisplay2->SetNbHits(nbhits);
  t2=clock();
  //	printf("\nCuts applied in....%f sec", ((double)t2-t1)/(10000*CLK_TCK));
  gAliDisplay2->SetNbParticles(GetNbActivePoints()); 
}

//_____________________________________________________________
Int_t AliDisplayFrame::GetNbActivePoints() const
{
  // Returns the number of active points
  Int_t ans=0;
  for(Int_t i=0;i<fPoints2->GetEntries();i++){
    if(fActivePoints[i]) ans++;
  }
  return ans;
}
//_____________________________________________________________
void AliDisplayFrame::DrawX3d()
{
  // Draws using X3d
  TPad *pad = dynamic_cast<TPad*>(gPad);
  pad->cd();
  TView *view = pad->GetView();
  if(!view) return;
  pad->x3d();
}

//_____________________________________________________________
void AliDisplayFrame::SavePadGIF(const char *file)
{
  // Save the current pad in a GIF file
  if(!gPad){
    printf("\nThere is no active pad");
    return;
  }
  gPad->SaveAs(file);
}

//_____________________________________________________________
void AliDisplayFrame::DrawGL()
{
  // Draws using GL
  TPad *pad = dynamic_cast<TPad*>(gPad);
  pad->cd();
  TView *view = pad->GetView();
  if(!view) return;
  pad->x3d("OPENGL");
}

//_____________________________________________________________
void AliDisplayFrame::ExecuteEvent(Int_t event, Int_t px,Int_t py,TObject *)
{
  static Float_t x0,y0,x1,y1;
  static Int_t pxold,pyold;
  static Int_t px0,py0;
  static Int_t linedrawn;
  Float_t temp;
  
  
  switch(event){
  case kMouseMotion:{
    
    AliPoints *p=dynamic_cast<AliPoints*> (gPad->GetSelected());
    if(p){
      gAliDisplay2->SetStatusBar(p->GetName(),1);
      gAliDisplay2->SetStatusBar(p->GetDetector()->GetName(),2);
    }
  }
    break;
  default:break;
  }	
  
  if((!gAliDisplay2->GetZoomMode())&&(gPad->GetView())){
    gPad->GetView()->ExecuteRotateView(event,px,py);
    return;
  }
  
  
  
  if(gAliDisplay2->GetZoomMode()==kTRUE){
    switch(event){
      
    case kButton1Down:{
      gVirtualX->SetLineColor(-1);
      gPad->TAttLine::Modify();
      x0 = gPad->AbsPixeltoX(px);
      y0 = gPad->AbsPixeltoY(py);
      px0 = px;
      py0 = py;
      pxold = px;
      pyold = py;
      linedrawn = 0;
    }
      break;
    case kButton1Motion:{
      if(linedrawn) gVirtualX->DrawBox(px0,py0,pxold,pyold,TVirtualX::kHollow);
      pxold = px;
      pyold = py;
      linedrawn = 1;
      gVirtualX->DrawBox(px0,py0,pxold,pyold,TVirtualX::kHollow);
    }
      break;
      
    case kButton1Up:{
      gPad->GetCanvas()->FeedbackMode(kFALSE);
      if(px == px0) break;
      if(py == py0) break;
      x1 = gPad->AbsPixeltoX(px);
      y1 = gPad->AbsPixeltoY(py);
      if(x1<x0) { 
	temp = x0;
	x0 = x1;
	x1 = temp;
      }
      if(y1<y0) {
	temp = y0;
	y0 = y1;
	y1 = temp;
      }
      printf("\nBox (%f,%f)-(%f,%f)",x0,y0,x1,y1);
      gPad->SetEditable(kTRUE);
      //gPad->Range(x0,y0,x1,y1);
      gPad->SetEditable(kFALSE);
      //gPad->Range(0.5,0.5,1,1);
      //gAliDisplay2->SetZoomFactor(1);
      gPad->Modified(kTRUE);
      gAliDisplay2->Draw();	
      gAliDisplay2->SetZoomMode(kFALSE);
      gPad->SetEditable(kTRUE);
    }
      break; 
    default: break;		
    }		
  }
}
